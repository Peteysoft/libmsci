#include <stdio.h>
#include <getopt.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <sys/timeb.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "peteys_tmpl_lib.h"
#include "parse_command_opts.h"
#include "randomize.h"
#include "linked.h"

#include "agf_lib.h"

//because we love arbitrary limits:
#define COMMANDBUF_SIZE 1000000
#define END_INDICATOR '\0'

using namespace std;
using namespace libpetey;
using namespace libagf;

int main (int argc, char **argv) {
  char *controlfile;		//name of control file	
  char *trainfile;		//name of file containing training data
  char *modelbase;		//base name for output model files
  char *outfile;		//name of output control file
  char *extra=NULL;		//extra options to pass to train command

  FILE *commandfs;		//write commands to this file stream
  FILE *infs;			//input control file stream
  FILE *outfs;			//output control file stream
  char *buffer=NULL;		//input file buffer
  long bufsize;

  char *train;			//file containing training data
  char *commandname;		//command to train binary classifiers
  char *partcom=NULL;		//command to partition classes
  char *concom=NULL;		//optional file conversion command
  char *precom=NULL;	//optional external binary classifier (predict command)

  cls_ta ncls;

  int32_t nsv;			//number of singular values to keep
  char *command=NULL;		//command for performing the normalization
  char *normfile=NULL;		//normalization file

  //parsing is now done primarily by the data structures themselves:
  multiclass_hier<real_a, cls_ta> *shell;

  //for naming temporary files:
  //timeb date;
  //pid_t pid;
  unsigned long session_id;           //for naming temporary files

  agf_command_opts opt_args;
  void *opt[30];
  int flag[30];
  int testflag;			//test run (simply prints out the commands)

  char *tempbase=NULL;		//base name of normalized training data

  int err;
  char **argv1;
  int argc1;

  //we want to parse out the options that only apply at the top level
  //while keeping the ones to pass to class_borders
  //--the former include options for normalization or pre-processing:
  opt[3]=&nsv;
  argc=parse_command_opts(argc, argv, "a0nS-+^MKOuZAx", "%s%%%d%s%s%s%%%s%%%%", opt, flag, 3);
  if (argc<0) {
    fprintf(stderr, "multi_borders: error parsing command line\n");
    exit(FATAL_COMMAND_OPTION_PARSE_ERROR);
  }


  if (flag[5]) {
    //extra training options:
    extra=(char *) opt[5];
  } else {
    //use null string since it will be added to commandname:
    extra=new char[1];
    extra[0]='\0';
  }

  //if we are using class_borders, check options for correctness:
  if (flag[4]) {
    //get rid of command name at head of argument list:
    argc1=argc-1;
    argv++;
    commandname=(char *) opt[4];
    opt_args.uflag=0;
  } else {
    argc1=argc;
    err=agf_parse_command_opts(argc1, argv, "h:i:I:k:l:N:q:r:s:t:v:V:W:", &opt_args);
    //argc1=parse_command_opts(argc, argv, "hiIklrstvVWu", "%d%d%d%d%g%g%d%g%g%g%g%", opt+6, flag+6, 1);
    if (err==FATAL_COMMAND_OPTION_PARSE_ERROR) {
      fprintf(stderr, "multi_borders: error parsing command line\n");
      exit(err);
    }
    commandname=NULL;
  }

  if (argc1<3) {
    FILE *helpfs=stdout;
    fprintf(helpfs, "\nsyntax: multi_borders [-n] [-a normfile] [-S nsv] [-O pcom] [trainopt] \\\n");
    fprintf(helpfs, "                          [-- command [-M] [-+ extra] [-^ fcom]] \\\n");
    fprintf(helpfs, "                          [control] train border out\n");
    fprintf(helpfs, "\n");
    fprintf(helpfs, "where:\n");
    fprintf(helpfs, "  control        input control file (reads from stdin if omitted)\n");
    fprintf(helpfs, "  train          base filename for binary training data:\n");
    fprintf(helpfs, "                   .vec for vectors;\n");
    fprintf(helpfs, "                   .cls for classes.\n");
    fprintf(helpfs, "                   .std for transformation/normalization matrix\n");
    fprintf(helpfs, "                        (unless -a specified)\n");
    fprintf(helpfs, "  border         base name for files containing border data\n");
    fprintf(helpfs, "                   (file names are recorded in output control file)\n");
    fprintf(helpfs, "  out            output control file to feed to classify_m\n\n");
    fprintf(helpfs, "\noptions:\n");
    fprintf(helpfs, "  trainopt       are a series of options to pass to class_borders\n");
    fprintf(helpfs, "                     (-k, -i, -I, -N, -r, -s, -t, -v, -V, -W  --for an\n");
    fprintf(helpfs, "                     explanation of each, run class_borders with no arguments)\n");
    fprintf(helpfs, "\n");
    fprintf(helpfs, "  -n             normalize training data with standard deviations\n");
    fprintf(helpfs, "  -S nsv         perform SVD, use nsv singular values\n");
    fprintf(helpfs, "  -u             store borders data in un-normalized coordinates\n");
    fprintf(helpfs, "  -a normfile    normalize/use normfile for normalization data\n");
    fprintf(helpfs, "                     (default is <train>.std\n");
    fprintf(helpfs, "\n");
    fprintf(helpfs, "  -A             training data is in ASCII (LVQ) format\n");
    fprintf(helpfs, "  -M             training data is in LIBSVM format (only with -- or -^)\n");
    fprintf(helpfs, "\n");
    fprintf(helpfs, "  -- command     command to use for training [%s%s%s] (data in ASCII format)\n", AGF_COMMAND_PREFIX, 
			AGF_BINARY_CLASSIFIER, AGF_OPT_VER);
    fprintf(helpfs, "  -+ extra       additional options to pass to training command\n");
    fprintf(helpfs, "                   if -- specified they are defaults, not extras\n");
    fprintf(helpfs, "  -^ fcom        command to convert data to appropriate format for command\n");
    fprintf(helpfs, "\n");
    //fprintf(helpfs, "  -0             trial run only: print out commands but do not execute them\n");
    fprintf(helpfs, "  -O pcom        \"accelerator\" mode: convert existing models to border samples\n");
    fprintf(helpfs, "                   pcom is an external command that returns probability estimates\n");
    fprintf(helpfs, "  -Z             accelerate LIBSVM model using in-house codes\n");
    fprintf(helpfs, "\n");
    fprintf(helpfs, "  -K             keep temporary files; commands to train the model are written\n");
    fprintf(helpfs, "                   to stdout\n");
    fprintf(helpfs, "  -x             run in the background\n");
    fprintf(helpfs, "\n");
    fprintf(helpfs, "The syntax of the control file is as follows:\n\n");
    fprintf(helpfs, "  <branch>         ::= <model> \"{\" <branch_list> \"}\" | <CLASS>\n");
    fprintf(helpfs, "  <model>          ::= <OPTIONS> [<CODE>] | <partition_list>\n");
    fprintf(helpfs, "  <branch_list>    ::= <branch> | <branch_list> <branch>\n");
    fprintf(helpfs, "  <partition_list> ::= <partition> | <partition_list> <partition>\n");
    fprintf(helpfs, "  <partition>      ::= <OPTIONS> <class_list> \" %c \" <class_list> \";\"\n", PARTITION_SYMBOL);
    fprintf(helpfs, "  <class_list>     ::= <CLASS> | <class_list> <CLASS>\n\n");
    fprintf(helpfs, "  <OPTIONS>        ::= \"\"\"\" | \"\".\"\" | \"\" \"\"\n");
    fprintf(helpfs, "                      | \"\"\"[\"-v\" <vmin>] [\"-V\" <vmax>] [\"-W\" <W>]\n");
    fprintf(helpfs, "                        [\"-s\" <n>] [\"-t\" <tol>] ...\"\"\"\n");
    fprintf(helpfs, "  <CLASS>          ::= 0 | 1 | 2 | 3 ... | <ncls-1>\n\n");
    fprintf(helpfs, "where:\n");
    fprintf(helpfs, "  <OPTIONS>      are a double-quoted list of options to pass to class borders\n");
    fprintf(helpfs, "                 - a null string (\"\") means take options from the command line\n");
    fprintf(helpfs, "                 - a period in quotes (\".\") means use the previous option list\n");
    fprintf(helpfs, "                   on the same level or lower\n");
    fprintf(helpfs, "                 - to make the options list empty, use a whitespace character\n");
    fprintf(helpfs, "                   (i.e., \" \")\n");
    fprintf(helpfs, "  <CLASS>        class number from zero (0) to the number of classes less one\n");
    fprintf(helpfs, "  <partition_list> describes a non-hierarchical multi-borders model\n");
    fprintf(helpfs, "  <CODE>         is the letter code for one of the \"direct\" classifiers:\n");
    fprintf(helpfs, "                 - A for direct AGF\n");
    fprintf(helpfs, "                 - K for KNN\n");
    fprintf(helpfs, "                 - G for an external classifier whose command is contained in\n");
    fprintf(helpfs, "                   <OPTIONS>\n");
    fprintf(helpfs, "\n");
    exit(0);
  }

  if (argc1==3) {
    controlfile=NULL;
    trainfile=argv[0];
    modelbase=argv[1];
    outfile=argv[2];
  } else {
    controlfile=argv[0];
    trainfile=argv[1];
    modelbase=argv[2];
    outfile=argv[3];
  }

  if (flag[0]) {
    normfile=new char[strlen((char *) opt[0])+1];
    strcpy(normfile, (char *) opt[0]);
  }
  //"un-normalize" flag:
  opt_args.uflag=flag[10];
  opt_args.Kflag=flag[8];			//keep temporary files
  opt_args.xflag=flag[13];

  //write statements to stdout:
  char *commandbuf;
  if (opt_args.Kflag) {
    commandfs=stdout;
  } else {
    commandbuf=new char[COMMANDBUF_SIZE];
    commandfs=fmemopen(commandbuf, COMMANDBUF_SIZE, "w");
  }

  //create a unique session id (or at least as close as we can come):
  session_id=seed_from_clock();

  //not a great idea--what if we want to throw the control file in there 
  //as well?
/*  //check for non-existent folders in output name:
  for (int i=0; modelbase[i]!='\0'; i++) {
    struct stat fstatus;
    if (modelbase[i]=='/') {
      modelbase[i]='\0';
      stat(modelbase, &fstatus);
      if (S_ISDIR(fstatus.st_mode)==0) {
        fprintf(commandfs, "mkdir %s\n", modelbase);
      }
      modelbase[i]='/';
    }
  }*/

  //pre-processing options:
  if (flag[0] || flag[2] || flag[3]) {
    //stick transformed training data in a temporary file:
    tempbase=new char [strlen(modelbase)+20];
    sprintf(tempbase, "%s.%u", modelbase, session_id);

    if (normfile==NULL) {
      normfile=new char [strlen(modelbase)+5];
      sprintf(normfile, "%s.std", modelbase);
    }
    command=new char [strlen(trainfile)+strlen(normfile)+strlen(tempbase)+
		strlen(AGF_COMMAND_PREFIX)+strlen(AGF_LTRAN_COM)+
		strlen(AGF_OPT_VER)+104];
    sprintf(command, "%s%s%s -a %s", AGF_COMMAND_PREFIX, AGF_LTRAN_COM,
		AGF_OPT_VER, normfile);
    if (flag[2]) {
      sprintf(command+strlen(command), " -n");
    }
    if (flag[3]) {
      sprintf(command+strlen(command), " -S %d", nsv);
    }

    //class_borders can now work with both native binary and ASCII formats:
    if (flag[4] || flag[12]) {
      char mc[3];
      mc[0]='\0';
      if (flag[7]) strcpy(mc, "-M");
      sprintf(command+strlen(command), " -A %s %s %s", mc, trainfile, tempbase);
      fprintf(commandfs, "%s\n", command);
      train=new char [strlen(tempbase)+1];
      sprintf(train, "%s", tempbase);
    } else {
      sprintf(command+strlen(command), " %s.vec %s.vec", trainfile, tempbase);
      fprintf(commandfs, "%s\n", command);
      if (opt_args.uflag) {
        train=new char [strlen(trainfile)+1];
        sprintf(train, "%s", trainfile);
      } else {
        train=new char [strlen(tempbase)+1];
        sprintf(train, "%s", tempbase);
	if (flag[12]==0) {
          fprintf(commandfs, "cp %s.cls %s.cls\n", trainfile, tempbase);
	}
      }
    }
  } else {
    train=new char [strlen(trainfile)+1];
    sprintf(train, "%s", trainfile);
  }

  //fucking disaster...
  if (flag[4]) {
    //we can pass options by tagging them to the end:
    argc=argc1-4;
    argv+=4;
    //for (int i=0; i<argc; i++) printf("%s\n", argv[i]);
  } else {
    //default options are whatever was parsed by agf_parse_command_opts:
    argv-=argc-argc1-1;
    argc-=argc1+1;
  }

  //-+ option is used to pass extra command options
  //if the training command has been explicity specified, they are default
  //options, otherwise they are extra options and we tack them onto the 
  //command name:
  if (flag[5] && flag[4]) {
    argc++;
    argv--;
    //"extra" parameters passed with -+ option:
    argv[0]=extra;
  }

  //optional parameters:
  if (flag[4]) {
    char mc[3];
    mc[0]='\0';
    if (flag[7]) strcpy(mc, "-M");		//duplicate code...
    partcom=new char[strlen(AGF_COMMAND_PREFIX)+strlen(AGF_OPT_VER)+100];
    sprintf(partcom, "%s%s%s -A %s", AGF_COMMAND_PREFIX, 
		AGF_PARTCOM, AGF_OPT_VER, mc);
    if (flag[6]) {
      concom=(char *) opt[6];
      strcat(partcom, " -1");
    }
  }

  //use this command for training
  if (commandname==NULL) {
    if (opt_args.uflag && normfile!=NULL) {
      //if we want to un-normalize the data, we have to pass the option up to 
      //class borders through the command name:
      commandname=new char[strlen(AGF_COMMAND_PREFIX)+strlen(normfile)+
		strlen(AGF_BINARY_CLASSIFIER)+strlen(AGF_OPT_VER)+
		strlen(extra)+200];
      sprintf(commandname, "%s%s%s %s -u -a %s", AGF_COMMAND_PREFIX, 
		AGF_BINARY_CLASSIFIER, AGF_OPT_VER, extra, normfile);
    } else {
      commandname=new char[strlen(AGF_COMMAND_PREFIX)+strlen(extra)+
		strlen(AGF_BINARY_CLASSIFIER)+strlen(AGF_OPT_VER)+200];
      sprintf(commandname, "%s%s%s %s", AGF_COMMAND_PREFIX, 
		AGF_BINARY_CLASSIFIER, AGF_OPT_VER, extra);
    }
    //class_borders can work with both native binary and ASCII files:
    if (flag[12]) strcat(commandname, " -A");
    if (flag[7]) strcat(commandname, " -M");
  }

  //input control file:
  if (controlfile!=NULL) {
    infs=fopen(controlfile, "r");
    if (infs==NULL) {
      fprintf(stderr, "multi_borders: unable to open file %s for reading\n",
		controlfile);
      exit(UNABLE_TO_OPEN_FILE_FOR_READING);
    }
  } else {
    //because the parser has so much look-ahead, we need to buffer the input:
    linked_list<char> chars;
    int c1;
    c1=getchar();
    while (c1!=EOF) {
      chars.add((char)c1);
      c1=getchar();
    }
    chars.add('\0');
    buffer=chars.make_array(bufsize);
    infs=fmemopen(buffer, bufsize, "w+");
  }

  //output control file (last argument):
  outfs=fopen(outfile, "w");
  if (outfs==NULL) {
    fprintf(stderr, "multi_borders: unable to open file %s for writing\n",
		outfile);
    exit(UNABLE_TO_OPEN_FILE_FOR_WRITING);
  }

  //done with set-up phase, now we can actually get to the real work!
  //"accelerator" function:
  if (flag[9] || flag[11]) {
    if (flag[9]) {
      precom=(char *) opt[9];
    } else {
      precom=new char[1];
      precom[0]='\0';
    }
    if (flag[4]==0) {
      //pass -Z option to class_borders:
      if (flag[11]) {
        strcat(commandname, " -Z ");
	//dammit:
	for (int i=0; i<argc; i++) {
          strcat(commandname, argv[i]);
          strcat(commandname, " ");
	}
      }
    }
    shell=new multiclass_hier<real_a, cls_ta>(infs, 0, precom, 
		    flag[7], opt_args.Kflag);
  } else {
    //parse the control file and build the data structure:
    shell=new multiclass_hier<real_a, cls_ta>(infs, argc, argv);
  }

  //print out the commands to generate the model:
  ncls=shell->generate_commands(commandfs, train, modelbase,
		commandname, partcom, concom, opt_args.Kflag, session_id);
  fclose(infs);

  //delete temporary files:
  if (normfile!=NULL && opt_args.Kflag==0) {
    if (flag[6]) {
      fprintf(commandfs, "rm %s\n", tempbase);
    } else {
      fprintf(commandfs, "rm %s.vec\n", tempbase);
      if (opt_args.uflag==0) {
        fprintf(commandfs, "rm %s.cls\n", tempbase);
      }
    }
    delete [] normfile;
  }

  //now actually do all the work:
  if (opt_args.Kflag!=1) {
    char *c;
    int err;
    fprintf(commandfs, "%c", END_INDICATOR);
    fclose(commandfs);
    do {
      for (c=commandbuf; *c!='\n'; c++);
      *c='\0';
      //run in the background:
      if (opt_args.xflag) {
        char *com=new char[strlen(commandbuf)+2];
        sprintf(com, "%s&", commandbuf);
        printf("%s\n", com);
        err=system(com);
        delete [] com;
      } else {
        printf("%s\n", commandbuf);
        err=system(commandbuf);
      }
      if (err!=0) exit(err);
      commandbuf=c+1;
    } while (*commandbuf!=END_INDICATOR);
  }

  //print out the new control file:
  shell->print(outfs, modelbase);
  fclose(outfs);

  //clean up:
  delete shell;

  if (command!=NULL) delete [] command;
  if (tempbase!=NULL) delete [] tempbase;
  delete [] train;

  delete [] commandname;
  delete [] partcom;

  delete [] extra;

  if (buffer != NULL) delete [] buffer;

}

