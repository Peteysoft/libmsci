#include <math.h>
#include <string.h>
#include <stdio.h>

#include "full_util.h"
#include "agf_lib.h"
#include "randomize.h"

using namespace std;
using namespace libagf;
using namespace libpetey;

int main(int argc, char *argv[]) {
  char *testfile;		//test data
  char *outfile;		//output classes
  char *confile;		//output confidences
  FILE *fs;			//general file stream
  FILE *logfs=stderr;		//for logging error messages

  dim_ta nvar;			//number of variables
  cls_ta ncls;			//number of classes
  nel_ta ntest;			//number of test data points

  cls_ta *clist;		//list of class labels

  real_a **brd;			//class border samples
  real_a **grd;			//gradient at the class border
  real_a **test;		//test data vectors
  cls_ta *result;		//results of classification
  real_a *con=NULL;		//estimated confidence
  cls_ta *cls0=NULL;		//classes read from the test file

  multiclass_hier<real_a, cls_ta> *classifier;	//multi-class classifier
  
  int errcode;			//returned error code

  agf_command_opts opt_args;	//command line options

  //normalization data:
  real_a **mat;		//transformation matrix
  real_a *ave;		//constant term
  dim_ta nvar1, nvar2;

  opt_args.Qtype=-1;
  opt_args.algtype=0;
  errcode=agf_parse_command_opts(argc, argv, "O:a:Q:w:c:nuMCHE:Ky:Z", &opt_args);
  if (opt_args.nt==NT_DEFAULT) opt_args.nt=2;
  if (errcode==FATAL_COMMAND_OPTION_PARSE_ERROR) return errcode;

  //parse the command line arguments:
  if (argc != 3) {
    FILE *helpfs=stdout;
    fprintf(helpfs, "\n");
    fprintf(helpfs, "purpose:  Performs statistical classification using a multi-borders model\n");
    fprintf(helpfs, "\n");
    fprintf(helpfs, "Syntax:   classify_m [-Q type] [-c funcode] [-a normfile [-u]] \\\n");
    fprintf(helpfs, "                       [-O command [-K] [-M [-E missing] | [-H] [-C]]] \\\n");
    fprintf(helpfs, "                       control test output\n");
    fprintf(helpfs, "\n");
    fprintf(helpfs, "where:\n");
    fprintf(helpfs, "  control     control file or all-in-one model\n");
    fprintf(helpfs, "  test        file containing vector data to be classified\n");
    fprintf(helpfs, "  output      files containing the results of the classification:\n");
    fprintf(helpfs, "                .cls for classes, .con for confidence ratings\n");
    fprintf(helpfs, "\n");
    fprintf(helpfs, "options:\n");
    //printf("  -n          option to normalise the data\n");
    fprintf(helpfs, "  -a normfile file containing normalization data (no default--\n");
    fprintf(helpfs, "                  must always be specified explicitly)\n");
    fprintf(helpfs, "  -u          model (borders) data is not normalized\n");
    fprintf(helpfs, "  -Q          how to solve for probabilities in non-hierarchical scheme: [0]\n");
    fprintf(helpfs, "                0 = constrained inverse [default]\n");
    fprintf(helpfs, "                      (may not produce an optimal estimate in all cases)\n");
    fprintf(helpfs, "                1 = linear least squares, no re-normalization\n");
    fprintf(helpfs, "                2 = voting from pdf, no re-normalization\n");
    fprintf(helpfs, "                3 = voting from class label, no re-normalization\n");
    fprintf(helpfs, "                4 = voting from pdf overrides least squares, conditional\n");
    fprintf(helpfs, "                    probabilities are adjusted and re-normalized\n");
    fprintf(helpfs, "                5 = voting from class overrides least squares, conditional\n");
    fprintf(helpfs, "                    probabilities are adjusted and re-normalized\n");
    fprintf(helpfs, "                6 = experimental\n");
    fprintf(helpfs, "                7 = constrained inverse\n");
    fprintf(helpfs, "                      (may be extremely inefficient (NP) for some cases)\n");
    fprintf(helpfs, "                8 = voting from pdf, corrected and normalized\n");
    fprintf(helpfs, "                      (designed for orthogonal coding matrices)\n");
    fprintf(helpfs, "  -C          no class data\n");
    fprintf(helpfs, "  -E missing  (in combination with -M) value for missing data\n");
    fprintf(helpfs, "  -H          no header\n");
    fprintf(helpfs, "  -K          keep temporary files\n");
    fprintf(helpfs, "  -M          use LIBSVM format\n");
    fprintf(helpfs, "  -O command  external command for estimating classes and probabilities\n");
    fprintf(helpfs, "  -y path     path to data files\n");
    fprintf(helpfs, "  -Z          use in-house SVM codes\n");
    fprintf(helpfs, "  -c funcode  sigmoid function for transforming decision values:\n");
    fprintf(helpfs, "  	            0 = tanh\n");
    fprintf(helpfs, "  	            1 = erf\n");
    fprintf(helpfs, "  	            2 = logistic function [f(x)=2/(1+exp(x))]\n");
    fprintf(helpfs, "\n");
    fprintf(helpfs, "The syntax of the control file is as follows:\n");
    fprintf(helpfs, "\n");
    fprintf(helpfs, "  <branch>         ::= <model> \"{\" <branch_list> \"}\" | <CLASS>\n");
    fprintf(helpfs, "  <model>          ::= <FNAME> [<CODE> <OPTIONS>] | <partition_list>\n");
    fprintf(helpfs, "  <branch_list>    ::= <branch> | <branch_list> <branch>\n");
    fprintf(helpfs, "  <partition_list> ::= <partition> | <partition_list> <partition>\n");
    fprintf(helpfs, "  <partition>      ::= <FNAME> <class_list> \" %c \" <class_list> \";\"\n", PARTITION_SYMBOL);
    fprintf(helpfs, "  <class_list>     ::= <CLASS> | <class_list> <CLASS>\n");
    fprintf(helpfs, "  <CLASS>          ::= 0 | 1 | 2 | 3 ... | <ncls-1>\n\n");
    fprintf(helpfs, "where:\n");
    fprintf(helpfs, "  <FNAME>        base-name of files containing class borders binary classifier\n");
    fprintf(helpfs, "                   (<FNAME>.brd, <FNAME>.bgd), or a model data for an external\n");
    fprintf(helpfs, "                   binary classifier or data for direct classifier\n");
    fprintf(helpfs, "                    -- name can contain any character except spaces\n");
    fprintf(helpfs, "                   but must start with a letter\n");
    fprintf(helpfs, "  <CLASS>        class number from zero (0) to the number of classes less one\n");
    fprintf(helpfs, "  <partition_list> describes a non-hierarchical multi-borders model\n");
    fprintf(helpfs, "  <CODE>         is the letter code for one of the \"direct\" classifiers:\n");
    fprintf(helpfs, "                 - A for direct AGF\n");
    fprintf(helpfs, "                 - K for KNN\n");
    fprintf(helpfs, "                 - G for an external classifier whose command is contained in\n");
    fprintf(helpfs, "                   <OPTIONS>\n");
    fprintf(helpfs, "  <OPTIONS>      is a quoted list of options for a direct classifier\n");
    fprintf(helpfs, "\n");

    return INSUFFICIENT_COMMAND_ARGS;
  }

  ran_init();			//random numbers resolve ties

  fs=fopen(argv[0], "r");
  if (fs==NULL) {
    fprintf(stderr, "classify_m: unable to open file, %s, for reading\n", argv[0]);
    exit(UNABLE_TO_OPEN_FILE_FOR_READING);
  }
  classifier=new multiclass_hier<real_a, cls_ta>();
  //check first for all-in-one ASCII model file:
  if (classifier->load(fs, opt_args.Qtype)==PARAMETER_OUT_OF_RANGE) {
    delete classifier;
    rewind(fs);
    //otherwise use control file:
    classifier=new multiclass_hier<real_a, cls_ta>(fs,		//file stream  
		opt_args.Qtype,			//method to solve for prob.
		opt_args.path, 			//path to dat files
		opt_args.multicommand, 		//external binary classifier
		opt_args.Mflag,			//external uses LIBSVM format
		opt_args.Kflag, 		//keep temporary files
		opt_args.algtype, 		//sig. fn extrap. bin prob
		opt_args.Zflag);		//in-house SVM codes
  }
  fclose(fs);

  //classifier->print(stdout);

  testfile=argv[1];

  if (opt_args.multicommand==NULL) {
    outfile=new char[strlen(argv[2])+5];
    strcpy(outfile, argv[2]);
    strcat(outfile, ".cls");

    confile=new char[strlen(argv[2])+5];
    strcpy(confile, argv[2]);
    strcat(confile, ".con");

    test=read_vecfile<real_a>(testfile, ntest, nvar);
    if (nvar == -1 || ntest==-1) {
      fprintf(stderr, "Error reading input file: %s\n", testfile);
      exit(FILE_READ_ERROR);
    }
    if (test == NULL) {
      fprintf(stderr, "Unable to open file for reading: %s\n", testfile);
      exit(UNABLE_TO_OPEN_FILE_FOR_WRITING);
    }
  } else {
    outfile=new char [strlen(argv[2])+1];
    strcpy(outfile, argv[2]);
    if (opt_args.Mflag) {
      ntest=read_svm(testfile, test, cls0, nvar, opt_args.missing, opt_args.Uflag);
      if (nvar == -1 || ntest==-1) {
        fprintf(stderr, "Error reading input file: %s\n", testfile);
        exit(FILE_READ_ERROR);
      }
    } else {
      ntest=read_lvq(testfile, test, cls0, nvar, opt_args.Hflag+2*opt_args.Cflag);
      if (nvar == -1 || ntest==-1) {
        fprintf(stderr, "Error reading input file: %s\n", testfile);
        exit(FILE_READ_ERROR);
      }
    }
  }

  fprintf(logfs, "%d test vectors found in file %s\n", ntest, testfile);

  //normalization:
  if ((opt_args.uflag || opt_args.normflag) && opt_args.normfile==NULL) {
    fprintf(stderr, "classify_m: please specify normalization file with -a\n");
    exit(FATAL_COMMAND_OPTION_PARSE_ERROR);
    //opt_args.normfile=new char [strlen(argv[0])+5];
    //sprintf(opt_args.normfile, "%s.std", argv[0]);
  }

  if (opt_args.normfile!=NULL) {
    mat=read_stats2(opt_args.normfile, ave, nvar1, nvar2);

    errcode=classifier->ltran(mat, ave, nvar1, nvar2, opt_args.uflag);
    if (errcode!=0) exit(errcode);

    if (classifier->n_feat_t() != nvar) {
      fprintf(stderr, "classify_m: Dimensions of classifier (%d) do not match those of test data (%d).\n",
		classifier->n_feat_t(), nvar);
      exit(DIMENSION_MISMATCH);
    }
  } else {
    if (classifier->n_feat_t() != nvar && opt_args.multicommand==NULL) {
      fprintf(stderr, "classify_m: Dimension of classifier (%d) do not match dimension of test data (%d).\n",
		classifier->n_feat_t(), nvar);
      exit(DIMENSION_MISMATCH);
    }
  }

  //begin the classification scheme:
  result=new cls_ta[ntest];

  //ncls : number of class labels
  ncls=classifier->n_class();

  if (classifier->max_depth() == 1 && ncls>2) {
    cls_ta ncls1=0;			//largest class label-1
    //in case there are missing classes:
    //(ncls1 : maximum value for class label plus 1)
    clist=new cls_ta[ncls];
    classifier->class_list(clist);
    for (cls_ta i=0; i<ncls; i++) if (clist[i]>=ncls1) ncls1=clist[i]+1;

    if (opt_args.multicommand==NULL) {
      real_a *pdf=NULL;			//conditional probabilities (if applicable)
      real_a *pdf1=NULL;		//pdfs before mapping

      //allocate space for probabilities:
      con=new real_a[ntest];
      pdf1=new real_a[ncls];
      pdf=new real_a[ncls1];
      for (cls_ta i=0; i<ncls1; i++) pdf[i]=0;

      //perform the classifications:
      for (nel_ta i=0; i<ntest; i++) {
        result[i]=classifier->classify_t(test[i], pdf1);
        for (cls_ta j=0; j<ncls; j++) pdf[clist[j]]=pdf1[j];
        //result[i]=clist[result[i]];

        con[i]=(ncls*pdf[result[i]]-1)/(ncls-1);

        //print results to standard out:
        for (cls_ta j=0; j<ncls; j++) printf(" %9.6f", pdf[j]);
        printf("% 4d", result[i]);
        printf("\n");
      }
      delete [] pdf;
      delete [] pdf1;
    } else {
      real_a **p;

      //allocate space for probabilities:
      p=new real_a *[ntest];
      p[0]=new real_a[ntest*ncls];
      for (nel_ta i=1; i<ntest; i++) p[i]=p[0]+i*ncls; 

      //perform the classifications:
      classifier->batch_classify_t(test, result, p, ntest, nvar);

      //output results in same format as LIBSVM:
      fs=fopen(outfile, "w");
      fprintf(fs, "labels");
      for (cls_ta j=0; j<ncls; j++) fprintf(fs, " %d", clist[j]);
      fprintf(fs, "\n");
      for (nel_ta i=0; i<ntest; i++) {
        fprintf(fs, "%d", result[i]);
        for (cls_ta j=0; j<ncls; j++) fprintf(fs, " %9.6f", p[i][j]);
        fprintf(fs, "\n");
      }
      fclose(fs);

      //calculate "confidence ratings" for class evaluation:
      //inverse mapping from class label to class number:
      cls_ta *clisti=new cls_ta[ncls1];
      con=new real_a[ntest];
      for (cls_ta j=0; j<ncls; j++) {
        clisti[clist[j]]=j;
      }
      for (nel_ta i=0; i<ntest; i++) {
        con[i]=(ncls*p[i][clisti[result[i]]]-1)/(ncls-1);
      }
      delete [] clisti;

      delete [] p[0];
      delete [] p;
    }
    delete [] clist;
  } else {
    con=new real_a[ntest];
    if (opt_args.multicommand==NULL) {
      for (nel_ta i=0; i<ntest; i++) {
        result[i]=classifier->classify_t(test[i], con[i]);
        con[i]=(ncls*con[i]-1)/(ncls-1);
      }
    } else {
      classifier->batch_classify_t(test, result, con, ntest, nvar);
      fs=fopen(outfile, "w");
      for (nel_ta i=0; i<ntest; i++) {
        fprintf(fs, "%d %g\n", result[i], con[i]);
      }
      fclose(fs);
      //calculate confidence ratings for class evaluations:
      if (cls0!=NULL && opt_args.Cflag==0) {
        for (nel_ta i=0; i<ntest; i++) con[i]=(ncls*con[i]-1)/(ncls-1);
      }
    }
  }

  //printf("\n");

  if (opt_args.multicommand==NULL) {
    //write the results to a file:
    fs=fopen(outfile, "w");
    if (fs == NULL) {
      fprintf(stderr, "Unable to open file, %s, for writing\n", outfile);
      return UNABLE_TO_OPEN_FILE_FOR_WRITING;
    }
    fwrite(result, sizeof(cls_ta), ntest, fs);
    fclose(fs);

    //write the results to a file:
    fs=fopen(confile, "w");
    if (fs == NULL) {
      fprintf(stderr, "Unable to open file, %s, for writing\n", confile);
      return UNABLE_TO_OPEN_FILE_FOR_WRITING;
    }
    fwrite(con, sizeof(real_a), ntest, fs);
    fclose(fs);

    delete [] confile;
  } else {
    if (cls0!=NULL && opt_args.Cflag==0) {
      real_a **actab;
      class_eval(cls0, result, ntest);
      actab=con_acc_table(cls0, result, con, ntest, NCONHIST);
      print_con_acc(actab, ncls, NCONHIST);
      delete [] actab[0];
      delete [] actab;
    }
  }

  //clean up:
  delete classifier;

  delete [] result;
  if (con!=NULL) delete [] con;
  if (cls0!=NULL) delete [] cls0;
  delete [] test[0];
  delete [] test;
  delete [] outfile;
  if (opt_args.normfile!=NULL) {
    delete [] opt_args.normfile;
    delete_matrix(mat);
    delete [] ave;
  }

  ran_end();

  return 0;

}


