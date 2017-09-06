#include <math.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>

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
  FILE *fs;
  FILE *logfs=stderr;

  dim_ta nvar;			//number of variables
  cls_ta ncls;		//number of classes
  nel_ta ntest;			//number of test data points

  cls_ta *clist;		//list of class labels

  real_a **brd;			//class border samples
  real_a **grd;			//gradient at the class border
  real_a **test;			//test data vectors
  real_a *result;		//results of prediction
  cls_ta *cls;			//class results 
  real_a *con=NULL;			//estimated confidence
  cls_ta *cls0=NULL;			//classes read from the test file

  multiclass_c<real_a, cls_ta> *classifier;
  real_a *abscissa;
  
  int errcode;

  agf_command_opts opt_args;

  //normalization data:
  real_a **mat;		//transformation matrix
  real_a *ave;		//constant term
  dim_ta nvar1, nvar2;

  opt_args.Qtype=0;
  errcode=agf_parse_command_opts(argc, argv, "a:E:O:Q:w:CeHKMnu", &opt_args);
  if (opt_args.nt==NT_DEFAULT) opt_args.nt=2;
  if (opt_args.W1==-1) opt_args.W1=1.;
  if (errcode==FATAL_COMMAND_OPTION_PARSE_ERROR) return errcode;

  //parse the command line arguments:
  if (argc != 3) {
    FILE *helpfs=stdout;
    fprintf(helpfs, "\n");
    fprintf(helpfs, "Syntax:   classify_c [-Q type] [-a normfile [-u]] \\\n");
    fprintf(helpfs, "                       [-O command [-K] [-M [-E missing] | [-H] [-C]]] \\\n");
    fprintf(helpfs, "                       control test output\n");
    fprintf(helpfs, "\n");
    fprintf(helpfs, "where:\n");
    fprintf(helpfs, "  control     control file\n");
    fprintf(helpfs, "  test        file containing vector data for prediction\n");
    fprintf(helpfs, "  output      files containing the results of the prediction:\n");
    fprintf(helpfs, "                .dat for predictions, .err for error tolerances\n");
    fprintf(helpfs, "\n");
    fprintf(helpfs, "options:\n");
    //printf("  -n          option to normalise the data\n");
    fprintf(helpfs, "  -a normfile file containing normalization data (no default--\n");
    fprintf(helpfs, "                  must always be specified explicitly)\n");
    fprintf(helpfs, "  -u          model (borders) data is not normalized\n");
    fprintf(helpfs, "  -Q          how to solve for probabilities in non-hierarchical scheme: [0]\n");
    fprintf(helpfs, "                0 = constrained matrix inverse [default]\n");
    fprintf(helpfs, "                1 = matrix inversion, no re-normalization\n");
    fprintf(helpfs, "                2 = voting from probabilities, no re-normalization\n");
    fprintf(helpfs, "                3 = voting from classes, no re-normalization\n");
    fprintf(helpfs, "                4 = voting overrides matrix inversion, conditional probabilities\n");
    fprintf(helpfs, "                    are adjusted and re-normalized\n");
    fprintf(helpfs, "  -O command  external command for estimating classes and probabilities\n");
    fprintf(helpfs, "\n");
    fprintf(helpfs, "The following options are in conjunction with the -O option only:\n");
    fprintf(helpfs, "  -M          use LIBSVM format\n");
    fprintf(helpfs, "  -E missing  (in combination with -M) value for missing data\n");
    fprintf(helpfs, "  -H          no header\n");
    fprintf(helpfs, "  -C          no class data\n");
    fprintf(helpfs, "  -K          keep temporary files\n");
    fprintf(helpfs, "\n");
    fprintf(helpfs, "The syntax of the control file is as follows:\n\n");
    fprintf(helpfs, "  <branch>         ::= <model> \"{\" <branch_list> \"}\" | <CLASS>\n");
    fprintf(helpfs, "  <model>          ::= <FNAME> | <partition_list>\n");
    fprintf(helpfs, "  <branch_list>    ::= <branch> | <branch_list> <branch>\n");
    fprintf(helpfs, "  <partition_list> ::= <partition> | <partition_list> <partition>\n");
    fprintf(helpfs, "  <partition>      ::= <FNAME> <class_list> \" %c \" <class_list> \";\"\n", PARTITION_SYMBOL);
    fprintf(helpfs, "  <class_list>     ::= <CLASS> | <class_list> <CLASS>\n");
    fprintf(helpfs, "  <CLASS>          ::= 0 | 1 | 2 | 3 ... | <ncls-1>\n\n");
    fprintf(helpfs, "where:\n");
    fprintf(helpfs, "  <FNAME>        base-name of of files containing class borders binary\n");
    fprintf(helpfs, "                   classifier (<FNAME>.brd, <FNAME>.bgd)\n");
    fprintf(helpfs, "  <CLASS>        class number from zero (0) to the number of classes less one\n");
    fprintf(helpfs, "  <partition_list> describes a non-hierarchical multi-borders model\n");
    fprintf(helpfs, "\n");

    return INSUFFICIENT_COMMAND_ARGS;
  }

  ran_init();			//random numbers resolve ties

  //read in class borders:
  classifier=new multiclass_c<real_a, cls_ta>(argv[0], opt_args.Qtype);

  //ncls : number of class labels
  ncls=classifier->n_class();

  //classifier->print(stdout);

  testfile=argv[1];

  if (opt_args.multicommand==NULL) {
    outfile=new char[strlen(argv[2])+5];
    strcpy(outfile, argv[2]);
    strcat(outfile, ".dat");

    confile=new char[strlen(argv[2])+5];
    strcpy(confile, argv[2]);
    strcat(confile, ".err");

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

    if (classifier->n_feat() != nvar) {
      fprintf(stderr, "classify_m: Dimensions of classifier (%d) do not match those of test data (%d).\n",
		classifier->n_feat(), nvar);
      exit(DIMENSION_MISMATCH);
    }
  } else {
    if (classifier->n_feat() != nvar && opt_args.multicommand==NULL) {
      fprintf(stderr, "classify_m: Dimension of classifier (%d) do not match dimension of test data (%d).\n",
		classifier->n_feat(), nvar);
      exit(DIMENSION_MISMATCH);
    }
  }

  //begin the classification scheme:
  result=new real_a[ntest];

  if (opt_args.multicommand==NULL) {
    //allocate space for probabilities:
    con=new real_a[ntest];
    for (nel_ta i=0; i<ntest; i++) {

      result[i]=classifier->ret(test[i], con[i]);

      printf("%g +/- %g\n", result[i], con[i]);
    }

  } else {
      fprintf(stderr, "classify_c: Continuum estimates with an external classifier not yet implemented.  Sorry\n");
      exit(PARAMETER_OUT_OF_RANGE);

  }

  //printf("\n");

  if (opt_args.multicommand==NULL) {
    //write the results to a file:
    fs=fopen(outfile, "w");
    if (fs == NULL) {
      fprintf(stderr, "Unable to open file, %s, for writing\n", outfile);
      return UNABLE_TO_OPEN_FILE_FOR_WRITING;
    }
    fwrite(result, sizeof(real_a), ntest, fs);
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
      class_eval(cls0, cls, ntest);
      actab=con_acc_table(cls0, cls, con, ntest, NCONHIST);
      print_con_acc(actab, ncls, NCONHIST);
      delete [] actab[0];
      delete [] actab;
    }
  }

  //clean up:
  delete classifier;
  delete [] clist;

  delete [] result;
  delete [] cls;
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


