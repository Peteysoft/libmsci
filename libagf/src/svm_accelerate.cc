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

//do calculations in double-precision:
typedef float calc_t;

int main(int argc, char *argv[]) {
  FILE *fs;

  //do all the calculations in double-precision:
  svm_multi<calc_t, cls_ta> *classifier;
  borders1v1<calc_t, cls_ta> *accel;

  //transformation matrix:
  real_a **mat=0;
  real_a *ave;
  dim_ta nvar1, nvar2;

  dim_ta nvar;			//number of variables
  cls_ta nclass;		//number of classes
  nel_ta ntest;			//number of test data points

  real_a **test;		//test data vectors
  cls_ta *cls;

  int err;

  //do calculations in double-precision:
  calc_t **test2;		//test data vectors
  calc_t **mat2;
  calc_t *b2;

  agf_command_opts opt_args;

  opt_args.algtype=0;
  //err=agf_parse_command_opts(argc, argv, "a:i:ns:t:uAEMN:U", &opt_args);
  err=agf_parse_command_opts(argc, argv, "i:s:t:AEMN:U", &opt_args);
  if (err==FATAL_COMMAND_OPTION_PARSE_ERROR) return err;

  //parse the command line arguments:
  if (argc != 3) {
    printf("Syntax:   svm_accelerate \\\n");
    printf("                  [-A [-M [-E missing]]] \\\n");
    //printf("                  [-n] [-u] [-a normfile] \\\n");
    printf("                  modelfile train output\n");
    printf("\n");
    printf("where:\n");
    printf("  modelfile   file containing LIBSVM classification model\n");
    printf("  train       files containing training data for drawing samples from\n");
    printf("                - .vec for vector data; .cls for class data\n");
    printf("                - need not be the same as originally used to train the model\n");
    printf("  border      name of acclerated model output file\n");
    //printf("                - .std contains normalization data (unless -a specified)\n");
    printf("  output      files containing the results of the classification:\n");
    printf("                - .cls for classes; .con for confidence ratings\n");
    printf("\n");
    printf("options:\n");
    //printf("  -n          option to normalise the data\n");
    //printf("  -u          borders data is stored in un-normalized coords\n");
    //printf("  -a normfile file containing normalization data\n");
    printf("  -i maxit1   maximum number of iterations when searching for class border (%d)\n", (int32_t) agf_global_borders_maxiter);
    printf("  -s n        number of times to sample the border (default=%d)\n", (int32_t) opt_args.n);
    printf("  -t tol      tolerance of border samples (default=%g)\n", (float) opt_args.tol);
    printf("  -A          ASCII format for test data and output\n");
    printf("  -M          LIBSVM format for test data and output\n");
    printf("  -E missing  missing value for LIBSVM features data\n");
    printf("  -N maxit3   maximum number of iterations in supernewton (%d, %d)\n", (int32_t) agf_global_weights_maxiter, (int32_t) agf_global_borders_maxiter);
    printf("\n");
    return INSUFFICIENT_COMMAND_ARGS;
  }

  ran_init();			//random numbers resolve ties

  //read in class borders:
  //"in-house" SVM predictor:
  classifier=new svm_multi<calc_t, cls_ta>(argv[0]);
  //printf("%d border vectors found: %s\n", ntrain, argv[0]);

  if (opt_args.asciiflag) {
    if (opt_args.Mflag) {
      ntest=read_svm(argv[1], test, cls, nvar, opt_args.missing, opt_args.Uflag);
      if (ntest==0) {
        fprintf(stderr, "An error occurred reading file, %s, in LIBSVM format.\n", argv[1]);
	exit(FILE_READ_ERROR);
      } else if (ntest<0) ntest=-ntest;
    } else {
      ntest=read_lvq(argv[1], test, cls, nvar);
      if (ntest<0) {
        fprintf(stderr, "An error occurred reading file, %s, in LVQ format.\n", argv[1]);
	exit(FILE_READ_ERROR);
      }
    }
    if (test==NULL) {
      fprintf(stderr, "An error occurred reading ASCII file, %s\n", argv[1]);
      exit(FILE_READ_ERROR);
    }
  } else {
    err=agf_read_train(argv[1], test, cls, ntest, nvar);
    if (err!=0) exit(err);
  }

  fprintf(stderr, "%d training vectors found in file %s\n", ntest, argv[1]);

  //normalization (currently disabled):
  if ((opt_args.uflag || opt_args.normflag) && opt_args.normfile==NULL) {
    opt_args.normfile=new char [strlen(argv[0])+5];
    sprintf(opt_args.normfile, "%s.std", argv[2]);
  }

  if (opt_args.normfile!=NULL) {
    mat=read_stats2(opt_args.normfile, ave, nvar1, nvar2);

    mat2=allocate_matrix<calc_t, int32_t>(nvar1, nvar2);
    b2=new calc_t[nvar1];

    for (int i=0; i<nvar1; i++) {
      b2[i]=ave[i];
      for (int j=0; j<nvar2; j++) mat2[i][j]=mat[i][j];
    }

    err=classifier->ltran(mat2, b2, nvar1, nvar2, opt_args.uflag);
    if (err!=0) exit(err);

    if (classifier->n_feat_t() != nvar) {
      fprintf(stderr, "svm_accelerate: Dimensions of classifier (%d) do not match those of training data (%d).\n",
                classifier->n_feat_t(), nvar);
      exit(DIMENSION_MISMATCH);
    }
  } else {
    if (classifier->n_feat() != nvar) {
      fprintf(stderr, "svm_accelerate: Dimension of classifier (%d) does not match dimension of training data (%d).\n",
                classifier->n_feat(), nvar);
      exit(DIMENSION_MISMATCH);
    }
  }

  //calculations are done in double precision:
  test2=allocate_matrix<calc_t, int32_t>(ntest, nvar);
  for (nel_ta i=0; i<ntest; i++) {
    for (dim_ta j=0; j<nvar; j++) test2[i][j]=test[i][j];
  }

  accel=new borders1v1<calc_t, cls_ta>(classifier, test2, cls, nvar, ntest,
		  opt_args.n, opt_args.tol);

  fs=fopen(argv[2], "w");
  if (fs == NULL) {
    fprintf(stderr, "Unable to open file, %s, for writing\n", argv[2]);
    return UNABLE_TO_OPEN_FILE_FOR_WRITING;
  }
  accel->save(fs);
  fclose(fs);
  
  //clean up:
  delete [] test[0];
  delete [] test;
  delete [] cls;
  if (mat!=NULL) {
    delete_matrix(mat);
    delete [] ave;
    delete_matrix(mat2);
    delete [] b2;
  }
  if (opt_args.normfile!=NULL) delete [] opt_args.normfile;
  delete classifier;
  delete accel;

  delete_matrix(test2);

  ran_end();

}


