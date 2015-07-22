#include <math.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>

#include "agf_lib.h"

using namespace std;
using namespace libagf;

int main(int argc, char *argv[]) {
  FILE *diagfs;

  real_a *ave;
  real_a *std;

  nel_ta ntrain;	//number of training data points
  dim_ta nvar;		//number of variables
  cls_ta nclass;	//number of classes
  nel_ta ntest;		//number of test data points
  dim_ta nvar1;

  real_a **train;	//training data vectors
  cls_ta *cls;		//training data classes
  real_a **test;	//test data vectors
  cls_ta *result;	//results of classification

  real_a *all;		//for memory de-allocation

  real_a *con;		//confidence rating for the classification

  nel_ta ntrue;
  nel_ta nfalse;
  real_a corr;
  nel_ta **acc_mat;	//accuracy matrix...

  agf_command_opts opt_args;	//parameters

  //diagnostics:
  real_a *pdf;

  int exit_code;
  int err_code;

  //set defaults and parse command line options:
  opt_args.k=K_DEFAULT_KNN;
  opt_args.ftest=F_DEFAULT;

  exit_code=0;
  exit_code=agf_parse_command_opts(argc, argv, "k:f:m:n", &opt_args);
  if (exit_code==FATAL_COMMAND_OPTION_PARSE_ERROR) return exit_code;

  if (argc < 1) {
    printf("\n");
    printf("syntax:  test_knn [-n] [-k k] [-f frac] train\n");
    printf("\n");
    printf("arguments:\n");
    printf("  train    = binary files containing locations of the samples:\n");
    printf("              .vec for vectors;\n");
    printf("              .cls for classes\n");
    printf("\n");
    printf("options:\n");
    printf("  -k k     = number of nearest neighbours to use in each calculation\n");
    printf("               --default=%d\n", K_DEFAULT_KNN);
    printf("  -f frac  = fraction of training data to use for testing (default=%g)\n", F_DEFAULT);
    printf("  -n       = option to normalise the data\n");
    printf("\n");
    return INSUFFICIENT_COMMAND_ARGS;
  }

  diagfs=stderr;

  ave=NULL;
  std=NULL;

  err_code=agf_read_train(argv[0], train, cls, ntrain, nvar);
  if (err_code!=0) return err_code;
  all=train[0];		//must do this because of the way the data is allocated...

  printf("%d training vectors found: %s\n", ntrain, argv[0]);

  if (opt_args.k <= 0 || opt_args.k >= ntrain) {
    fprintf(stderr, "Error: parameter k=%d out of range.\n", opt_args.k);
    return PARAMETER_OUT_OF_RANGE;
  }

  //count the number of classes:
  nclass=1;
  for (nel_ta i=0; i<ntrain; i++) if (cls[i]>=nclass) nclass++;

  //randomize the training data:
  fprintf(diagfs, "Randomizing the data...\n");
  randomize_vec(train, nvar, ntrain, cls);

  //calculate the number of test vectors:
  ntest=(nel_ta) (opt_args.ftest*ntrain);
  fprintf(diagfs, "Using %d vectors for testing\n", ntest);

  if (opt_args.normflag == 1) {
    //calculate the averages and standard deviations:
    std=new real_a[nvar];
    ave=new real_a[nvar];
    calc_norm(train, nvar, ntrain, ave, std);

    fprintf(diagfs, "Statistics:\n");
    fprintf(diagfs, "dim    average  std. dev.\n");
    for (dim_ta i=0; i<nvar; i++) {
      fprintf(diagfs, "%3d %10.6g %10.6g\n", i, ave[i], std[i]);
    }
    fprintf(diagfs, "\n");

    //normalize the data:
    norm_vec(train, nvar, ntrain, ave, std);
  }

  //begin the classification scheme:
  fprintf(diagfs, "Beginning classification...\n");
  result=new cls_ta[ntest];
  con=new real_a[ntest];
  pdf=new real_a[nclass];

  ntrue=0;
  nfalse=0;

  printf("  number cls ret   con nfalse  ntrue\n");
  for (nel_ta i=0; i<ntest; i++) {
    result[i]=knn(global_metric2_pointer_list[opt_args.metrictype], 
		    train+ntest, nvar, ntrain-ntest, cls+ntest,   
		    nclass, train[i], opt_args.k, pdf);
    con[i]=(nclass*pdf[result[i]]-1)/(nclass-1);
    if (cls[i] == result[i]) ntrue++; else nfalse++;

      //print results to standard out:
    printf("%8d %3d %3d %5.2f %6d %6d\n", i, cls[i], result[i], con[i], 
		    nfalse, ntrue);
  }

  uncertainty_coefficient(cls, result, ntest, nclass, stdout);
  printf("\n");

  printf("Accuracy: %f\n", 1.*ntrue/ntest);

  check_confidence(cls, result, con, nclass, ntest);
  
  //clean up:
  delete [] result;
  delete [] con;
  delete [] all;
  delete [] train;
  delete [] cls;
  delete [] pdf;

  if (std != NULL) delete [] std;
  if (ave != NULL) delete [] ave;

}


