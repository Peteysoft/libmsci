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
  real_a vart;

  nel_ta ntrain;		//number of training data points
  dim_ta nvar;		//number of variables
  cls_ta nclass;		//number of classes
  nel_ta ntest;		//number of test data points
  dim_ta nvar1;

  real_a **train;	//training data vectors
  cls_ta *cls;		//training data classes
  real_a **test;		//test data vectors
  cls_ta *result;		//results of classification

  real_a *all;		//for memory de-allocation

  real_a *con;		//confidence rating for the classification

  nel_ta ntrue;
  nel_ta nfalse;
  real_a corr;
  nel_ta **acc_mat;	//accuracy matrix...

  int exit_code;
  int err_code;

  agf_command_opts opt_args;	//parameters

  //diagnostics:
  real_a *pdf;
  agf_diag_param diag_param;
  iter_ta min_nd, max_nd, total_nd;
  real_a min_f, max_f, total_f;
  real_a min_W, max_W, total_W;

  //set defaults and parse command line options:
  opt_args.k=-1;
  opt_args.W2=W_DEFAULT;
  opt_args.ftest=F_DEFAULT;

  exit_code=0;
  exit_code=agf_parse_command_opts(argc, argv, "i:I:l:v:V:k:W:f:n", &opt_args);
  if (exit_code==FATAL_COMMAND_OPTION_PARSE_ERROR) return exit_code;

  if (argc < 1) {
    printf("\n");
    printf("syntax:  test_classify [-n] [-V var1] [-V var2] [-k k] [-W Wc] [-f frac] train\n");
    printf("\n");
    printf("arguments:\n");
    printf("  train    = binary files containing locations of the samples:\n");
    printf("              .vec for vectors;\n");
    printf("              .cls for classes\n");
    printf("\n");
    printf("options:\n");
    printf("  -f frac  = fraction of training data to use for testing (default=%g)\n", opt_args.ftest);
    printf("  -v var1  = lower filter variance bracket\n");
    printf("               --default is to use the total variance of the data/n^(2/D)\n");
    printf("  -V var2  = upper filter variance bracket\n");
    printf("               --default is to use the total variance of the data\n");
    printf("  -k k     = number of nearest neighbours to use in each calculation\n");
    printf("               --default is to use all the data\n");
    printf("  -W Wc    = objective total weight (default=%g)\n", opt_args.W2);
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

  fprintf(diagfs, "%d training vectors found: %s\n", ntrain, argv[0]);

  //count the number of classes:
  nclass=1;
  for (nel_ta i=0; i<ntrain; i++) if (cls[i]>=nclass) nclass++;

  //randomize the training data:
  fprintf(diagfs, "Randomizing the data...\n");
  randomize_vec(train, nvar, ntrain, cls);

  //calculate the number of test vectors:
  ntest=(nel_ta) (opt_args.ftest*ntrain);
  fprintf(diagfs, "Using %d vectors for testing\n", ntest);

  if (opt_args.normflag == 1 || opt_args.var[0] <= 0 || opt_args.var[1] <= 0) {
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

    if (opt_args.normflag == 1) {
      //normalize the data:
      norm_vec(train, nvar, ntrain, ave, std);
      //if the initial filter variance is not set, set it to the total
      //variance of the data:
      if (opt_args.var[1] <= 0) opt_args.var[1]=1;
      if (opt_args.var[0] <= 0) opt_args.var[0]=1/pow(ntrain, 2./nvar);

    } else {
      //if the initial filter variance is not set, set it to the total
      //variance of the data:
      vart=0;
      for (dim_ta i=0; i<nvar; i++) vart+=std[i]*std[i];
      if (opt_args.var[1] <= 0) {
        opt_args.var[1]=vart;
        fprintf(diagfs, "Using %10.3g for second filter variance bracket\n\n", opt_args.var[1]);
      }
      if (opt_args.var[0] <= 0) {
        opt_args.var[0]=vart/pow(ntrain, 2./nvar);
        fprintf(diagfs, "Using %10.3g for first filter variance bracket\n\n", opt_args.var[0]);
      }

    }
  }

  if (opt_args.k <= opt_args.W2 || opt_args.k >= ntrain) {
    if (opt_args.k != -1) {
      fprintf(stderr, "Parameter k=%d out of range.  Using all the training data.\n", opt_args.k);
      opt_args.k=-1;
      exit_code=PARAMETER_OUT_OF_RANGE;
    }
  }

  //begin the classification scheme:
  fprintf(diagfs, "Beginning classification...\n");
  result=new cls_ta[ntest];
  con=new real_a[ntest];
  pdf=new real_a[nclass];

  //initialize diagnostic values:
  min_nd=1000;
  max_nd=0;
  total_nd=0;

  min_f=1;
  max_f=0;
  total_f=0;

  min_W=1000*opt_args.W2;
  max_W=0;
  total_W=0;

  ntrue=0;
  nfalse=0;

  printf("  number cls ret   con nfalse  ntrue\n");
  if (opt_args.k == -1) {
    for (nel_ta i=0; i<ntest; i++) {
      result[i]=agf_classify(train+ntest, nvar, cls+ntest, ntrain-ntest,  
		    nclass, train[i], opt_args.var, opt_args.W2, pdf, &diag_param);
      con[i]=(nclass*pdf[result[i]]-1)/(nclass-1);
      if (cls[i] == result[i]) ntrue++; else nfalse++;

      //print results to standard out:
      printf("%8d %3d %3d %5.2f %6d %6d\n", i, cls[i], result[i], con[i], 
		    nfalse, ntrue);

      //calculate diagnostics:
      if (diag_param.nd < min_nd) min_nd=diag_param.nd;
      else if (diag_param.nd > max_nd) max_nd=diag_param.nd;
      total_nd+=diag_param.nd;

      if (diag_param.W < min_W) min_W=diag_param.W;
      else if (diag_param.W > max_W) max_W=diag_param.W;
      total_W+=diag_param.W;
    }
  } else {
    for (nel_ta i=0; i<ntest; i++) {
      result[i]=agf_classify(train+ntest, nvar, cls+ntest, ntrain-ntest,  
		    nclass, train[i], opt_args.var, opt_args.k, opt_args.W2, pdf, &diag_param);
      con[i]=(nclass*pdf[result[i]]-1)/(nclass-1);
      if (cls[i] == result[i]) ntrue++; else nfalse++;

      //print results to standard out:
      printf("%8d %3d %3d %5.2f %6d %6d\n", i, cls[i], result[i], con[i], 
		    nfalse, ntrue);

      //calculate diagnostics:
      if (diag_param.nd < min_nd) min_nd=diag_param.nd;
      else if (diag_param.nd > max_nd) max_nd=diag_param.nd;
      total_nd+=diag_param.nd;

      if (diag_param.f < min_f) min_f=diag_param.f;
      else if (diag_param.f > max_f) max_f=diag_param.f;
      total_f+=diag_param.f;

      if (diag_param.W < min_W) min_W=diag_param.W;
      else if (diag_param.W > max_W) max_W=diag_param.W;
      total_W+=diag_param.W;
    }
  }
  printf("\n");

  //print out diagnostics:
  fprintf(diagfs, "\n");
  fprintf(diagfs, "diagnostic parameter          %8s   %8s   %8s\n",
                      "min", "max", "average");
  fprintf(diagfs, "iterations in agf_calc_w:      %8d   %8d %10.3g\n",
                   min_nd, max_nd, (real_a) total_nd/(real_a) ntest);
 if (opt_args.k != -1) {
    fprintf(diagfs, "value of f:                    %10.3g %10.3g %10.3g\n",
                   min_f, max_f, total_f/ntest);
  }
  fprintf(diagfs, "value of W:                    %10.2f %10.2f %10.2f\n",
                   min_W, max_W, total_W/ntest);
  fprintf(diagfs, "\n");

  //evaluate results:
  class_eval(cls, result, ntest);
  check_confidence(cls, result, con, nclass, NCONHIST);
  
  //clean up:
  delete [] result;
  delete [] all;
  delete [] train;
  delete [] cls;
  delete [] pdf;

  if (ave != NULL) delete [] ave;
  if (std != NULL) delete [] std;

}


