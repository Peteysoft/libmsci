
//Copyright (C) 2007 Peter Mills.  All rights reserved.

#include <math.h>
#include <string.h>
#include <stdio.h>

#include "agf_lib.h"

using namespace std;
using namespace libpetey;
using namespace libagf;

int main(int argc, char *argv[]) {
  char *testfile;		//test data
  char *trainfile;		//training data
  char *outfile;		//output classes
  char *errfile;		//output error estimates
  FILE *fs;
  FILE *diagfs;			//print diagnostics to this file stream

  nel_ta ntrain;		//number of training data points
  dim_ta nvar;		//number of variables
  real_a *parm;		//parameters
  nel_ta ntest;		//number of test data points
  dim_ta nvar1;

  int flags=0;		//flags to pass to AGF optimal subroutine

  real_a **train;	//training data vectors
  real_a **test;		//test data vectors
  real_a *result;	//results of pdf estimation
  real_a *err;		//error estimate

  agf_command_opts opt_args;

  int exit_value;

  real_a *std, *ave;
  real_a vart;

  //set defaults and parse command line options:
  opt_args.W1=-1;
  opt_args.W2=-1;
  opt_args.var[0]=-1;
  opt_args.Qtype=-1;

  exit_value=0;
  exit_value=agf_parse_command_opts(argc, argv, "w:W:q:v:V:Q:nzR", &opt_args);
  if (exit_value==FATAL_COMMAND_OPTION_PARSE_ERROR) return exit_value;

  //parse the command line arguments:
  if (argc < 3) {
    printf("\n");
    printf("Syntax:	pdf_agf_opt [-n] [-v var0] [-q nt] [-Q alg] [-w wmin] [-W wmax] train test output\n");
    printf("\n");
    printf("arguments:\n");
    printf("    train    file containing training vectors\n");
    printf("    test     file containing test vectors\n");
    printf("    output   binary output file\n");
    printf("\n");
    printf("options:\n");
    printf("    -Q alg   algorithm for choosing filter variances:\n");
    printf("               0 = initial variance specified, half it each subsequent trial\n");
    printf("               1 = minimum and maximum variance specified\n");
    printf("               2 = minimum and maximum total weights specified\n");
    printf("               * options 1 and 2 use geometric progression for filter variance\n");
    printf("    -v var1  lower bracket for filter variance\n");
    printf("               --default is to use total variance of the data/n^(4/D)\n");
    printf("    -V var2  upper bracket for filter variance\n");
    printf("               --default is to use total variance of the data\n");
    printf("    -q nt    number of trials (default=%3d)\n", opt_args.nt);
    printf("    -w wmin  minimum value for W\n");
    printf("    -W wmax  maximum value for W\n");
    printf("    -n       normalize the data\n");
    printf("    -R       use sorted (unless -z specified) random filter variances\n");
    printf("               (exponential distribution as per -o 1 and 2)\n");
    printf("    -z       randomize the order of trials (filter variances not sorted)\n");
    printf("\n");
    return INSUFFICIENT_COMMAND_ARGS;
  }

  ran_init();

  diagfs=stderr;

  trainfile=argv[0];
  testfile=argv[1];
  outfile=new char[strlen(argv[2])+5];
  sprintf(outfile, "%s.dat", argv[2]);
  errfile=new char[strlen(argv[2])+5];
  sprintf(errfile, "%s.err", argv[2]);

  //read in the training data:
  train=read_vecfile<real_a>(trainfile, ntrain, nvar);
  if (nvar == -1) {
    fprintf(stderr, "Error reading file: %s\n", trainfile);
    return FILE_READ_ERROR;
  }
  if (ntrain == -1) {
    fprintf(stderr, "Error reading file: %s\n", trainfile);
    return ALLOCATION_FAILURE;
  }
  if (train==NULL) {
    fprintf(stderr, "Unable to open file for reading: %s\n", trainfile);
    return UNABLE_TO_OPEN_FILE_FOR_READING;
  }
  fprintf(diagfs, "%d training vectors found in file: %s\n", ntrain, trainfile);

  test=read_vecfile<real_a>(testfile, ntest, nvar1);
  if (nvar1 == -1) {
    fprintf(stderr, "Error reading file: %s\n", testfile);
    return FILE_READ_ERROR;
  }
  if (ntest == -1) {
    fprintf(stderr, "Error reading file: %s\n", testfile);
    return ALLOCATION_FAILURE;
  }
  if (test == NULL) {
    fprintf(stderr, "Unable to open input file: %s\n", testfile);
    return UNABLE_TO_OPEN_FILE_FOR_READING;
  }
  if (nvar != nvar1) {
    fprintf(stderr, "Error: dimensions of training and test data do not agree:\n");
    fprintf(stderr, "       %s: D=%d; %s: D=%d\n", argv[0], nvar, argv[1], nvar1);
    return DIMENSION_MISMATCH;
  }

  fprintf(diagfs, "%d test vectors found in file %s\n", ntest, testfile);

  std=NULL;
  ave=NULL;
  if (opt_args.normflag == 1 || opt_args.var[0] <= 0 || opt_args.var[1] <= 0) {
    //calculate the averages and standard deviations:
    std=new real_a[nvar];
    ave=new real_a[nvar];
    calc_norm(train, nvar, ntrain, ave, std);

    fprintf(diagfs, "\nStatistics:\n");
    fprintf(diagfs, "dim    average  std. dev.\n");
    for (dim_ta i=0; i<nvar; i++) {
      fprintf(diagfs, "%3d %10.6g %10.6g\n", i, ave[i], std[i]);
    }
    fprintf(diagfs, "\n");

    if (opt_args.normflag == 1) {
      //fprintf(diagfs, "Normalizing the data.  ");
      //normalize the data:
      norm_vec(train, nvar, ntrain, ave, std);
      norm_vec(test, nvar, ntest, ave, std);
      //if the initial filter variance is not set, set it to the total
      //variance of the data:
      if (opt_args.var[1] <= 0) {
        opt_args.var[1]=1;
        //fprintf(diagfs, "Using 1 for initial filter variance.\n");
        //fprintf(diagfs, "\n");
      }
      if (opt_args.var[0] <= 0) {
        opt_args.var[0]=1/pow(ntrain, 2./nvar);
        if (opt_args.Qtype==-1) opt_args.Qtype=1;
      }
    } else {
      //if the initial filter variance is not set, set it to the total
      //variance of the data:
      vart=0;
      for (dim_ta i=0; i<nvar; i++) vart+=std[i]*std[i];
      if (opt_args.var[1] <= 0) {
        opt_args.var[1]=vart;
        fprintf(diagfs, "Using %10.3g for upper filter variance\n\n", opt_args.var[1]);
        //fprintf(diagfs, "Using 1 for initial filter variance.\n");
        //fprintf(diagfs, "\n");
      }
      if (opt_args.var[0] <= 0) {
        opt_args.var[0]=vart/pow(ntrain, 2./nvar);
        if (opt_args.Qtype==-1) opt_args.Qtype=1;
        fprintf(diagfs, "Using %10.3g for lower filter variance\n\n", opt_args.var[0]);
      }
    }
  }

  //fprintf(diagfs, "Objective total weight: Wc=%f\n", opt_args.Wc);
  //decide which algorithm to use:
  if (opt_args.Qtype == -1 && (opt_args.W1 >0 || opt_args.W2 >0)) {
    opt_args.Qtype=2; 
  } else {
    opt_args.Qtype=1;
  }

  if (opt_args.W1 < 0) opt_args.W1=WMIN_DEFAULT;
  if (opt_args.W2 < 0) opt_args.W2=WMAX_DEFAULT;

  flags=opt_args.Qtype+10*opt_args.Rflag+100*opt_args.zflag;
  printf("%d %d\n", opt_args.Qtype, flags);
  
  //begin the classification scheme:
  result=new real_a[ntest];
  err=new real_a[ntest];

  //set parameters::
  parm=new real_a[4];
  parm[0]=opt_args.var[0];
  parm[1]=opt_args.var[1];
  parm[2]=opt_args.W1;
  parm[3]=opt_args.W2;

  printf("var0= %g\n", opt_args.var[0]);
  printf("varf= %g\n", opt_args.var[1]);
  printf("wmin= %g\n", opt_args.W1);
  printf("wmax= %g\n", opt_args.W2);

  real_a v2[2];
  v2[0]=opt_args.var[0];
  v2[1]=opt_args.var[1];

  real_a pcor=1;
  if (opt_args.normflag==1) {
    for (dim_ta j=0; j<nvar; j++) pcor=pcor*std[j];
  }

  for (nel_ta i=0; i<ntest; i++) {
    agf_diag_param diag_param;
    printf("x=[");
    for (dim_ta j=0; j<nvar; j++) printf("%g ", test[i][j]);
    printf("]\n");

    //get classification result
    result[i]=agf_calc_pdf_opt(train, nvar, ntrain, test[i],
		      parm, opt_args.nt, err[i], flags);
    result[i]/=pcor;

/*
    result[i]=agf_calc_pdf_opt2(train, nvar, ntrain, test[i],
		      parm, err[i], flags);
    real_a pdum=agf_calc_pdf(train, nvar, ntrain, test[i],
                      v2, opt_args.W2, &diag_param);
    printf("\n");
    printf("%g %g\n", result[i], pdum);
    printf("\n");
    result[i]=pdum;
*/

    /*result[i]=agf_calc_pdf_opt(train, nvar, ntrain, test[i],
		      parm, opt_args.nt, 3);
    result[i]=agf_calc_pdf_opt(train, nvar, ntrain, test[i],
		      parm, opt_args.nt, 2);*/
    if (opt_args.normflag==1) {
      for (dim_ta j=0; j<nvar; j++) result[i]/=std[j];
    }

    //print results to standard out:
    printf("%8d %12.6g %8.3g\n", i, result[i], err[i]);

  }

  //write the results to a file:
  fs=fopen(outfile, "w");
  if (fs==NULL) {
    fprintf(stderr, "Unable to open file for writing: %s\n", outfile);
    return UNABLE_TO_OPEN_FILE_FOR_WRITING;
  }
  fwrite(result, sizeof(real_a), ntest, fs);
  fclose(fs);

  fs=fopen(errfile, "w");
  if (fs==NULL) {
    fprintf(stderr, "Unable to open file for writing: %s\n", errfile);
    return UNABLE_TO_OPEN_FILE_FOR_WRITING;
  }
  fwrite(err, sizeof(real_a), ntest, fs);
  fclose(fs);

  //clean up:
  delete [] result;
  delete [] err;
  delete [] errfile;
  delete [] outfile;

  delete [] train[0];
  delete [] train;
  delete [] test[0];
  delete [] test;

  delete [] parm;

  if (std != NULL) delete std;
  if (ave != NULL) delete ave;

  return exit_value;

}


