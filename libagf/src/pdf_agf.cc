
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
  FILE *fs;
  FILE *diagfs;			//print diagnostics to this file stream

  nel_ta ntrain;		//number of training data points
  dim_ta nvar;		//number of variables
  nel_ta ntest;		//number of test data points
  dim_ta nvar1;

  real_a **train;	//training data vectors
  real_a **test;		//test data vectors
  real_a *result;	//results of pdf estimation

  agf_command_opts opt_args;

  int exit_value;

  real_a *std, *ave;
  real_a vart;

  //diagnostics:
  agf_diag_param diag_param;
  iter_ta min_nd, max_nd, total_nd;
  real_a min_f, max_f, total_f;
  real_a min_W, max_W, total_W;

  //set defaults and parse command line options:
  opt_args.k=-1;
  opt_args.W2=W_DEFAULT_BORDERS;

  exit_value=agf_parse_command_opts(argc, argv, "i:I:l:k:W:v:V:n", &opt_args);
  if (exit_value==FATAL_COMMAND_OPTION_PARSE_ERROR) return exit_value;

  //parse the command line arguments:
  if (argc < 3) {
    printf("\n");
    printf("Syntax:	pdf_agf [-n] [-k k] [-W Wc] train test output\n");
    printf("\n");
    printf("arguments:\n");
    printf("    train    file containing training vectors\n");
    printf("    test     file containing test vectors\n");
    printf("    output   binary output file\n");
    printf("\n");
    printf("options:\n");
    printf("    -v var1  lower filter variance\n");
    printf("               --default is to use total variance of the data/n^(2./D)\n");
    printf("    -V var2  upper filter variance bracket\n");
    printf("               --default is to use total variance of the data\n");
    printf("    -k k     number of nearest neighbours to use in each estimate\n");
    printf("               --default is to use all of the data\n");
    printf("    -W Wc    objective total weight (default=%6.1f)\n", opt_args.W2);
    printf("    -n       normalize the data\n");
    printf("\n");
    return INSUFFICIENT_COMMAND_ARGS;
  }

  diagfs=stderr;

  trainfile=argv[0];
  testfile=argv[1];
  outfile=argv[2];

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
        opt_args.var[1]=nvar;
        //fprintf(diagfs, "Using %g for initial filter variance.\n", 1.*nvar);
        //fprintf(diagfs, "\n");
      }
      if (opt_args.var[0] <= 0) {
        opt_args.var[0]=nvar/pow(ntrain, 2./nvar);
      }
    } else {
      //if the initial filter variance is not set, set it to the total
      //variance of the data:
      vart=0;
      for (dim_ta i=0; i<nvar; i++) vart+=std[i]*std[i];
      if (opt_args.var[1] <= 0) {
        opt_args.var[1]=vart;
        fprintf(diagfs, "Using %10.3g for upper filter variance bracket\n\n", opt_args.var[1]);
      }
      if (opt_args.var[0] <= 0) {
        opt_args.var[0]=vart/pow(ntrain, 2./nvar);
        fprintf(diagfs, "Using %10.3g for initial filter variance bracket\n\n", opt_args.var[0]);
      }
    }
  }

  //check range of k:
  if (opt_args.k <= opt_args.W2 || opt_args.k >= ntrain) {
    if (opt_args.k != -1) {
      fprintf(stderr, "Warning: parameter k=%d out of range.  Using all the training data\n", opt_args.k);
      opt_args.k=-1;
      exit_value=PARAMETER_OUT_OF_RANGE;
    }
    //fprintf(diagfs, "Using all the training data");
    //fprintf(diagfs, "\n");
  } else {
    //fprintf(diagfs, "Number of nearest neighbours to use: k=%ld\n", opt_args.k);
  }

  //fprintf(diagfs, "Objective total weight: Wc=%f\n", opt_args.Wc);
  
  //begin the classification scheme:
  result=new real_a[ntest];

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

  real_a pcor=1;
  if (opt_args.normflag==1) {
    for (dim_ta j=0; j<nvar; j++) pcor=pcor*std[j];
  }

  if (opt_args.k == -1) {
    for (nel_ta i=0; i<ntest; i++) {
      //get classification result
      result[i]=agf_calc_pdf(train, nvar, ntrain, test[i],
		      opt_args.var, opt_args.W2, &diag_param);
      result[i]/=pcor;

      //print results to standard out:
      printf("%8d %12.6g\n", i, result[i]);

      //calculate diagnostics:
      if (diag_param.nd < min_nd) min_nd=diag_param.nd;
      if (diag_param.nd > max_nd) max_nd=diag_param.nd;
      total_nd+=diag_param.nd;

      if (diag_param.W < min_W) min_W=diag_param.W;
      if (diag_param.W > max_W) max_W=diag_param.W;
      total_W+=diag_param.W;
    }
  } else {
    for (nel_ta i=0; i<ntest; i++) {
      //get classification result
      result[i]=agf_calc_pdf(train, nvar, ntrain, test[i], 
		      opt_args.var, opt_args.k, opt_args.W2, &diag_param);
      result[i]/=pcor;

      //print results to standard out:
      printf("%8d %12.6g\n", i, result[i]);

      //calculate diagnostics:
      if (diag_param.nd < min_nd) min_nd=diag_param.nd;
      if (diag_param.nd > max_nd) max_nd=diag_param.nd;
      total_nd+=diag_param.nd;

      if (diag_param.f < min_f) min_f=diag_param.f;
      if (diag_param.f > max_f) max_f=diag_param.f;
      total_f+=diag_param.f;

      if (diag_param.W < min_W) min_W=diag_param.W;
      if (diag_param.W > max_W) max_W=diag_param.W;
      total_W+=diag_param.W;
    }
  }

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

  //write the results to a file:
  fs=fopen(outfile, "w");
  if (fs==NULL) {
    fprintf(stderr, "Unable to open file for writing: %s\n", outfile);
    return UNABLE_TO_OPEN_FILE_FOR_WRITING;
  }
  fwrite(result, sizeof(real_a), ntest, fs);
  fclose(fs);

  //clean up:
  delete [] result;
  delete [] train[0];
  delete [] train;
  delete [] test[0];
  delete [] test;

  if (std != NULL) delete [] std;
  if (ave != NULL) delete [] ave;

  return exit_value;

}


