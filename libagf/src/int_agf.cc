//
// This software is released under the following terms:
//
// 1. No commercial use.
// 2. Copies and derivative works are free to use and modify.
// 3. Attribution must be given to all contributors of both original and derivative works.
//
// Authors:
//
// 2017-07-16 Peter Mills: added license information 
//

#include <math.h>
#include <string.h>
#include <stdio.h>

#include "agf_lib.h"

#define WC_DEFAULT 5.

using namespace std;
using namespace libagf;

int main(int argc, char *argv[]) {
  char *vecfile;		//coordinate file
  char *ordfile;		//ordinate file
  char *testfile;		//test data
  char *outfile;		//output classes
  char *errfile;		//output confidence ratings
  FILE *fs;
  FILE *diagfs;			//print diagnostics to this file stream

  nel_ta ntrain;		//number of training data points
  dim_ta nvar;		//number of variables
  cls_ta nclass;		//number of classes
  nel_ta ntest;		//number of test data points
  dim_ta nvar1;
  nel_ta n1;

  real_a **train;	//training data vectors
  real_a *ord;		//training data classes
  real_a **test;		//test data vectors
  real_a *result;	//results of classification

  agf_command_opts opt_args;

  real_a *err;		//confidence rating for the classification

  int errcode;
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
  opt_args.W2=WC_DEFAULT;

  exit_value=0;
  exit_value=agf_parse_command_opts(argc, argv, "i:I:l:k:W:v:V:ne", &opt_args);
  if (exit_value==FATAL_COMMAND_OPTION_PARSE_ERROR) return exit_value;

  //parse the command line arguments:
  if (argc < 3) {
    printf("\n");
    printf("Syntax:      classify [-n] [-e] [-v var1] [-V var2] [-k k] [-W Wc]\n");
    printf("                     train test output\n");
    printf("\n");
    printf("arguments:\n");
    printf("    train    files containing training data:\n");
    printf("               .vec for coordinate vectors\n");
    printf("               .dat for ordinates\n");
    printf("    test     file containing vector data to be interpolated\n");
    printf("    output   binary output files:\n");
    printf("               .dat for interpolation results\n");
    printf("               .err for error estimates (if applicable)\n");
    printf("\n");
    printf("options:\n");
    printf("    -v var1  lower filter variance bracket\n");
    printf("               --default is to use the total variance of the data/n^(2./D)\n");
    printf("    -V var2  upper filter variance bracket\n");
    printf("               --default is to use the total variance of the data\n");
    printf("    -k k     number of nearest neighbours to use in each estimate\n");
    printf("               --default is to use all of the data\n");
    printf("    -W Wc    objective total weight (default=%6.1f)\n", opt_args.W2);
    printf("    -n       normalize the data\n");
    printf("    -e       return error estimates\n");
    printf("\n");
    return INSUFFICIENT_COMMAND_ARGS;
  }

  diagfs=stderr;

  //read in the training data:
  vecfile=new char[strlen(argv[0])+5];
  strcpy(vecfile, argv[0]);
  strcat(vecfile, ".vec");

  ordfile=new char[strlen(argv[0])+5];
  strcpy(ordfile, argv[0]);
  strcat(ordfile, ".dat");

  //read in the training data:
  train=read_vecfile(vecfile, ntrain, nvar);
  if (nvar == -1) {
    fprintf(stderr, "Error reading file: %s\n", vecfile);
    return FILE_READ_ERROR;
  }
  if (ntrain == -1) {
    fprintf(stderr, "Error reading file: %s\n", vecfile);
    return ALLOCATION_FAILURE;
  }
  if (train==NULL) {
    fprintf(stderr, "Unable to open input file: %s\n", vecfile);
    return UNABLE_TO_OPEN_FILE_FOR_READING;
  }

  ord=read_datfile(ordfile, n1);

  if (n1 == -1) {
    fprintf(stderr, "Error reading file: %s\n", ordfile);
    delete [] train[0];
    delete [] train;
    return FILE_READ_ERROR;
  }
  if (ord == NULL) {
    fprintf(stderr, "Unable to open input file: %s\n", ordfile);
    delete [] train[0];
    delete [] train;
    return UNABLE_TO_OPEN_FILE_FOR_READING;
  }
  if (n1!=ntrain) {
    fprintf(stderr, "Error: number of samples in coordinate and ordinate files do not agree:\n");
    fprintf(stderr, "       %d in %s vs. %d in %s\n", ntrain, vecfile, n1, ordfile);
    delete [] train[0];
    delete [] train;
    delete [] ord;
    return SAMPLE_COUNT_MISMATCH;
  }

  fprintf(diagfs, "%d training vectors found: %s\n", ntrain, argv[0]);

  testfile=argv[1];

  outfile=new char[strlen(argv[2])+5];
  strcpy(outfile, argv[2]);
  strcat(outfile, ".dat");

  errfile=new char[strlen(argv[2])+5];
  strcpy(errfile, argv[2]);
  strcat(errfile, ".err");

  test=read_vecfile(testfile, ntest, nvar1);
  if (nvar == -1) {
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
        opt_args.var[0]=1./pow(ntrain, 2./nvar);
      }
    } else {
      //if the initial filter variance is not set, set it to the total
      //variance of the data:
      vart=0;
      for (dim_ta i=0; i<nvar; i++) vart+=std[i]*std[i];
      if (opt_args.var[1] <= 0) {
        opt_args.var[1]=vart;
        fprintf(diagfs, "Using %10.3g for upper filter variance bracket\n\n", opt_args.var[1]);
        //fprintf(diagfs, "Using 1 for initial filter variance.\n");
        //fprintf(diagfs, "\n");
      }
      if (opt_args.var[0] <= 0) {
        opt_args.var[0]=vart/pow(ntrain, 2./nvar);
        fprintf(diagfs, "Using %10.3g for lower filter variance bracket\n\n", opt_args.var[0]);
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
  err=new real_a[ntest];

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

  if (opt_args.k == -1) {
    if (opt_args.errflag == 1) {
      for (nel_ta i=0; i<ntest; i++) {
        //get interpolation result
        result[i]=adgaf_err(train, nvar, ord, ntrain,  
		      test[i], opt_args.var, opt_args.W2, err[i], &diag_param);

        //print results to standard out:
        printf("%8d %12.6g %8.2g\n", i, result[i], err[i]);

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
        //get interpolation result
        result[i]=adgaf(train, nvar, ord, ntrain,  
		      test[i], opt_args.var, opt_args.W2, &diag_param);

        //print results to standard out:
        printf("%8d %12.6g\n", i, result[i]);

        //calculate diagnostics:
        if (diag_param.nd < min_nd) min_nd=diag_param.nd;
        else if (diag_param.nd > max_nd) max_nd=diag_param.nd;
        total_nd+=diag_param.nd;

        if (diag_param.W < min_W) min_W=diag_param.W;
        else if (diag_param.W > max_W) max_W=diag_param.W;
        total_W+=diag_param.W;
      }
    }
  } else {
    if (opt_args.errflag == 1) {
      for (nel_ta i=0; i<ntest; i++) {
        //get interpolation result
        result[i]=adgaf_err(train, nvar, ord, ntrain,  
		      test[i], opt_args.var, opt_args.W2, err[i], &diag_param);

        //print results to standard out:
        printf("%8d %12.6g %8.2g\n", i, result[i], err[i]);

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
    } else {
      for (nel_ta i=0; i<ntest; i++) {
        //get interpolation result
        result[i]=adgaf(train, nvar, ord, ntrain,  
		      test[i], opt_args.var, opt_args.W2, &diag_param);

        //print results to standard out:
        printf("%8d %12.6g\n", i, result[i]);

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
  if (fs == NULL) {
    fprintf(fs, "Unable to open file for writing: %s\n", outfile);
    return UNABLE_TO_OPEN_FILE_FOR_WRITING;
  }
  fwrite(result, sizeof(real_a), ntest, fs);
  fclose(fs);

  fs=fopen(errfile, "w");
  if (fs == NULL) {
    fprintf(fs, "Unable to open file for writing: %s\n", errfile);
    return UNABLE_TO_OPEN_FILE_FOR_WRITING;
  }
  fwrite(err, sizeof(real_a), ntest, fs);
  fclose(fs);

  //clean up:
  delete [] result;
  delete [] train[0];
  delete [] train;
  delete [] test[0];
  delete [] test;
  delete [] ord;

  if (std != NULL) delete std;
  if (ave != NULL) delete ave;

  return exit_value;

}
