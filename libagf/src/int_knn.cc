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

#define K_DEFAULT 5

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
  dim_ta nvar1, n1;

  real_a **train;	//training data vectors
  real_a *ord;		//training data classes
  real_a **test;		//test data vectors
  real_a *result;	//results of classification

  agf_command_opts opt_args;

  real_a *err;		//confidence rating for the classification

  int errcode;
  int exit_value;

  real_a *std, *ave;

  //set defaults and parse command line options:
  opt_args.k=K_DEFAULT;

  exit_value=0;
  exit_value=agf_parse_command_opts(argc, argv, "k:m:ne", &opt_args);
  if (exit_value==FATAL_COMMAND_OPTION_PARSE_ERROR) return exit_value;

  //parse the command line arguments:
  if (argc < 3) {
    printf("\n");
    printf("Syntax:	classify [-n] [-e] [-k k] train test output\n");
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
    printf("    -k k     number of nearest neighbours to use in each estimate\n");
    printf("               (default = %d)\n", opt_args.k);
    printf("    -n       normalize the data\n");
    printf("    -e       return error estimates\n");
    printf("  -m metric   type of metric to use\n");
    printf("       (%s)\n", global_metric_type);
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
  if (opt_args.normflag == 1) {
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
    }
  }

  //check range of k:
  if (opt_args.k <= 0 || opt_args.k >= ntrain) {
    fprintf(stderr, "Warning: parameter k=%d out of range.\n", opt_args.k);
    exit_value=PARAMETER_OUT_OF_RANGE;
    exit(exit_value);
  }

  //fprintf(diagfs, "Objective total weight: Wc=%f\n", opt_args.Wc);
  
  //begin the classification scheme:
  result=new real_a[ntest];

  if (opt_args.errflag == 1) {
    err=new real_a[ntest];
    for (nel_ta i=0; i<ntest; i++) {
      //get interpolation result
      result[i]=int_knn(global_metric2_pointer_list[opt_args.metrictype], 
		      train, nvar, ntrain, 
		      ord, test[i], opt_args.k, err[i]);

      //print results to standard out:
      printf("%8d %12.6g %8.2g\n", i, result[i], err[i]);

    }
  } else {
    err=NULL;
    for (nel_ta i=0; i<ntest; i++) {
      //get interpolation result
      result[i]=int_knn(global_metric2_pointer_list[opt_args.metrictype], 
		      train, nvar, ntrain, 
		      ord, test[i], opt_args.k);

      //print results to standard out:
      printf("%8d %12.6g\n", i, result[i]);

    }
  }

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
  if (err != NULL) delete [] err;

  if (std != NULL) delete std;
  if (ave != NULL) delete ave;

  return exit_value;

}
