
//Copyright (C) 2007 Peter Mills.  All rights reserved.

#include <math.h>
#include <string.h>
#include <stdio.h>

#include "agf_lib.h"

#define KDEFAULT 20

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

  //set defaults and parse command line options:
  opt_args.k=KDEFAULT;

  exit_value=agf_parse_command_opts(argc, argv, "k:n", &opt_args);
  if (exit_value==FATAL_COMMAND_OPTION_PARSE_ERROR) return exit_value;

  //parse the command line arguments:
  if (argc < 3) {
    printf("\n");
    printf("Syntax:   pdf_knn [-n] [-k k] train test output\n");
    printf("\n");
    printf("arguments:\n");
    printf("    train    file containing training vectors\n");
    printf("    test     file containing test vectors\n");
    printf("    output   binary output file\n");
    printf("\n");
    printf("options:\n");
    printf("    -k k     number of nearest neighbours (default=%d)\n", KDEFAULT);
    printf("    -n       normalize the data\n");
    printf("\n");
    return INSUFFICIENT_COMMAND_ARGS;
  }

  diagfs=stderr;

  trainfile=argv[0];
  testfile=argv[1];
  outfile=argv[2];

  //read in the training data:
  train=read_vecfile(trainfile, ntrain, nvar);
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

  //check range of k:
  if (opt_args.k <= 0 || opt_args.k >= ntrain) {
    fprintf(stderr, "Error, pdf_knn: parameter k=%d out of range.\n", opt_args.k);
    return PARAMETER_OUT_OF_RANGE;
  }

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

    norm_vec(train, nvar, ntrain, ave, std);
    norm_vec(test, nvar, ntest, ave, std);
  }

  real_a pcor=1;
  if (opt_args.normflag==1) {
    for (dim_ta j=0; j<nvar; j++) pcor=pcor*std[j];
  }

  //begin the pdf estimation scheme:
  result=new real_a[ntest];

  for (nel_ta i=0; i<ntest; i++) {
    //get the estimate:
    result[i]=knn_pdf(train, nvar, ntrain, test[i], opt_args.k);
    result[i]/=pcor;

    //print results to standard out:
    printf("%8d %12.6g\n", i, result[i]);
  } 

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


