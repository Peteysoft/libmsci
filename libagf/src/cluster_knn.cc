
//Copyright (C) 2007 Peter Mills.  All rights reserved.

#include <math.h>
#include <string.h>
#include <stdio.h>

#include "agf_lib.h"

#define KDEFAULT 20
#define PMIN_DEFAULT 0.1

using namespace std;
using namespace libagf;

int main(int argc, char *argv[]) {
  char *trainfile;		//training data
  char *outfile;		//output classes
  FILE *fs;
  FILE *diagfs;			//print diagnostics to this file stream

  nel_ta ntrain;		//number of training data points
  dim_ta nvar;			//number of variables

  real_a **train;		//training data vectors
  cls_ta *result;		//results of cluster analysis

  agf_command_opts opt_args;

  int exit_value;

  //set defaults and parse command line options:
  opt_args.k=KDEFAULT;

  exit_value=agf_parse_command_opts(argc, argv, "p:k:nS:", &opt_args);
  if (exit_value==FATAL_COMMAND_OPTION_PARSE_ERROR) return exit_value;

  //parse the command line arguments:
  if (argc < 2) {
    printf("\n");
    printf("Syntax:   cluster_knn [-n] [-k k] [-p pt] train output\n");
    printf("\n");
    printf("arguments:\n");
    printf("    train    file containing training vectors\n");
    printf("    output   binary output file containing classes\n");
    printf("\n");
    printf("options:\n");
    printf("    -k k     number of nearest neighbours (default=%d)\n", KDEFAULT);
    printf("    -p pt    threshold density\n");
    printf("    -n       normalize the data\n");
    printf("    -S nsv   perform SVD: number of singular values to keep\n");
    printf("\n");
    return INSUFFICIENT_COMMAND_ARGS;
  }

  diagfs=stderr;

  trainfile=argv[0];
  outfile=argv[1];

  //read in the training data:
  //if we need a normalization file and one hasn't been named, 
  //construct the name:
  if ((opt_args.svd>0 || opt_args.normflag) && opt_args.normfile == NULL) {
    opt_args.normfile=new char[strlen(argv[1])+5];
    sprintf(opt_args.normfile, "%s.std", argv[1]);
  }

  //get the training co-ordinate data, pre-process if necessary
  train=agf_get_features(trainfile, &opt_args, ntrain, nvar, 1);

  fprintf(diagfs, "%d training vectors found in file: %s\n", ntrain, trainfile);

  //check range of k:
  if (opt_args.k <= 0 || opt_args.k >= ntrain) {
    fprintf(stderr, "Error, pdf_knn: parameter k=%d out of range.\n", opt_args.k);
    return PARAMETER_OUT_OF_RANGE;
  }

  //begin the pdf estimation scheme:
  result=new cls_ta[ntrain];

  cluster_knn(train, nvar, ntrain, opt_args.k, opt_args.pmin, result);

  for (nel_ta i=0; i<ntrain; i++) {
    //print results to standard out:
    printf("%8d %6d\n", i, result[i]);
  } 

  //write the results to a file:
  fs=fopen(outfile, "w");
  if (fs==NULL) {
    fprintf(stderr, "Unable to open file for writing: %s\n", outfile);
    return UNABLE_TO_OPEN_FILE_FOR_WRITING;
  }
  fwrite(result, sizeof(cls_ta), ntrain, fs);
  fclose(fs);

  //clean up:
  delete [] result;
  delete [] train[0];
  delete [] train;

  return exit_value;

}


