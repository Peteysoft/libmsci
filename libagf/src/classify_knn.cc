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

//
// Command line utility for KNN statistical classification.
//
// See help screen for details.
//

#include <math.h>
#include <string.h>
#include <stdio.h>

#include "agf_lib.h"

#define KDEFAULT 11 

using namespace std;
using namespace libagf;
using namespace libpetey;

int main(int argc, char *argv[]) {
  char *testfile;		//test data
  char *outfile;		//output classes
  char *confile;		//output confidence ratings
  FILE *fs;
  FILE *diagfs;			//print diagnostics to this file stream

  nel_ta ntrain;		//number of training data points
  dim_ta nvar;			//number of variables
  cls_ta nclass;		//number of classes
  nel_ta ntest;			//number of test data points
  dim_ta nvar1;

  real_a **train;	//training data vectors
  cls_ta *cls;		//training data classes
  real_a **test;		//test data vectors
  cls_ta *result;	//results of classification

  char pdformat[10];

  agf_command_opts opt_args;

  real_a *con;		//confidence rating for the classification

  int errcode;
  int exit_value;

  real_a *std, *ave;

  //diagnostics:
  real_a *pdf;		//the returned pdfs

  //set defaults and parse command line options:
  opt_args.k=KDEFAULT;

  exit_value=0;
  exit_value=agf_parse_command_opts(argc, argv, "a:k:m:nj", &opt_args);
  if (exit_value==FATAL_COMMAND_OPTION_PARSE_ERROR) return exit_value;

  //parse the command line arguments:
  if (argc < 3) {
    printf("\n");
    printf("Syntax:       classify_knn [-n] [-k k] train test output\n");
    printf("\n");
    printf("arguments:\n");
    printf("    train     files containing training data:\n");
    printf("               .vec for vectors\n");
    printf("               .cls for classes\n");
    printf("    test      file containing vector data to be classified\n");
    printf("    output    files containing the results of the classification:\n");
    printf("               .cls for classification results\n");
    printf("               .con for confidence ratings\n");
    printf("\n");
    printf("options:\n");
    printf("  -k k        number of nearest neighbours to use in each estimate\n");
    printf("               (default=%d)\n", opt_args.k);
    printf("  -n          normalize the data\n");
    printf("  -j          print joint instead of cond. prob. to stdout\n");
    printf("  -m metric   type of metric to use\n");
    printf("       (%s)\n", global_metric_type);
    printf("\n");
    return INSUFFICIENT_COMMAND_ARGS;
  }

  diagfs=stderr;

  if (opt_args.jointflag) strcpy(pdformat, "%12.6g "); else strcpy(pdformat, "%8.6f ");

  //read in the training data:
  errcode=agf_read_train(argv[0], train, cls, ntrain, nvar);
  if (errcode != 0) return errcode;
  fprintf(diagfs, "%d training vectors found: %s\n", ntrain, argv[0]);

  if (opt_args.k <= 0 || opt_args.k >= ntrain) {
    fprintf(stderr, "Error: parameter k=%d out of range.\n", opt_args.k);
    return PARAMETER_OUT_OF_RANGE;
  }

  //count the number of classes:
  nclass=1;
  for (nel_ta i=0; i<ntrain; i++) if (cls[i]>=nclass) nclass=cls[i]+1;
  fprintf(diagfs, "Performing a %d class classification\n", nclass);

  testfile=argv[1];

  outfile=new char[strlen(argv[2])+5];
  strcpy(outfile, argv[2]);
  strcat(outfile, ".cls");

  confile=new char[strlen(argv[2])+5];
  strcpy(confile, argv[2]);
  strcat(confile, ".con");

  test=read_vecfile(testfile, ntest, nvar1);
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

  std=new real_a[nvar];
  ave=new real_a[nvar];
  if (opt_args.normflag == 1) {
    //calculate the averages and standard deviations:
    calc_norm(train, nvar, ntrain, ave, std);

    print_stats(stdout, ave, std, nvar);
    printf("\n");

    //fprintf(diagfs, "Normalizing the data.  ");
    //normalize the data:
    norm_vec(train, nvar, ntrain, ave, std);
    norm_vec(test, nvar, ntest, ave, std);

    if (opt_args.normfile != NULL) {
      fs=fopen(opt_args.normfile, "w");
      print_stats(fs, ave, std, nvar);
    }
  } else if (opt_args.normfile != NULL) {
    errcode=read_stats(opt_args.normfile, ave, std, nvar);
    if (errcode == 0) {
      norm_vec(train, nvar, ntrain, ave, std);
      norm_vec(test, nvar, ntest, ave, std);
    }
  }

  //begin the classification scheme:
  result=new cls_ta[ntest];
  con=new real_a[ntest];
  pdf=new real_a[nclass];

  for (nel_ta i=0; i<ntest; i++) {
    real_a p_x;
    //get classification result
    result[i]=knn(global_metric2_pointer_list[opt_args.metrictype], 
			train, nvar, ntrain, cls, nclass, 
			test[i], opt_args.k, pdf, opt_args.jointflag);
    if (opt_args.jointflag) {
      p_x=0;
      for (cls_ta j=0; j<nclass; j++) p_x+=pdf[j];
      con[i]=(nclass*pdf[result[i]]/p_x-1)/(nclass-1);
    } else {
      con[i]=(nclass*pdf[result[i]]-1)/(nclass-1);
    }


    //print results to standard out:
    printf("%8d %4d   ", i, result[i]);
    for (cls_ta j=0; j<nclass; j++) printf(pdformat, pdf[j]);
    printf("\n");

  }

  //write the results to a file:
  fs=fopen(outfile, "w");
  if (fs == NULL) {
    fprintf(fs, "Unable to open file for writing: %s\n", outfile);
    return UNABLE_TO_OPEN_FILE_FOR_WRITING;
  }
  fwrite(result, sizeof(cls_ta), ntest, fs);
  fclose(fs);

  fs=fopen(confile, "w");
  if (fs == NULL) {
    fprintf(fs, "Unable to open file for writing: %s\n", confile);
    return UNABLE_TO_OPEN_FILE_FOR_WRITING;
  }
  fwrite(con, sizeof(real_a), ntest, fs);
  fclose(fs);

  //clean up:
  delete [] result;
  delete [] train[0];
  delete [] train;
  delete [] test[0];
  delete [] test;
  delete [] cls;
  delete [] pdf;

  delete std;
  delete ave;

  return exit_value;
}


