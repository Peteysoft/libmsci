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
  dim_ta nvar;		//number of variables
  cls_ta nclass;		//number of classes
  nel_ta ntest;		//number of test data points
  dim_ta nvar1;

  real_a **train;	//training data vectors
  cls_ta *cls;		//training data classes
  real_a **test;		//test data vectors
  cls_ta *result;		//results of classification

  real_a vart;		//for calculating filter brackets

  char pdformat[10];

  agf_command_opts opt_args;

  real_a *con;		//confidence rating for the classification

  int errcode;
  int exit_value;

  real_a *std, *ave;

  //diagnostics:
  real_a *pdf;		//the returned pdfs
  real_a p_x;		//total probability
  agf_diag_param diag_param;
  iter_ta min_nd, max_nd, total_nd;
  real_a min_f, max_f, total_f;
  real_a min_W, max_W, total_W;

  //set defaults and parse command line options:
  opt_args.k=-1;
  opt_args.W2=W_DEFAULT;

  exit_value=0;
  exit_value=agf_parse_command_opts(argc, argv, "i:I:l:a:k:W:v:V:nj", &opt_args);
  if (exit_value==FATAL_COMMAND_OPTION_PARSE_ERROR) return exit_value;

  //parse the command line arguments:
  if (argc < 3) {
    printf("\n");
    printf("Syntax:      classify_a [-n] [-v var1] [-V var2] [-k k] [-W Wc] train test output\n");
    printf("\n");
    printf("arguments:\n");
    printf("    train    files containing training data:\n");
    printf("               .vec for vectors\n");
    printf("               .cls for classes\n");
    printf("    test     file containing vector data to be classified\n");
    printf("    output   files containing the results of the classification:\n");
    printf("               .cls for classification results\n");
    printf("               .con for confidence ratings\n");
    printf("\n");
    printf("options:\n");
    printf("    -v var1  first bracket of filter variance\n");
    printf("               --default is to use the total variance/n^(2/D)\n");
    printf("    -V var2  second bracket of filter variance/initial filter variance\n");
    printf("               --default is to use the total variance of the data\n");
    printf("    -k k     number of nearest neighbours to use in each estimate\n");
    printf("               --default is to use all of the data\n");
    printf("    -W Wc    objective total weight (default=%6.1f)\n", opt_args.W2);
    printf("    -n       normalize the data\n");
    printf("    -j       print joint instead of cond. prob. to stdout\n");
    printf("    -l tol   tolerance of W (default=%g)\n", WEIGHTS_TOL);
    printf("    -I/-i maxiter   maximum number of iterations in supernewton (%d)\n", WEIGHTS_MAXITER);
    printf("\n");
    return INSUFFICIENT_COMMAND_ARGS;
  }

  diagfs=stderr;

  if (opt_args.jointflag) strcpy(pdformat, "%12.6g "); else strcpy(pdformat, "%8.6f ");

  //read in the training data:
  errcode=agf_read_train(argv[0], train, cls, ntrain, nvar);
  if (errcode != 0) exit(errcode);
  fprintf(diagfs, "%d training vectors found: %s\n", ntrain, argv[0]);

  //count the number of classes:
  nclass=1;
  for (nel_ta i=0; i<ntrain; i++) if (cls[i]>=nclass) nclass=cls[i]+1;
  if (nclass < 2) {
    fprintf(stderr, "Cannot perform classifications with less than two classes!\n");
    return PARAMETER_OUT_OF_RANGE;
  }
  fprintf(diagfs, "Performing a %d class classification\n", nclass);

  testfile=argv[1];

  outfile=new char[strlen(argv[2])+5];
  strcpy(outfile, argv[2]);
  strcat(outfile, ".cls");

  confile=new char[strlen(argv[2])+5];
  strcpy(confile, argv[2]);
  strcat(confile, ".con");

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

  if (opt_args.normflag == 0 && opt_args.normfile != NULL) {
    errcode=read_stats(opt_args.normfile, ave, std, nvar);
    if (errcode == 0) {
      norm_vec(train, nvar, ntrain, ave, std);
      norm_vec(test, nvar, ntest, ave, std);
    } 
  }

  if (opt_args.normflag == 1 || opt_args.var[1] <= 0 || opt_args.var[0] <= 0) {
    //calculate the averages and standard deviations:
    std=new real_a[nvar];
    ave=new real_a[nvar];
    calc_norm(train, nvar, ntrain, ave, std);

    print_stats(stdout, ave, std, nvar);

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
      if (opt_args.var[0] <= 0) opt_args.var[0]=1./pow(ntrain, 2./nvar);
      if (opt_args.normfile != NULL) {
        fs=fopen(opt_args.normfile, "w");
        print_stats(fs, ave, std, nvar);
      }
    } else {
      //if the initial filter variance is not set, set it to the total
      //variance of the data:
      vart=0;
      for (dim_ta i=0; i<nvar; i++) vart+=std[i]*std[i];
      if (opt_args.var[0]<=0) opt_args.var[0]=vart/pow(ntrain, 2./nvar);
      if (opt_args.var[1]<=0) opt_args.var[1]=vart;
      //fprintf(diagfs, "Using %10.3g for initial filter variance\n\n", opt_args.var_0);
      fprintf(diagfs, "Bracketing filter variance with: [%10.3g, %10.3g]\n\n", opt_args.var[0], opt_args.var[1]);
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

  if (opt_args.k == -1) {
    for (nel_ta i=0; i<ntest; i++) {
      //get classification result
      result[i]=agf_classify(train, nvar, cls, ntrain, nclass, 
		      test[i], opt_args.var, opt_args.W2, pdf, &diag_param,
		      opt_args.jointflag);
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
      result[i]=agf_classify(train, nvar, cls, ntrain, nclass, 
		      test[i], opt_args.var, opt_args.k, opt_args.W2, 
		      pdf, &diag_param, opt_args.jointflag);
      if (opt_args.jointflag) {
        p_x=0;
        for (cls_ta j=0; j<nclass; j++) p_x+=pdf[j];
        con[i]=(nclass*pdf[result[i]]/p_x-1)/(nclass-1);
      } else {
        con[i]=(nclass*pdf[((cls_ta *)result)[i]]-1)/(nclass-1);
      }

      //print results to standard out:
      printf("%8d %4d   ", i, result[i]);
      for (cls_ta j=0; j<nclass; j++) printf(pdformat, pdf[j]);
      printf("\n");

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

  if (std != NULL) delete std;
  if (ave != NULL) delete ave;

  return exit_value;

}


