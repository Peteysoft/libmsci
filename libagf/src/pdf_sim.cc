
//Copyright (C) 2007 Peter Mills.  All rights reserved.

#include <math.h>
#include <string.h>
#include <stdio.h>

#include "agf_lib.h"
#include "randomize.h"

using namespace std;
using namespace libpetey;
using namespace libagf;

int main(int argc, char *argv[]) {
  char *trainfile;		//training data
  char *outfile;		//output classes
  FILE *fs;
  FILE *diagfs;			//print diagnostics to this file stream

  nel_ta ntrain;		//number of training data points
  dim_ta nvar;		//number of variables

  real_a **train;	//training data vectors
  real_a **result;	//results of pdf estimation

  agf_command_opts opt_args;

  //diagnostics:
  agf_diag_param diag_param;
  iter_ta min_nd, max_nd, total_nd;
  real_a min_f, max_f, total_f;
  real_a min_W, max_W, total_W;

  int exit_value;

  real_a *std, *ave;
  real_a vart;

  real_a pold, pnew;
  real_a rann;
  nel_ta ntest;

  real_a scale=1.0;		//scaling paramter for random walk

  //set defaults and parse command line options:
  opt_args.k=-1;
  opt_args.W2=W_DEFAULT_BORDERS;

  exit_value=0;
  exit_value=agf_parse_command_opts(argc, argv, "i:I:l:k:W:v:V:s:", &opt_args);
  if (exit_value==FATAL_COMMAND_OPTION_PARSE_ERROR) return exit_value;

  //parse the command line arguments:
  if (argc < 2) {
    printf("\n");
    printf("Syntax:	pdf_sim [-v var1] [-V var2] [-k k] [-W Wc] [-s nsamp] train output\n");
    printf("\n");
    printf("arguments:\n");
    printf("    train    file containing training vectors\n");
    printf("    output   binary output file\n");
    printf("\n");
    printf("options:\n");
    printf("    -s nsamp number of samples to generate [%d]\n", NBORD_DEFAULT);
    printf("    -v var1  upper filter variance bracket\n");
    printf("               --default is to use total variance of the data/n^(2/D)\n");
    printf("    -V var2  upper filter variance bracket\n");
    printf("               --default is to use total variance of the data\n");
    printf("    -k k     number of nearest neighbours to use in each estimate\n");
    printf("               --default is to use all of the data\n");
    printf("    -W Wc    objective total weight (default=%6.1f)\n", opt_args.W2);
    printf("\n");
    return INSUFFICIENT_COMMAND_ARGS;
  }

  diagfs=stderr;

  trainfile=argv[0];
  outfile=argv[1];

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

  std=NULL;
  ave=NULL;
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

  //fprintf(diagfs, "Normalizing the data.  ");
  //always normalize the data:
  norm_vec(train, nvar, ntrain, ave, std);

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

  //check range of k:
  if (opt_args.k <= opt_args.W2 || opt_args.k >= ntrain) {
    if (opt_args.k != -1) {
      fprintf(stderr, "Warning: parameter k=%d out of range.  Using all the training data\n", opt_args.k);
      opt_args.k=-1;
      exit_value=PARAMETER_OUT_OF_RANGE;
    }
    //fprintf(diagfs, "Using all the training data");
    //fprintf(diagfs, "\n");
  }

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

  ran_init();
  //begin the classification scheme:
  result=new real_a*[opt_args.n];
  result[0]=new real_a[opt_args.n*nvar];

  for (dim_ta j=0; j<nvar; j++) result[0][j]=rang()*scale;
  if (opt_args.k<0) {
    pold=agf_calc_pdf(train, nvar, ntrain, result[0],
                      opt_args.var, opt_args.W2, &diag_param);
  } else {
    pold=agf_calc_pdf(train, nvar, ntrain, result[0],
                      opt_args.var, opt_args.k, opt_args.W2, &diag_param);

    //calculate diagnostics:
    if (diag_param.f < min_f) min_f=diag_param.f;
    if (diag_param.f > max_f) max_f=diag_param.f;
    total_f+=diag_param.f;
  }
  //calculate diagnostics:
  if (diag_param.nd < min_nd) min_nd=diag_param.nd;
  if (diag_param.nd > max_nd) max_nd=diag_param.nd;
  total_nd+=diag_param.nd;

  if (diag_param.W < min_W) min_W=diag_param.W;
  if (diag_param.W > max_W) max_W=diag_param.W;
  total_W+=diag_param.W;

  ntest=1;
  for (nel_ta i=1; i<opt_args.n; i++) {
    result[i]=result[0]+nvar*i;
    do {
      for (dim_ta j=0; j<nvar; j++) result[i][j]=result[i-1][j]+rang()*scale;
      if (opt_args.k<0) {
        pnew=agf_calc_pdf(train, nvar, ntrain, result[i],
                      opt_args.var, opt_args.W2, &diag_param);
      } else {
        pnew=agf_calc_pdf(train, nvar, ntrain, result[i],
                      opt_args.var, opt_args.k, opt_args.W2, &diag_param);
        //calculate diagnostics:
        if (diag_param.f < min_f) min_f=diag_param.f;
        if (diag_param.f > max_f) max_f=diag_param.f;
        total_f+=diag_param.f;
      }
      //calculate diagnostics:
      if (diag_param.nd < min_nd) min_nd=diag_param.nd;
      if (diag_param.nd > max_nd) max_nd=diag_param.nd;
      total_nd+=diag_param.nd;

      if (diag_param.W < min_W) min_W=diag_param.W;
      if (diag_param.W > max_W) max_W=diag_param.W;
      total_W+=diag_param.W;

      rann=ranu();
      ntest++;
    } while (rann>pnew/pold);
  }

  for (nel_ta i=0; i<opt_args.n; i++) {
    for (dim_ta j=0; j<nvar; j++) {
      result[i][j]=result[i][j]*std[j]+ave[j];
      printf("%g ", result[i][j]);
    }
    printf("\n");
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
  fwrite(&nvar, sizeof(nel_ta), 1, fs);
  fwrite(result[0], sizeof(real_a), opt_args.n*nvar, fs);
  fclose(fs);

  //clean up:
  ran_end();
  delete [] result[0];
  delete [] result;
  delete [] train[0];
  delete [] train;

  if (std != NULL) delete std;
  if (ave != NULL) delete ave;

  return exit_value;

}


