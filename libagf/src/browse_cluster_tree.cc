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
//
// My attempt at hierarchical cluster analysis. Allows you to browse the tree 
// and add clusters.
//

#include <math.h>
#include <string.h>
#include <stdio.h>

#include "agf_lib.h"

using namespace std;
using namespace libagf;

int main(int argc, char *argv[]) {
  char *trainfile;		//training data
  char *outfile;		//output classes
  FILE *fs;
  FILE *diagfs;

  nel_ta ntrain;		//number of training data points
  dim_ta nvar;			//number of variables
  dim_ta nvar1;

  real_a **train;		//training data vectors
  cls_ta *result;		//results of cluster analysis

  cluster_tree<real_a, cls_ta> ctree;

  agf_command_opts opt_args;

  int exit_value;

  exit_value=0;
  exit_value=agf_parse_command_opts(argc, argv, "nS:o:", &opt_args);
  if (exit_value==FATAL_COMMAND_OPTION_PARSE_ERROR) return exit_value;

  //parse the command line arguments:
  if (argc < 2) {
    printf("\n");
    printf("Syntax:   browse_cluster_tree [-o log] [-n] [-S nsv] train output\n");
    printf("\n");
    printf("arguments:\n");
    printf("    train    file containing training vectors\n");
    printf("    output   binary output file containing classes\n");
    printf("\n");
    printf("options:\n");
    printf("    -n       normalize the data\n");
    printf("    -S nsv   perform SVD: number of singular values to keep\n");
    printf("    -o log   ASCII file logging steps\n");
    printf("\n");
    exit(0);
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
  train=agf_get_features<real_a>(trainfile, &opt_args, ntrain, nvar, 1);

  fprintf(diagfs, "%d training vectors found in file: %s.vec\n", ntrain, trainfile);

  //build the dendrogram:
  exit_value=ctree.build_all(train, ntrain, nvar);

  if (exit_value!=0) exit(exit_value);

  //start the interactive classification scheme:
  result=new cls_ta[ntrain];
  if (opt_args.ofile!=NULL) {
    fs=fopen(opt_args.ofile, "w");
    ctree.browse(result, fs);
    fclose(fs);
  } else {
    ctree.browse(result);
  }

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


