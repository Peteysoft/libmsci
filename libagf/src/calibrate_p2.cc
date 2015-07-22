#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>

#include <gsl/gsl_linalg.h>

#include "full_util.h"

#include "agf_lib.h"

using namespace std;
using namespace libagf;
using namespace libpetey;

//reads in a multi-borders control file and a set of training data and uses 
//the training data to determine the mapping from the binary classifications
//to the multi-class classifications

//setting this as a stand-alone utility until I clean up the multi_borders/
//classify_m complex enough to figure out where to fit it in...
int main(int argc, char ** argv) {
  FILE *fs;

  char *clsfile;
  char *confile;
  size_t slen;

  cls_ta *class1;		//true classes
  real_a **x;
  dim_ta nvar;
  nel_ta n1, n2;

  multiclass_hier<real_a, cls_ta> *classifier;

  int exit_code=0;
  char c;

  while ((c = getopt(argc, argv, "r:")) != -1) {
    switch (c) {
      case ('r'):
             fprintf(stderr, "optimal_r0: -r not implemented yet!\n");
	     break;
      case ('?'):
             fprintf(stderr, "Unknown option: -%c -- ignored\n", optopt);
	     exit_code=COMMAND_OPTION_PARSE_ERROR;
	     break;
      default:
	     fprintf(stderr, "Error parsing command line\n");
	     return FATAL_COMMAND_OPTION_PARSE_ERROR;
	     break;
    }
  }

  argc-=optind;
  argv+=optind;

  if (argc < 1) {
    printf("\nDevelops a mapping from a set of binary classifiers to the multi-class classifier\n\n");
    printf("usage:  optimize_p2 control train mapfile\n\n");
    printf("where:\n");
    printf("  control   = input control file\n");
    printf("  train     = binary files containing training data\n");
    printf("             .vec for vector data\n");
    printf("             .cls for class labels\n\n");
    printf("  mapfile   = name of control file containing the mapping\n");
    printf("\n");
    return INSUFFICIENT_COMMAND_ARGS;
  }

  classifier=new multiclass_hier<real_a, cls_ta>(argv[0]);

  slen=strlen(argv[1]);
  clsfile=new char[slen+5];
  strcpy(clsfile, argv[1]);
  strcat(clsfile, ".cls");
  confile=new char[slen+5];
  strcpy(confile, argv[1]);
  strcat(confile, ".vec");

  //read in the classes:
  class1=read_clsfile(clsfile, n1);
  if (n2 < 0) {
    fprintf(stderr, "Error reading input file: %s\n", clsfile);
    return ALLOCATION_FAILURE;
  }
  if (class1 == NULL) {
    fprintf(stderr, "Unable to open file for reading: %s\n", clsfile);
    return UNABLE_TO_OPEN_FILE_FOR_READING;
  }

  //read in vector data:
  x=read_vecfile(confile, n2, nvar);
  if (n2 < 0) {
    fprintf(stderr, "Error reading input file: %s\n", confile);
    return ALLOCATION_FAILURE;
  }
  if (x == NULL) {
    fprintf(stderr, "Unable to open file for reading: %s\n", confile);
    return UNABLE_TO_OPEN_FILE_FOR_READING;
  }
  if (n1 != n2) {
    fprintf(stderr, "Data elements in files, %s and %s, do not agree: %d vs. %d\n", 
		    clsfile, confile, n1, n2);
    return SAMPLE_COUNT_MISMATCH;
  }

  classifier->train_map(x, class1, n1);

  fs=fopen(argv[2], "w");
  classifier->print(fs);
  fclose(fs);

  delete [] class1;

  delete_matrix(x);

  return exit_code;

}

