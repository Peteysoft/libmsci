#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>

#include "peteys_tmpl_lib.h"
#include "randomize.h"
#include "agf_lib.h"

using namespace std;
using namespace libagf;
using namespace libpetey;

//compares two sets of class and outputs statistics on them:
int main(int argc, char ** argv) {
  char *tfile;
  char *pfile;
  char *ofile;
  size_t slen;

  cls_ta *class1;
  cls_ta *class2;
  nel_ta n1, n2;

  real_a **p;

  flag_a bflag=0;
  flag_a Hflag=0;
  int exit_code=0;

  char c;

  while ((c = getopt(argc, argv, "HC")) != -1) {
    switch (c) {
      case ('C'):
             Cflag=1;
	     break;
      case ('H'):
	     Hflag=1;
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

  if (argc < 2) {
    printf("\nPoint-by-point validation of class conditional probabilities\n");
    printf("\n");
    printf("usage:  validate_probabilities class prob [output]\n\n");
    printf("where:\n");
    printf("  class  = binary file containing ture classes\n");
    printf("  prob   = ASCII file containing estimated probabilities\n");
    printf("  output = output file for plotting\n");
    printf("  -C     = no class data in ASCII file\n");
    printf("  -H     = no header in ASCII file\n");
    printf("\n");
    return INSUFFICIENT_COMMAND_ARGS;
  }

  ran_init();

  tfile=argv[0];
  pfile=argv[1];
  if (argc > 2) ofile=argv[2]; else ofile=NULL;

  class1=read_clsfile<cls_ta>(tfile, n1);
  if (n1 < 0) {
    fprintf(stderr, "Error reading input file: %s\n", tfile);
    return ALLOCATION_FAILURE;
  }
  if (class1 == NULL) {
    fprintf(stderr, "Unable to open file for reading: %s\n", tfile);
    return UNABLE_TO_OPEN_FILE_FOR_READING;
  }
 
  n2=read_lvq(pfile, p, class2, ncls, Hflag+2*Cflag); 
  if (n2 < 0) {
    fprintf(stderr, "Error reading input file: %s\n", tfile);
    return ALLOCATION_FAILURE;
  }
  if (p == NULL) {
    fprintf(stderr, "Unable to open file for reading: %s\n", clsfile);
    return UNABLE_TO_OPEN_FILE_FOR_READING;
  }
  if (n1 != n2) {
    fprintf(stderr, "Data elements in files, %s and %s, do not agree: %d vs. %d\n", 
		    tfile, clsfile, n1, n2);
    return SAMPLE_COUNT_MISMATCH;
  }

/*
  for (nel_ta i=0; i<n1; i++) {
    printf("%d %d\n", class1[i], class2[i]);
  }
*/

  if (bflag == 1) {
    class_eval_basic(class1, class2, n1, stdout, Hflag);
    //printf("%f %f\n", uc, (real_a) nt/(real_a) n1);
  } else {
    class_eval(class1, class2, n1);

    //check accuracy of confidence ratings:
    con=read_datfile<real_a>(confile, n2);
    if (con!=NULL && nhist > 0) {
      if (n1!=n2) {
        fprintf(stderr, "cls_comp_stats: number of confidence ratings (%d) does not match the number of classes (%d)\n", n2, n1);
        exit(SAMPLE_COUNT_MISMATCH);
      }
      check_confidence(class1, class2, con, n1, nhist);

      delete [] con;
    }
  }

  delete [] class1;
  delete [] class2;

  delete [] clsfile;
  delete [] confile;

  ran_end();

  return exit_code;

}

