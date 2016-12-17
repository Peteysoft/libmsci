#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>

#include "agf_lib.h"

using namespace std;
using namespace libagf;

//compares two sets of class and outputs statistics on them:
int main(int argc, char ** argv) {
  char *tfile;
  char *clsfile;
  char *confile;
  size_t slen;

  cls_ta *class1;
  cls_ta *class2;
  nel_ta n1, n2;

  real_a *con;

  int nhist=NCONHIST;

  flag_a bflag=0;
  flag_a Hflag=0;
  int exit_code=0;

  char c;

  while ((c = getopt(argc, argv, "Hbq:")) != -1) {
    switch (c) {
      case ('b'):
             bflag=1;
	     break;
      case ('q'):
             sscanf(optarg, "%d", &nhist);
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
    printf("\nCompares two sets of classes.  Prints out a set of statistics\n");
    printf("Optionally checks the accuracy of estimated conditional probabilities\n");
    printf("by comparing with calculated accuracies\n\n");
    printf("usage:  cls_comp_stats [-b] file1 file2\n\n");
    printf("where:\n");
    printf("  file1  = binary file containing first set of classes ('truth')\n");
    printf("  file2  = binary file(s) containing second set of classes ('retrieved')\n");
    printf("             .cls for classification results\n");
    printf("             .con for confidence ratings\n\n");
    printf("  -b     = short output (uncertainty coefficient and accuracy)\n");
    printf("  -H     = for short output, don't print header\n");
    printf("  -q     = number of histogram bins for evaluating confidence ratings [%d]\n", NCONHIST);
    printf("\n");
    return INSUFFICIENT_COMMAND_ARGS;
  }

  tfile=argv[0];
  slen=strlen(argv[1]);
  clsfile=new char[slen+5];
  strcpy(clsfile, argv[1]);
  strcat(clsfile, ".cls");
  confile=new char[slen+5];
  strcpy(confile, argv[1]);
  strcat(confile, ".con");

  class1=read_clsfile<cls_ta>(tfile, n1);
  if (n1 < 0) {
    fprintf(stderr, "Error reading input file: %s\n", tfile);
    return ALLOCATION_FAILURE;
  }
  if (class1 == NULL) {
    fprintf(stderr, "Unable to open file for reading: %s\n", tfile);
    return UNABLE_TO_OPEN_FILE_FOR_READING;
  }
  
  class2=read_clsfile<cls_ta>(clsfile, n2);
  if (n2 < 0) {
    fprintf(stderr, "Error reading input file: %s\n", clsfile);
    return ALLOCATION_FAILURE;
  }
  if (class2 == NULL) {
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
    if (con!=NULL) {
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

  return exit_code;

}

