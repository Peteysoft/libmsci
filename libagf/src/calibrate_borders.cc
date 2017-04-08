#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>

#include "randomize.h"
#include "agf_lib.h"

using namespace std;
using namespace libagf;
using namespace libpetey;

//compares two sets of class and outputs statistics on them:
int main(int argc, char ** argv) {
  FILE *fs;
  classifier_obj *classifier;
  int nhist=NCONHIST;
  int order=CALIBRATION_ORDER;
  real_a **train;
  cls_ta *cls;
  nel_ta ntrain;
  int exit_code=0;

  char c;

  while ((c = getopt(argc, argv, "O:q:")) != -1) {
    switch (c) {
      case ('O'):
             bflag=1;
	     break;
      case ('q'):
             sscanf(optarg, "%d", &nhist);
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
    printf("\nCalibrates a borders classifier.\n");
    printf("\n");
    printf("usage:  calibrate_borders classifier train\n\n");
    printf("where:\n");
    printf("  classifier = borders classifier\n");
    printf("  train      = training data\n");
    printf("  -q         = number of histogram bins for fitting [%d]\n", NCONHIST);
    printf("  -O         = order of fitting [%d]\n", CALIBRATION_ORDER);
    printf("\n");
    return INSUFFICIENT_COMMAND_ARGS;
  }

  fs=fopen(argv[0], "r");
  classifier=new multiclass_hier<real_a, cls_ta>();
  if (classifier->load(argv[0])==PARAMETER_OUT_OF_RANGE) {
    delete classifier;
    try {
      classifier=new multiclass_hier<real_a, cls_ta>(fs);
      fclose (fs);
    } catch (int) {
      delete classifier;
      fclose (fs);
      classifier=new borders_calibrated<real_a, cls_ta>();
      classifier->init(argv[0]);
    }
  } else {
    fclose(fs);
  }

  exit_code=agf_read_train(argv[1], train, cls, ntrain, nvar);

  if (exit_code!=0) exit(exit_code);

  classifier->calibrate(

  return exit_code;

}

