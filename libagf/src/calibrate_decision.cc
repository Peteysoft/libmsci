
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>

#include "full_util.h"
#include "linked.h"
#include "agf_lib.h"

using namespace std;
using namespace libagf;
using namespace libpetey;

int main(int argc, char *argv[]) {
  char *confile;	//back to the old way...
  char *classfile=NULL;
  char *truthfile;	//optional model data to pass to external command
  char *cpcommand;	//copy input borders files to output
  int cpnormfile=0;	//copy the normalization file?

  agf_command_opts opt_args;

  cls_ta *truth;
  real_a *con;		//training vectors
  cls_ta *cls;		//training data classes
  real_a *r;		//decision values

  nel_ta ntrain;	//number of training data points
  nel_ta n1;		//number of class labels
  nel_ta n2;

  real_a *param;	//derived coefficients

  real_a r0;		//threshold probability

  int exit_code=0;
  int err_code;

  //set defaults:
  opt_args.nt=NCONHIST;
  opt_args.Qtype=CALIBRATION_ORDER;
  opt_args.algtype=0;
  opt_args.rthresh=0;

  //normalization options: -n -S -a
  //supernewton iteration: -h -i -I
  exit_code=agf_parse_command_opts(argc, argv, "c:Q:r:", &opt_args);
  if (exit_code==FATAL_COMMAND_OPTION_PARSE_ERROR) return exit_code;

  if (argc < 3) {
    printf("\n");
    printf("purpose:  calibrates a binary decision function\n");
    printf("\n");
    printf("syntax:   calibrate_decision [-Q order] truth results\\\n");
    printf("                   [-A [-M]]  border-in train border-out \\\n");
    printf("\n");
    printf("arguments:\n");
    printf("  truth       binary file containing true class values\n");
    printf("  result      results from a classification run:\n");
    printf("                .cls for classes\n");
    printf("                .con for \"confidence ratings\"\n");
    printf("\n");
    printf("options:\n");
    printf("  -c          method: [0]\n");
    printf("                0 = free fit\n");
    printf("                1 = fix threshold (r0) based on -r\n");
    printf("                2 = optimize accuracy\n");
    printf("                3 = optimize uncertainty coefficient\n");
    printf("                4 = optimize correlation coefficient\n");
    printf("  -Q order    order of fitted polynomial [%d]\n", CALIBRATION_ORDER);
    printf("  -r r0       threshold probability\n");
    printf("\n");
    return INSUFFICIENT_COMMAND_ARGS;
  }

  truthfile=argv[0];
  classfile=new char[strlen(argv[1])+1];
  sprintf(classfile, "%s.cls", argv[1]);
  confile=new char[strlen(argv[1])+5];
  sprintf(confile, "%s.con", argv[1]);

  //read in the data:
  truth=read_clsfile(truthfile, ntrain);
  if (truth==NULL) exit(FILE_READ_ERROR);
  cls=read_clsfile(classfile, n1);
  if (cls==NULL) exit(FILE_READ_ERROR);
  con=read_datfile(confile, n2);
  if (con==NULL) exit(FILE_READ_ERROR);

  if (ntrain != n1 || ntrain != n2) exit(SAMPLE_COUNT_MISMATCH);

  //convert back to decision values:
  r=new real_a[ntrain];
  for (int i=0; i<ntrain; i++) {
    if (cls[i]==0) {
      r[i]=-con[i];
    } else if (cls[i]==1) {
      r[i]=con[i];
    } else {
      fprintf(stderr, "calibrate_decision: only valid for binary classifiers\n");
      exit(PARAMETER_OUT_OF_RANGE);
    }
  }

  //we might add more calibration methods in later but for now...
  param=calibrate_decision(truth, r, ntrain, opt_args.Qtype);

  //print out the fitted coefficients:
  printf_matrix(stdout, param, 1, opt_args.Qtype+1);

  //delete character strings containing file names:
  delete [] confile;
  delete [] classfile;
  delete [] truthfile;

  //delete integer and real_aing point arrays:
  delete[] truth;
  delete [] cls;
  delete [] con;
  delete [] param;

  return 0;

}


