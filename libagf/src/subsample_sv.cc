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
// For directly sub-sampling LIBSVM support vectors.
//
// Seems to work OK for binary classifiers but not at all for multi-class
// classifiers.
//

#include <math.h>
#include <string.h>
#include <stdio.h>

#include <gsl/gsl_linalg.h>

#include "randomize.h"
#include "peteys_tmpl_lib.h"
#include "full_util.h"
#include "roots_mins.h"
#include "agf_lib.h"

using namespace std;
using namespace libagf;
using namespace libpetey;

void ffromk(real_a k, void *param, real_a *f) {
  void **p2=(void **) param;
  cls_ta ncls=*(cls_ta *) p2[0];
  nel_ta *cind=(nel_ta *) p2[1];
  real_a C=pow(*(nel_ta *) p2[2], k);
  real_a f0=*(real_a *) p2[3];
  real_a sum=0;


  for (cls_ta i=0; i<ncls; i++) {
    nel_ta ni=cind[i+1]-cind[i];
    sum+=C*pow(ni, -k)*ni;
  }
  //printf("C=%g; sum=%g; nt=%d\n", C, sum, cind[ncls]);

  *f=sum/cind[ncls]-f0;
}

int main(int argc, char *argv[]) {
  FILE *diagfs;
  FILE *fs;
  svm_multi<real_a, cls_ta> *svmmod;

  int exit_value;

  agf_command_opts opt_args;

  exit_value=0;

  //parse the command line arguments:
  //exit_value=agf_parse_command_opts(argc, argv, "c:d:f:AMRz", &opt_args);
  exit_value=agf_parse_command_opts(argc, argv, "d:f:CRz", &opt_args);
  if (exit_value==FATAL_COMMAND_OPTION_PARSE_ERROR) return exit_value;

  //print help page if there are not enough mandatory arguments:
  if (argc<2) {
    printf("\n");
    printf("Syntax:	subsample_sv [-d ndiv] [-f frac] \n");
    printf("                          input output test\n");
    printf("\n");
    printf("arguments:\n");
    printf("  input        name LIBSVM model file\n");
    printf("  output       name sub-sampled LIBSVM model file\n");
    printf("\n");
    //printf("file options:\n");
    //printf("  -A           operate on ASCII files\n");
    //printf("  -M           specifies LIBSVM format\n");
    //printf("\n");
    printf("  -C           keep relative class numbers constant\n");
    printf("  -z           randomly permute data\n");
    printf("  -R           data separation works by random selection rather than\n");
    printf("                 permutation (output files are not of definite size)\n");
    printf("  -d ndiv      number of separate output files\n");
    printf("  -f frac      separate into test and training (over-rides -d)\n");
    printf("\n");
    printf("\n");
    return INSUFFICIENT_COMMAND_ARGS;
  }

  diagfs=stderr;
  ran_init();

  svmmod=new svm_multi<real_a, cls_ta>(argv[0]);

  svmmod->subsample(opt_args.ftest, opt_args.Cflag);

  fs=fopen(argv[1], "w");

  svmmod->save(fs);

  ran_end();

  fclose(fs);

  return exit_value;

}


