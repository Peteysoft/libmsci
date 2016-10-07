#include <stdio.h>
#include <getopt.h>
#include <stdlib.h>

#include "sparse_test.h"
#include "error_codes.h"
#include "randomize.h"

using namespace libpetey;
using namespace libsparse;

int main(int argc, char **argv) {
  int m, n, p;
  float sparsity;
  float res;
  float res2;
  int ntrial=1;
  int maxsize=20;
  int minsize=1;
  float minsparse=0.;
  float maxsparse=0.9;
  float eps=1e-7;
  double eps2=1e-15;
  int dflag=0;
  int gflag=0;
  int hflag=0;
  int errcode=0;
  char c;

  while ((c = getopt(argc, argv, "dGh")) != -1) {
    switch (c) {
      case ('d'):
             dflag=1;
             break;
      case ('G'):
             gflag=1;
             break;
      case ('h'):
             hflag=1;
             break;
      case ('?'):
             fprintf(stderr, "Unknown option: -%c -- ignored\n", optopt);
             errcode=COMMAND_OPTION_PARSE_ERROR;
             break;
      default:
             fprintf(stderr, "Error parsing command line\n");
             exit(FATAL_COMMAND_OPTION_PARSE_ERROR);
    }
  }

  argc-=optind;
  argv+=optind;

  if (hflag) {
    printf("  usage: test_sparse [-d] [-G] [-h] [ntrial [maxsize]]\n");
    printf("\n    where:\n");
    printf("  ntrial     = number of trials [%d]\n", ntrial);
    printf("  maxsize    = maximum size of system [%d]\n", maxsize);
    printf("  -d         = doube precision\n");
    printf("  -G         = use Gaussian deviates\n");
    printf("  -h         = print this help screen\n");
    return 0;
  }

  if (argc>0) {
    if (argc>1) {
      maxsize=atoi(argv[1]);
    }
    ntrial=atoi(argv[0]);
  }

  ran_init();

  for (int i=0; i<ntrial; i++) {
    m=ranu()*(maxsize-minsize)+minsize;
    n=ranu()*(maxsize-minsize)+minsize;
    p=ranu()*(maxsize-minsize)+minsize;
    sparsity=ranu()*(maxsparse-minsparse)+minsparse;
    printf("%dx%dx%d; %g sparsity\n", m, n, p, sparsity);
    if (dflag) {
      res2=test_sparse_arithmetic(m, n, p, (double) sparsity, eps2, gflag);
      printf("%lg\n", res2);
    } else {
      res=test_sparse_arithmetic(m, n, p, sparsity, eps, gflag);
      printf("%g\n", res);
    }
  }

  ran_end();

  return errcode;

}

