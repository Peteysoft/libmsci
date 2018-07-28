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
  //maximum tolerated residual before throwing an error:
  float eps=1e-7;
  double eps2=1e-15;
  //maximum residual for test_sparse_mem:
  float tol=1e-5;

  //option flags
  int dflag=0;
  int gflag=0;
  int hflag=0;
  int errcode=0;	//returned error
  //type of test:
  int type=0;
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
    printf("  usage: test_sparse [-d] [-G] [-h] type [ntrial [maxsize]]\n");
    printf("\n    where:\n");
    printf("  type       = type of test:\n");
    printf("                 0 = sparse arithmetic\n");
    printf("                 1 = sparse memory management\n");
    printf("  ntrial     = number of trials [%d]\n", ntrial);
    printf("  maxsize    = maximum size of system [%d]\n", maxsize);
    printf("  -d         = doube precision\n");
    printf("  -G         = use Gaussian deviates\n");
    printf("  -h         = print this help screen\n");
    return 0;
  }

  type=atoi(argv[0]);
  if (argc>1) {
    if (argc>2) {
      maxsize=atoi(argv[2]);
    }
    ntrial=atoi(argv[1]);
  }

  ran_init();

  for (int i=0; i<ntrial; i++) {
    m=ranu()*(maxsize-minsize+1)+minsize;
    n=ranu()*(maxsize-minsize+1)+minsize;
    p=ranu()*(maxsize-minsize+1)+minsize;
    sparsity=ranu()*(maxsparse-minsparse)+minsparse;
    switch (type) {
      case (0):
        printf("%dx%dx%d; %g sparsity\n", m, n, p, sparsity);
        if (dflag) {
          res2=test_sparse_arithmetic(m, n, p, (double) sparsity, eps2, gflag);
          printf("%lg\n", res2);
        } else {
          res=test_sparse_arithmetic(m, n, p, sparsity, eps, gflag);
          printf("%g\n", res);
        }
	break;
      case(1):
        printf("%dx%d; %g sparsity\n", m, n, sparsity);
        if (dflag) {
          errcode=test_sparse_mem(m, n, (double) sparsity, (double) tol, gflag);
          printf("%d\n", errcode);
        } else {
          errcode=test_sparse_mem(m, n, sparsity, tol, gflag);
          printf("%d\n", errcode);
        }
	break;
    }
  }

  ran_end();

  return errcode;

}

