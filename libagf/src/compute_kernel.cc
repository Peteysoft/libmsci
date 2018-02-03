#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#include "getopt.h"

#include "full_util.h"
#include "agf_lib.h"

using namespace std;
using namespace libagf;
using namespace libpetey;

real_a kappa_chi(real_a x) {
  return 1/cosh(x);
}

real_a kappa_int(real_a x) {
  return 2/(1+4*x*x)/M_PI;
}

int main(int argc, char ** argv) {

  FILE *fs;
  char *testfile;
  char *testoutfile;

  real_a (* kfunc) (real_a *, real_a *, dim_ta, void *);

  real_a (* kappa) (real_a);

  nel_ta ntrain;		//number of training samples
  dim_ta nvar;			//number of features

  real_a **vec;
  cls_ta *cls;
  real_a min=0;

  agf_command_opts opts;

  int errcode=0;

  agf_parse_command_opts(argc, argv, "Q:q:w:", &opts);

//  long *clind;
//  long ncl;
//
  
  if (argc < 2 ) {
    printf("Pre-compute additive kernels for LIBSVM\n");
    printf("\n");
    printf("usage: compute_kernel [-Q type] train [trainout] [test testout]\n\n");
    printf("    where:\n");
    printf("-Q       = type of kernel:\n");
    printf("             0: Hellinger's\n");
    printf("             1: Chi-squared\n");
    printf("             2: intersection\n");
    printf("train    = original training data\n");
    printf("trainout = transformed training data\n");
    printf("test     = original test data\n");
    printf("testout  = transformed test data\n");
    exit(1);
  }

  //read in the training data:
  fs=fopen(argv[0], "r");
  ntrain=read_svm(fs, vec, cls, nvar);
  if (nvar==-1 || ntrain==-1) {
    fprintf(stderr, "compute_kernel: error reading file: %s\n", argv[0]);
    exit(FILE_READ_ERROR);
  }
  if (vec==NULL) {
    fprintf(stderr, "compute_kernel: could not open file, %s, for reading\n", argv[0]);
    exit(UNABLE_TO_OPEN_FILE_FOR_READING);
  }
  fclose(fs);

  switch (opts.Qtype) {
    case (0): kfunc=&Hellingers_kernel<real_a>;
	      break;
    case (1): kfunc=&Chi_squared_kernel<real_a>;
	      break;
    case (2): kfunc=&intersection_kernel<real_a>;
	      break;
    case (3): kappa=&kappa_chi;
	      break;
    case (4): kappa=&kappa_int;
	      break;
    default: kfunc=&Hellingers_kernel<real_a>;
	      break;
  }

  //some of the kernels don't work with negative training data:
  if (opts.Qtype==0 || opts.Qtype == 3 || opts.Qtype == 4 || opts.Qtype==6) {
    for (int i=0; i<ntrain; i++) {
      for (dim_ta j=0; j<nvar; j++) {
        if (vec[i][j]<min) min=vec[i][j];
      }
    }
    for (int i=0; i<ntrain; i++) {
      for (dim_ta j=0; j<nvar; j++) {
        vec[i][j]-=1.0001*min;			//stupid hacks...
      }
    }
  }

  printf("argc=%d\n", argc);


  if (argc!=3) {

    printf("argc!=3\n");
    fs=fopen(argv[1], "w");
    if (opts.Qtype<3) {
      for (int i=0; i<ntrain; i++) {
        fprintf(fs, "%d", cls[i]);
        fprintf(fs, " 0:%d", i+1);
        for (int j=0; j<ntrain; j++) {
          fprintf(fs, " %d:%g", j+1, (*kfunc)(vec[i], vec[j], nvar, NULL));
        }
        fprintf(fs, "\n");
      }
    } else if (opts.Qtype==3 || opts.Qtype==4) {
      real_a f1, f2;
      dim_ta cnt;
      for (int i=0; i<ntrain; i++) {
        fprintf(fs, "%d", cls[i]);
        cnt=1;
        for (int j=0; j<nvar; j++) {
          f1=sqrt(vec[i][j]*opts.W1)*sqrt((*kappa)(vec[i][j]));
          fprintf(fs, " %d:%g", cnt, f1);
          cnt++;
          for (int k=1; k<opts.nt; k++) {
            f1=sqrt(vec[i][j]*opts.W1)*sqrt((*kappa)(opts.W1*(k+1)/2))*cos((k+1)*opts.W1*log(vec[i][j])/2);
            f2=sqrt(vec[i][j]*opts.W1)*sqrt((*kappa)(opts.W1*k/2))*sin(k*opts.W1*log(vec[i][j])/2);
            fprintf(fs, " %d:%g %d:%g", cnt, f1, cnt+1, f2);
	    cnt+=2;
          }
        }
        fprintf(fs, "\n");
      }
    } else if (opts.Qtype==5) {
      dim_ta cnt;
      for (int i=0; i<ntrain; i++) {
        fprintf(fs, "%d", cls[i]);
        cnt=1;
        for (int j=0; j<nvar; j++) {
          fprintf(fs, " %d:%g", cnt, vec[i][j]);
          cnt++;
        }
        for (int j=0; j<nvar; j++) {
          fprintf(fs, " %d:%g", cnt, vec[i][j]*vec[i][j]);
          cnt++;
        }
        for (int j=0; j<nvar; j++) {
          for (int k=j+1; k<nvar; k++) {
            fprintf(fs, " %d:%g", cnt, vec[i][j]*vec[i][k]);
            cnt++;
          }
        }
        fprintf(fs, "\n");
      }
    } else if (opts.Qtype==6) {
      for (int i=0; i<ntrain; i++) {
        fprintf(fs, "%d", cls[i]);
        for (int j=0; j<nvar; j++) {
          fprintf(fs, " %d:%g", j+1, sqrt(vec[i][j]));
        }
        fprintf(fs, "\n");
      }
    }
    fclose(fs);
    testfile=argv[2];
    testoutfile=argv[3];
  } else {
    testfile=argv[1];
    testoutfile=argv[2];
  }

  printf("%s %s\n", testfile, testoutfile);

  if (argc>2 && opts.Qtype < 3) {
    nel_ta ntest;			//number of test samples
    dim_ta nvar1;			//features in test data
    real_a **testvec;
    cls_ta *testcls;

    fs=fopen(testfile, "r");
    ntest=read_svm(fs, testvec, testcls, nvar1);
    if (nvar1==-1 || ntest==-1) {
      fprintf(stderr, "compute_kernel: error reading file: %s\n", argv[2]);
      exit(FILE_READ_ERROR);
    }
    if (testvec==NULL) {
      fprintf(stderr, "compute_kernel: could not open file, %s, for reading\n", argv[2]);
      exit(UNABLE_TO_OPEN_FILE_FOR_READING);
    }
    fclose(fs);
    if (nvar!=nvar1) {
      fprintf(stderr, "compute_kernel: dimensions in test and train do not match\n");
      exit(DIMENSION_MISMATCH);
    }

    //have to do the same stupid thing to the test data:
    if (opts.Qtype==0 || opts.Qtype == 3) {
      for (int i=0; i<ntest; i++) {
        for (dim_ta j=0; j<nvar; j++) {
          testvec[i][j]-=1.0001*min;			//stupid hacks...
        }
      }
    }
    fs=fopen(testoutfile, "w");
    for (int i=0; i<ntest; i++) {
      fprintf(fs, "%d", testcls[i]);
      fprintf(fs, " 0:%d", i+1);
      for (int j=0; j<ntrain; j++) {
        fprintf(fs, " %d:%g", j+1, (*kfunc)(testvec[i], vec[j], nvar, NULL));
      }
      fprintf(fs, "\n");
    }
    fclose(fs);

    delete [] testcls;
    delete [] testvec[0];
    delete [] testvec;

  }

  delete [] cls;
  delete [] vec[0];
  delete [] vec;

  return errcode;

}

