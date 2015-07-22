#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include "agf_lib.h"

int main(int argc, char *argv[]) {
  char *vecfile;
  char *clsfile;

  float **vec;
  long *cls;
  long n, n1;
  long nvar;

  float *ave, *std;
  int normflag;

  if (argc < 2) {
    printf("Converts vector data in libAGF compatible binary format\n");
    printf("to LIBSVM compatible ASCII format.\n\n");
    printf("usage: vecfile2svm [-n] vecfile [clsfile]\n\n");
    printf("    where:\n");
    printf("-n      = option to normalise the data\n");
    printf("vecfile = binary file containing vector data\n");
    printf("clsfile = optional binary file containing class values\n");
    exit(FATAL_COMMAND_OPTION_PARSE_ERROR);
  }

  normflag=0;
  if (argv[1][0] == '-') {
    if (argv[1][1] == 'n') {
      normflag=1;
    }
    argc--;
    argv++;
  }

  vecfile=argv[1];
  vec=read_vecfile(vecfile, n, nvar);

  if (normflag) {
    float diff;

    ave=new float[nvar];
    std=new float[nvar];

    for (long i=0; i<nvar; i++) {
      ave[i]=0;
      std[i]=0;
      for (long j=0; j<n; j++) ave[i]+=vec[j][i];
      ave[i]/=n;
      for (long j=0; j<n; j++) {
        diff=vec[j][i]-ave[i];
        vec[j][i]=diff;
        std[i]+=diff*diff;
      }
      std[i]=sqrt(std[i]/(n-1));
      for (long j=0; j<n; j++) vec[j][i]/=std[i];
    }
  }

  if (argc > 2) {
    clsfile=argv[2];
    cls=read_clsfile(clsfile, n1);
    assert(n==n1);
    printf("%d\n", nvar);
    for (long i=0; i<n; i++) {
      for (long j=0; j<nvar; j++) printf("%f ", vec[i][j]);
      printf("%d\n", cls[i]);
    }
    delete [] cls;
  } else {
    printf("%d\n", nvar);
    for (long i=0; i<n; i++) {
      for (long j=0; j<nvar; j++) printf("%f ", vec[i][j]);
      printf("\n");
    }
  }

  delete [] vec[0];
  delete [] vec;

}

