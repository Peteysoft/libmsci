#include <stdint.h>
#include <stdio.h>
#include <math.h>

#include "error_codes.h"
#include "randomize.h"
#include "bit_array.h"
#include "peteys_tmpl_lib.h"

using namespace libpetey;

int main(int argc, char **argv) {
  int n;
  long *trial;
  bit_array *tobits;
  int **coding_matrix;
  int nfilled;
  int nperm;

  ran_init();

  n=atoi(argv[1]);
  nperm=pow(2, n);

  coding_matrix=new int *[n];
  coding_matrix[0]=new int[n*n];
  for (int i=1; i<n; i++) coding_matrix[i]=coding_matrix[0]+n*i;

  if (n>8*sizeof(long)-1) {
    fprintf(stderr, "Size must be %d or less\n", 8*sizeof(long)-1);
    exit(PARAMETER_OUT_OF_RANGE);
  }

  trial=randomize(nperm);

  nfilled=0;
  for (int i=0; i<nperm; i++) {
    int dprod;
    //printf("%d\n", trial[i]);
    tobits=new bit_array((word *) (trial+i), (sizeof(long)+1)/sizeof(word), n);
    for (long j=0; j<n; j++) {
      coding_matrix[nfilled][j]=2*(long) (*tobits)[j]-1;
      //printf("%d ", (int32_t) (*tobits)[j]);
    }
    //printf("\n");
    //for (int j=0; j<n; j++) printf("%3d", coding_matrix[nfilled][j]);
    //printf("\n");
    dprod=0;
    for (int j=0; j<nfilled; j++) {
      dprod=0;
      for (int k=0; k<n; k++) {
        dprod+=coding_matrix[nfilled][k]*coding_matrix[j][k];
      }
      if (dprod!=0) break;
    }
    if (dprod==0) {
      for (int j=0; j<n; j++) printf("%3d", coding_matrix[nfilled][j]);
      printf("\n");
      nfilled++;
    }
    if (nfilled>=n) break;
    delete tobits;
  }

  delete [] trial;
  delete [] coding_matrix[0];
  delete [] coding_matrix;

  ran_end();
}

