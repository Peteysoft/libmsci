#include <assert.h>

#include "error_codes.h"
#include "full_util.h"
#include "peteys_tmpl_lib.h"

#include "agf_lib.h"

#include "../../libpetey/nnls.h"

using namespace libagf;
using namespace libpetey;

int main(int argc, char **argv) {
  FILE *fs;

  //as read from files:
  nel_ta n1, n2;		//number of samples
  real_a *g;			//decision values
  cls_ta *cls;			//"true" class values

  //as passed to nnls:
  int32_t n;			//number of samples
  double *f;			//decision values
  double *y;			//class values
  double **at;			//matrix mapping delta to p
  double *delta;		//difference between each probability

  //incidental values:
  double rnorm;
  double *work;
  double *zz;
  int32_t *index;
  int32_t mode;

  //more variables:
  long *sind;
  double p;
  double sum;

  if (argc != 3) {
    printf("syntax:     calibrate_probabilities class decision\n");
    exit(0);
  }

  cls=read_clsfile<cls_ta>(argv[1], n1);
  g=read_datfile<real_a>(argv[2], n2);

  assert(n1==n2);

  n=n1;
  at=allocate_matrix<double>(n, n);
  f=new double[n];
  y=new double[n];
  delta=new double[n];

  sind=heapsort(g, n);

  sum=0;
  for (int32_t i=0; i<n; i++) {
    f[i]=g[sind[i]];
    y[i]=1.*cls[sind[i]];
    //sum+=cls[sind[i]];
    //y[i]=sum;
    for (int32_t j=0; j<=i; j++) {
      at[j][i]=1;
      //at[j][i]=n-j;
    }
    for (int32_t j=i+1; j<n; j++) at[j][i]=0;
  }

  work=new double[n];
  zz=new double[n];
  index=new int32_t[n];

  FORTRAN_FUNC(nnls) (at[0], &n, &n, &n, y, delta, 
		  &rnorm, work, zz, index, &mode);

  /*
  p=0;
  for (int32_t i=0; i<n; i++) {
    p+=delta[i];
    printf("%lg %lg\n", f[i], 2*p-1);
  }
  printf("\n");
  */

  p=delta[0];
  printf("%lg %lg\n", f[0], 2*p-1);
  for (int32_t i=1; i<n-1; i++) {
    p+=delta[i];
    //printf("%d %lg %lg %lg %lg\n", cls[sind[i]], y[i], f[i], delta[i], p);
    if (delta[i]!=0) {
      if (delta[i-1]==0) printf("%lg %lg\n", f[i-1], 2*(p-delta[i])-1);
      printf("%lg %lg\n", f[i], 2*p-1);
    }
  }
  p+=delta[n-1];
  printf("%lg %lg\n", f[n-1], 2*p-1);

  //delete unneeded variables:
  delete [] cls;
  delete [] g;
  delete [] sind;

  delete [] at[0];
  delete [] at;
  delete [] f;
  delete [] y;
  delete [] work;
  delete [] zz;
  delete [] index;

  switch (mode) {
    case (1): return 0;
    case (2): return DIMENSION_MISMATCH;
    case (3): return MAX_ITER_EXCEEDED;
    default: return OTHER_ERROR;
  }

}
