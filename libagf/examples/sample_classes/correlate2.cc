#include <math.h>
#include <assert.h>

#include "../agf_util.h"

int main(int argc, char **argv) {
  float *con1;
  float *con2;
  long n1, n2;

  float ave1, ave2;
  float var1, var2;
  float cov;
  float diff1, diff2;
  float r;

  con1=(float *) read_clsfile(argv[1], n1);
  con2=(float *) read_clsfile(argv[2], n2);
  assert(n1 == n2);

  ave1=0;
  ave2=0;
  for (long i=0; i<n1; i++) {
    ave1+=con1[i];
    ave2+=con2[i];
  }
  ave1/=n1;
  ave2/=n1;

  var1=0;
  var2=0;
  cov=0;
  for (long i=0; i<n1; i++) {
    diff1=con1[i]-ave1;
    diff2=con2[i]-ave2;
    var1+=diff1*diff1;
    var2+=diff2*diff2;
    cov+=diff1*diff2;
  }

  r=cov/sqrt(var1*var2);

  printf("r =%g\n", r);

  delete [] con1;
  delete [] con2;

}

