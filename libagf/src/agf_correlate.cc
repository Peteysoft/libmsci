#include <math.h>
#include <assert.h>

#include "agf_lib.h"

using namespace std;
using namespace libagf;

int main(int argc, char **argv) {
  real_a *con1;
  real_a *con2;
  nel_ta n1, n2;

  real_a ave1, ave2;
  real_a var1, var2;
  real_a cov;
  real_a diff1, diff2;
  real_a r;

  if (argc < 3) {
    printf("Usage: correlate file1 file2\n");
    printf("\nCorrelates real data from two binary files\n");
    exit(-1);
  }

  con1=read_datfile(argv[1], n1);
  con2=read_datfile(argv[2], n2);
  assert(n1 == n2);

  ave1=0;
  ave2=0;
  for (nel_ta i=0; i<n1; i++) {
    ave1+=con1[i];
    ave2+=con2[i];
  }
  ave1/=n1;
  ave2/=n1;

  var1=0;
  var2=0;
  cov=0;
  for (nel_ta i=0; i<n1; i++) {
    diff1=con1[i]-ave1;
    diff2=con2[i]-ave2;
    var1+=diff1*diff1;
    var2+=diff2*diff2;
    cov+=diff1*diff2;
  }

  r=cov/sqrt(var1*var2);

  printf("r =%7.3f\n", r);

  delete [] con1;
  delete [] con2;

}

