#include <math.h>
#include <assert.h>

#include "agf_util.h"

int main(int argc, char **argv) {
  float *data1;
  float *data2;
  long n1, n2;
  long n;
  float ave1, ave2;
  float std1, std2;
  float diff1, diff2;
  float cov;

  data1=(float *) read_clsfile(argv[1], n1);
  data2=(float *) read_clsfile(argv[2], n2);

  assert(n1==n2);

  ave1=0;
  ave2=0;
  for (long i=0; i<n1; i++) {
    if (finite(data1[i])) {
      ave1+=data1[i];
    } else {
      n1--;
    }
    if (finite(data2[i])) {
      ave2+=data2[i];
    } else {
      n2--;
    }
  }

  ave1/=n1;
  ave2/=n2;

  std1=0;
  std2=0;
  cov=0;
  n=0;
  for (long i=0; i<n1; i++) {
    if (finite(data1[i])) {
      diff1=data1[i]-ave1;
      std1+=diff1*diff1;
    }
    if (finite(data2[i])) {
      diff2=data2[i]-ave2;
      std2+=diff2*diff2;
    }
    if (finite(data1[i]) && finite(data2[i])) {
      n++;
      cov+=diff1*diff2;
    }
  }

  std1=sqrt(std1/(n1+1));
  std2=sqrt(std2/(n2+1));

  printf("%g\n", cov/std1/std2/(n-1));

}

