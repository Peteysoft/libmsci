#include <math.h>

int bv_field(double t, float *x, float *v, void *param) {
  float b, x1, y1, r;
  if (((int) t) % 2 == 0) b=-1; else b=1;
  x1=x[0]-b;
  y1=x[1];
  r=sqrt(x1*x1+y1*y1);
  v[0]=y1/r;
  v[1]=-x1/r;

  return 0;
}

