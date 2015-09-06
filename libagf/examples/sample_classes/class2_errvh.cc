#include <math.h>
#include <stdio.h>

#include "parse_command_opts.h"

#include "sample_class2_obj.h"

using namespace libagf;
using namespace libpetey;

typedef float real;

int main(int argc, char **argv) {

  int n;		//number of samples
  real *x, *y;		//location of samples
  sample_class2_obj<real> pdf;	//object for generating samples and calculating
  				//actual pdf
  real x0, y0;		//test point
  real h1, h2;		//range for filter variance
  int ntest;
  real *h;		//array of filter variances (bandwidths)
  real *f;		//estimated pdf
  real p;		//actual pdf
  real *d2;		//distances squared

  void *opt_arg[10];
  int flag[10];

  opt_arg[0]=&h1;
  opt_arg[1]=&h2;
  opt_arg[2]=&ntest;

  h1=0.001;
  h2=1.;
  ntest=20;

  int err=parse_command_opts(argc, argv, "vVqg", "%g%g%d%", opt_arg, flag, 1);

  n=atoi(argv[1]);
  x0=atof(argv[2]);
  y0=atof(argv[3]);

  x=new real[n];
  y=new real[n];

  for (int i=0; i<n; i++) {
    pdf.sample(x[i], y[i]);
  }

  h=new real[ntest];
  if (flag[3]) {
    for (int i=0; i<ntest; i++) h[i]=h1*pow(h2/h1, 1.*i/(ntest-1));
  } else {
    for (int i=0; i<ntest; i++) h[i]=h1+i*(h2-h1)/(ntest-1);
  }

  d2=new real[n];
  for (int i=0; i<n; i++) {
    real xdiff=x0-x[i];
    real ydiff=y0-y[i];
    d2[i]=xdiff*xdiff+ydiff*ydiff;
  }

  f=new real[ntest];
  for (int i=0; i<ntest; i++) {
    f[i]=0;
    for (int j=0; j<n; j++) f[i]+=exp(-d2[j]/h[i]/2);
    f[i]/=h[i]*M_PI*2*n;
  }
  p=pdf.pdf(x0, y0);

  printf("p=%g\n", p);
  for (int i=0; i<ntest; i++) {
    printf("%g %g %g %g\n", h[i], f[i], p-f[i], (p-f[i])*(p-f[i]));
  }

  delete [] d2;
  delete [] h;
  delete [] f;
  delete [] x;
  delete [] y;

  return 0;

}
  
