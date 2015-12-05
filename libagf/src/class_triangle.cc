#include "agf_lib.h"

using namespace libagf;

template <class real, class cls_t>
real tprob(real **x, int n, int D, 
		int clind,
		real *x0,
		real W,
		real var[2],
		real *dcond, 			//gradient of cond. prob
		real *dprob) {			//grad. of total prob
  real d2[n];			//distances
  real **dwdx;			//gradients of weights
  real w[n];			//weights
  real varf;			//returned variance (bandwidth)
  real wplus, wminus;
  real dwplus[D], dwminus[D];
  real tw;			//total of the weights
  real p;			//density
  real norm;			//normalization coefficient
  dwdx=new real *[n];
  dwdx[0]=new real[n*D];
  for (int i=0; i<n; i++) {
    d2[i]=metric2(x[i], x0, D);
    dwdx[i]=dwdx[0]+i*D;
  }

  AGF_CALC_W_FUNC(d2, n, W, var, w, varf);

  agf_grad_w(x, D, x0, w, d2, n, varf, dwdx);

  wplus=0;
  wminus=0;
  for (int j=0; j<D; j++) {
    dwplus[j]=0;
    dwminus[j]=0;
    for (int i=0; i<clind; i++) dwminus[i][j]+=w[i];
    for (int i=0; i<n-clind; i++) dwplus[i][j]+=w[i];
  }

  for (int i=0; i<clind; i++) wminus+=w[i];
  for (int i=0; i<n-clind; i++) wplus+=w[i];

  norm=pow(sqrt(varf*M_PI*2), D);
  tw=wminus+wplus;
  p=(wminus+wplus)/norm/n;

  for (int j=0; j<D; j++) {
    dcond[j]-=dwminus[j];
    dcond[j]+=dwplus[j];
    dcond[j]/=tw;
    dprob[j]+=dwminus[j];
    dprob[j]+=dwplus[j];
    dprob[j]/=norm/n;
  }

  delete [] dwdx[0];
  delete [] dwdx;

  return p;
}

int main(int argc, char **argv) {
  int err;
  real_a **x;
  cls_ta *cls;
  int n;
  int D;
  real_a *d;
  int ncls;
  cluster_tree<real_a, cls_ta> dg;


  err=agf_read_train(argv[1], x, cls, n, D);
  d=class_triangle(x, cls, n, D);
  ncls=0;
  for (int i=0; i<n; i++) if (cls[i]>=ncls) ncls=cls[i]+1;

  dg.build_all(d, ncls, D);
  dg.print(stdout);

  return err;

}
