#include <sys/timeb.h>
#include <math.h>

//#include "nr.h"
//#include "nrutil.h"
#include <gsl/gsl_randist.h>
#include "randomize.h"

#include "agf_lib.h"
#include "sample_class2_obj.h"

using namespace std;
using namespace libpetey;

namespace libagf {

template <class real>
sample_class2_obj<real>::sample_class2_obj(real w) {
  real ds;

  double xdiff;
  double ydiff;

  rann=gsl_rng_alloc(agf_gsl_rng_type);
  gsl_rng_set(rann, seed_from_clock());

  //set the spine:
  set_spine();

  //pre-fit the spline (spine ;-)):

  //calculate the path:
  s=new double[ns];
  s[0]=0;
  for (int i=0; i<ns-1; i++) {
    xdiff=xs[i+1]-xs[i];
    ydiff=ys[i+1]-ys[i];
    ds=sqrt(xdiff*xdiff+ydiff*ydiff);
    s[i+1]=s[i]+ds;		//***this is approximate!
  }

  //fit the spline:
  xinterp=gsl_interp_alloc(gsl_interp_cspline, ns);
  gsl_interp_init(xinterp, s, xs, ns); 

  yinterp=gsl_interp_alloc(gsl_interp_cspline, ns);
  gsl_interp_init(yinterp, s, ys, ns);

  saccel=gsl_interp_accel_alloc();

  width=w;

}

//delete the "spine samples":
template <class real>
sample_class2_obj<real>::~sample_class2_obj() {
  delete [] xs;
  delete [] ys;
  delete [] s;

  gsl_interp_free(xinterp);
  gsl_interp_free(yinterp);

  gsl_interp_accel_free(saccel);

  gsl_rng_free(rann);

}

template <class real>
void sample_class2_obj<real>::set_spine() {
  ns=9;

  xs=new double[ns];
  ys=new double[ns];

  xs[0]=0.17;    ys[0]=0.79;
  xs[1]=0.36;    ys[1]=0.7;
  xs[2]=0.51;    ys[2]=0.84;
  xs[3]=0.7;     ys[3]=0.86;
  xs[4]=0.76;    ys[4]=0.68;
  xs[5]=0.77;    ys[5]=0.48;
  xs[6]=0.68;    ys[6]=0.32;
  xs[7]=0.46;    ys[7]=0.28;
  xs[8]=0.24;    ys[8]=0.26;

  ddx0=1;
  ddxf=1;
  ddy0=0;
  ddyf=0;

}

template <class real>
void sample_class2_obj<real>::sample(real &x, real &y) {
  real x1, y1;
  real d;
  real theta;
  real offset;

  //random distance along the spine:
  d=gsl_rng_uniform(rann)*s[ns-1];

  x1=gsl_interp_eval(xinterp, s, xs, d, saccel);
  y1=gsl_interp_eval(yinterp, s, ys, d, saccel);

  //offset the location a random distance in a random direction:
  /*
  offset=gsl_ran_gaussian(rann, width);
  theta=gsl_rng_uniform(rann)*M_PI*2;

  x=x1+offset*cos(theta);
  y=y1+offset*sin(theta);
  */

  x=x1+gsl_ran_gaussian(rann, width);
  y=y1+gsl_ran_gaussian(rann, width);

}

template <class real>
real sample_class2_obj<real>::pdf(real x, real y, real *dpdx, real *dpdy, real *Lap, nel_ta nsamples) {
  real ds;		//step size
  real x1, y1;		//point along spine
  real ss;		//distance along spine
  real dx, dy;		//for calculating distance, below
  real d2;		//squared distance from point along spine
  real p1, p2;		//local probability
  real p;		//integrated probability
  real w2=width*width;	//squared width

  //mid-point method integration along the length of the spine:
  ds=s[ns-1]/nsamples;

  dx=x-xs[0];
  dy=y-ys[0];
  d2=dx*dx+dy*dy;
  p1=exp(-d2/2/w2)/2/M_PI/w2;
  p=p1/2;

  if (dpdx!=NULL && dpdy!=NULL) {
    //derivatives:
    *dpdx=-p1*dx/w2/2;
    *dpdy=-p1*dy/w2/2;
  }
  if (Lap!=NULL) {
    //Laplacian:
    *Lap=p1*d2/w2/2;
  }

  for (nel_ta i=1; i<nsamples; i++) {
    ss=i*ds;
    x1=gsl_interp_eval(xinterp, s, xs, ss, saccel);
    y1=gsl_interp_eval(yinterp, s, ys, ss, saccel);
    dx=x-x1;
    dy=y-y1;
    d2=dx*dx+dy*dy;
    p1=exp(-d2/2/w2)/2/M_PI/w2;
    p+=p1;

    if (dpdx!=NULL && dpdy!=NULL) {
      //derivatives:
      *dpdx-=p1*dx/w2;
      *dpdy-=p1*dy/w2;
    }
    if (Lap!=NULL) {
      //Laplacian:
      *Lap+=p1*d2/w2;
    }
  }

  dx=x-xs[ns-1];
  dy=y-ys[ns-1];
  d2=dx*dx+dy*dy;
  p1=exp(-d2/2/w2)/2/M_PI/w2;

  p+=p1/2;

  if (dpdx!=NULL && dpdy!=NULL) {
    //derivatives:
    *dpdx-=p1*dx/w2/2;
    *dpdy-=p1*dy/w2/2;
    *dpdx/=nsamples;
    *dpdy/=nsamples;
  }
  if (Lap!=NULL) {
    //Laplacian:
    *Lap+=p1*d2/2/M_PI/w2;
    *Lap=(*Lap-p)/nsamples;
  }

  return p/nsamples;

}
  
template <class real>
real sample_class2_obj<real>::pdf(real x, real y, nel_ta nsamples) {
  return pdf(x, y, NULL, NULL, NULL, nsamples);
}

template <class real>
real sample_class2_obj<real>::pdf(real x, real y, real &dpdx, real &dpdy, nel_ta nsamples) {
  return pdf(x, y, &dpdx, &dpdy, NULL, nsamples);
}

template <class real>
real sample_class2_obj<real>::pdf(real x, real y, real &Lap, nel_ta nsamples) {
  return pdf(x, y, NULL, NULL, &Lap, nsamples);
}

template class sample_class2_obj<float>;
template class sample_class2_obj<double>;

} //end namespace libagf

