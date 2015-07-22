#include <math.h>
#include <sys/timeb.h>

//#include "nr.h"

#include <gsl/gsl_randist.h>
#include "randomize.h"

#include "agf_lib.h"
#include "sample_class1_obj.h"

using namespace libpetey;

namespace libagf {

template <class real>
sample_class1_obj<real>::sample_class1_obj() {

  //initialize random number generator:
  rann=gsl_rng_alloc(agf_gsl_rng_type);
  gsl_rng_set(rann, seed_from_clock());

  //set default parameters:
  x0=0.4;
  y0=0.5;

  w1=0.13;
  w2=0.07;

  cosang=cos(M_PI/7);
  sinang=sin(M_PI/7);

}

template <class real>
sample_class1_obj<real>::~sample_class1_obj() {
  gsl_rng_free(rann);
}

//user can set parameters of class:
template <class real>
sample_class1_obj<real>::sample_class1_obj(real x0u,
					real y0u,
					real w1u, 
					real w2u, 
					real angle) {

  //initialize random number generator:
  rann=gsl_rng_alloc(agf_gsl_rng_type);
  gsl_rng_set(rann, seed_from_clock());

  //set user-defined parameters:
  x0=x0u;
  y0=y0u;

  w1=w1u;
  w2=w2u;

  cosang=cos(angle);
  sinang=sin(angle);

}

template <class real>
void sample_class1_obj<real>::sample(real &x, real &y) {
  real x1, y1;
  real d;
  real theta;
  real offset;

  x1=gsl_ran_gaussian(rann, w1);
  y1=gsl_ran_gaussian(rann, w2);

  x=cosang*x1-sinang*y1+x0;
  y=sinang*x1+cosang*y1+y0;

}

template <class real>
real sample_class1_obj<real>::pdf(real x, real y) {
  real x1, y1;
  real d2;
  real p;

  //transform to normalized distribution:
  x=x-x0;
  y=y-y0;

  x1= cosang*x + sinang*y;
  y1=-sinang*x + cosang*y;

  x1/=w1;
  y1/=w2;

  d2=x1*x1+y1*y1;

  p=exp(-d2/2)/M_PI/2/w1/w2;

  return p;

}


template <class real>
real sample_class1_obj<real>::pdf(real x, real y, real &dpdx, real &dpdy) {
  real x1, y1;
  real d2;
  real p;

  //transform to normalized distribution:
  x=x-x0;
  y=y-y0;

  x1= cosang*x + sinang*y;
  y1=-sinang*x + cosang*y;

  x1/=w1;
  y1/=w2;

  d2=x1*x1+y1*y1;

  p=exp(-d2/2)/M_PI/2/w1/w2;

  dpdx=-p*(cosang*x1/w1-sinang*y1/w2);
  dpdy=-p*(sinang*x1/w1+cosang*y1/w2);

  return p;

}

template class sample_class1_obj<float>;
template class sample_class1_obj<double>;

} //end namespace libagf

