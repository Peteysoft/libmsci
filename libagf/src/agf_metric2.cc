#include <math.h>

#include "agf_metric2.h"

using namespace std;
using namespace libagf;

const enum_a libagf::global_nmetric=6;

real_a (* libagf::global_metric2) (real_a *, real_a *, dim_ta m) = &metric2<real_a>;

real_a (* libagf::global_metric2_pointer_list[global_nmetric]) (real_a *, real_a *, dim_ta m)
	={&metric2<real_a>, 
 	  &manhattan2<real_a>, 
	  &quad_metric2<real_a>, 
	  &max_metric2<real_a>,
	  &quad_metric4<real_a>,
	  &manhattan<real_a>};

char libagf::global_metric_type[200]="0=Cart.^2, 1=Manhattan^2, 2=quad.^2, 3=max. element^2, 4=quad.^4, 5=Manhattan";

//since we rarely stray from a Cartesian metric,
//we maintain the convention of returning the squared 
//metric, rather than the value itself
template <class real>
real libagf::metric2(real *v1, real *v2, dim_ta m)
{
  real d;
  real diff;

  d=0;
  for (dim_ta i=0; i<m; i++) {
    diff=v2[i]-v1[i];
    d+=diff*diff;
  }

  return d;
}

template <class real>
real libagf::manhattan2(real *v1, real *v2, dim_ta m)
{
  real d;
  real diff;

  d=0;
  for (dim_ta i=0; i<m; i++) {
    diff=v2[i]-v1[i];
    d+=fabs(diff);
  }

  return d*d;
}


template <class real>
real libagf::quad_metric2(real *v1, real *v2, dim_ta m)
{
  real d;
  real diff;

  d=0;
  for (dim_ta i=0; i<m; i++) {
    diff=v2[i]-v1[i];
    diff=diff*diff;
    d+=diff*diff;
  }

  return sqrt(d);
}

template <class real>
real libagf::max_metric2(real *v1, real *v2, dim_ta m)
{
  real d;
  real diff;

  d=fabs(v2[0]-v1[0]);
  for (dim_ta i=1; i<m; i++) {
    diff=fabs(v2[i]-v1[i]);
    if (diff > d) d=diff;
  }

  return d*d;
}

template <class real>
real libagf::quad_metric4(real *v1, real *v2, dim_ta m)
{
  real d;
  real diff;

  d=0;
  for (dim_ta i=0; i<m; i++) {
    diff=v2[i]-v1[i];
    diff=diff*diff;
    d+=diff*diff;
  }

  return d;
}

template <class real>
real libagf::manhattan(real *v1, real *v2, dim_ta m)
{
  real d;
  real diff;

  d=0;
  for (dim_ta i=0; i<m; i++) {
    diff=v2[i]-v1[i];
    d+=fabs(diff);
  }

  return d;
}

template float metric2<float>(float *v1, float *v2, dim_ta m);
template double metric2<double>(double *v1, double *v2, dim_ta m);

template float manhattan2<float>(float *v1, float *v2, dim_ta m);
template double manhattan2<double>(double *v1, double *v2, dim_ta m);

template float quad_metric2<float>(float *v1, float *v2, dim_ta m);
template double quad_metric2<double>(double *v1, double *v2, dim_ta m);

template float max_metric2<float>(float *v1, float *v2, dim_ta m);
template double max_metric2<double>(double *v1, double *v2, dim_ta m);

//these two probably aren't all that useful (if any of them are...)
template float manhattan<float>(float *v1, float *v2, dim_ta m);
template double manhattan<double>(double *v1, double *v2, dim_ta m);

template float quad_metric4<float>(float *v1, float *v2, dim_ta m);
template double quad_metric4<double>(double *v1, double *v2, dim_ta m);

