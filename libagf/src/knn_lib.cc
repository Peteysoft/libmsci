
//Copyright (C) 2007 Peter Mills.  All rights reserved.

#include <stdio.h>
#include <math.h>
#include <assert.h>

#include "gsl/gsl_sf.h"

#include "agf_lib.h"
#include "kextreme.h"
#include "quicksort.h"
#include "peteys_tmpl_lib.h"

using namespace std;
using namespace libpetey;
using namespace libagf;

template <class real, class cls_t>
cls_t libagf::knn(real (* metric) (real *, real *, dim_ta),
		real **mat, dim_ta m, nel_ta n, cls_t *cl, cls_t ncl, real *vec,
		nel_ta k, real *pdf, flag_a joint) {
  real *d2;			//the distances (squared)
  real *knearest;		//distances of k nearest neighbours
  long *ind;			//indices of k nearest neighbours
  cls_t cls;			//estimated class

  //for joint probabilities:
  real r, V;			//radius, volume of hypersphere
  real sqrtpir;		//intermediate value

  //first we calculate all the distances:
  d2=new real[n];
  for (nel_ta i=0; i<n; i++) {
    //d2[i]=(*metric) (vec, mat[i], m);
    d2[i]=metric2 (vec, mat[i], m);
  }

  knearest=new real[k+joint];
  ind=new long[k+joint];
  KLEAST_FUNC(d2, n, k+joint, knearest, ind);
 
  for (cls_t i=0; i<ncl; i++) pdf[i]=0;
  for (nel_ta i=0; i<k; i++) {
    if (cl[ind[i]]==-1) continue;		//exclude class values of -1
    pdf[cl[ind[i]]]++;
  }

  cls=choose_class(pdf, ncl);

  //calculate normalization coefficient:
  if (joint) {
    r=(sqrt(knearest[k-1])+sqrt(knearest[k]))/2;
    sqrtpir=sqrt(M_PI)*r;
    V=1;
    for (dim_ta i=0; i<m; i++) V*=sqrtpir;
    V=V/gsl_sf_gamma(m/2.+1)*(real) n;
  } else {
    V=(real) k;
  }

  for (cls_t i=0; i<ncl; i++) pdf[i]/=V;

  //clean up:
  delete [] d2;
  delete [] knearest;
  delete [] ind;

  return cls;

}

template <class real>
real libagf::knn_pdf(real **mat, dim_ta m, nel_ta n, real *vec, nel_ta k) {
  real *d2;			//the distances (squared)
  real *knearest;		//distances of k nearest neighbours
  //real *kn2;		//distances of k nearest neighbours
  real r, V;			//radius, volume of hypersphere
  real sqrtpir;		//intermediate value
  real pdf;

  //first we calculate all the distances:
  d2=new real[n];
  for (nel_ta i=0; i<n; i++) {
    d2[i]=metric2(vec, mat[i], m);
  }

  //select out the k nearest:
  knearest=new real[k+1];
  KLEAST_FUNC(d2, n, k+1, knearest);
  //kn2=new real[k+1];
  //kleast(d2, n, k+1, kn2);
/*
  heapsort_inplace(knearest, k+1);
  for (nel_ta i=0; i<=k; i++) {
    printf("%g %g\n", knearest[i], kn2[i]);
    assert(knearest[i]==kn2[i]);
  }
  delete [] kn2;
*/

  //pdf is k divided by volume of hypersphere of radius intermediate
  //between farthest of the k samples and the next farthest:
  //calculate volume of hypersphere:
  r=(sqrt(knearest[k-1])+sqrt(knearest[k]))/2;
  sqrtpir=sqrt(M_PI)*r;
  V=1;
  for (dim_ta i=0; i<m; i++) V*=sqrtpir;
  V=V/gsl_sf_gamma(m/2.+1);

  //calculate pdf:
  pdf=(real) k/V/(real) n;
  //printf("r=%g, V=%g\n", r, V);

  //clean up:
  delete [] d2;
  delete [] knearest;
  
  return pdf;

}

template <class real>
real libagf::int_knn(real (* metric) (real *, real *, dim_ta),
	real **mat, dim_ta D, nel_ta n, real *ord, real *vec, nel_ta k) {
  real *d2;
  real *knearest;
  long *ind;
  real result;

  d2=new real[n];
  knearest=new real[k];
  ind=new long[k];

  for (nel_ta i=0; i<n; i++) d2[i]=(*metric) (mat[i], vec, D);
  KLEAST_FUNC(d2, n, k, knearest, ind);
  result=0;
  for (nel_ta i=0; i<k; i++) result+=ord[ind[i]];

  delete [] d2;
  delete [] knearest;
  delete [] ind;

  return result/k;
}

template <class real>
real libagf::int_knn(real (* metric) (real *, real *, dim_ta),
	real **mat, dim_ta D, nel_ta n, real *ord, real *vec, nel_ta k, real &rms) {
  real *d2;
  real *knearest;
  long *ind;
  real diff;
  real result;

  d2=new real[n];
  knearest=new real[k];
  ind=new long[k];

  for (nel_ta i=0; i<n; i++) d2[i]=(*metric) (mat[i], vec, D);

  KLEAST_FUNC(d2, n, k, knearest, ind);
  result=0;
  for (nel_ta i=0; i<k; i++) result+=ord[ind[i]];
  result/=k;

  rms=0;
  for (nel_ta i=0; i<k; i++) {
    diff=ord[ind[i]]-result;
    rms+=diff*diff;
  }
  rms=sqrt(rms/(k-1));

  delete [] d2;
  delete [] knearest;
  delete [] ind;

  return result/k;
}

template cls_ta knn<float, cls_ta>(float (*metric) (float *, float *, dim_ta),
		float **mat, dim_ta m, nel_ta n, cls_ta *cl, 
		cls_ta ncl, float *vec, nel_ta k, float *pdf, flag_a joint);

template float knn_pdf<float>(float **mat, dim_ta m, nel_ta n, float *vec, nel_ta k);

template float int_knn<float>(float (*metric) (float *, float *, dim_ta),
	float **mat, dim_ta D, nel_ta n, float *ord, float *vec, nel_ta k);

template float int_knn<float>(float (*metric) (float *, float *, dim_ta),
	float **mat, dim_ta D, nel_ta n, float *ord, float *vec, nel_ta k, float &rms);

template cls_ta knn<double, cls_ta>(double (*metric) (double *, double *, dim_ta),
		double **mat, dim_ta m, nel_ta n, cls_ta *cl, 
		cls_ta ncl, double *vec, nel_ta k, double *pdf, flag_a joint);

template double knn_pdf<double>(double **mat, dim_ta m, nel_ta n, double *vec, nel_ta k);

template double int_knn<double>(double (*metric) (double *, double *, dim_ta),
	double **mat, dim_ta D, nel_ta n, double *ord, double *vec, nel_ta k);

template double int_knn<double>(double (*metric) (double *, double *, dim_ta),
	double **mat, dim_ta D, nel_ta n, double *ord, double *vec, nel_ta k, double &rms);

