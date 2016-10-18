
//Copyright (C) 2007 Peter Mills.  All rights reserved.

#include <stdio.h>
#include <math.h>
#include <sys/timeb.h>

#include "agf_lib.h"
#include "kextreme.h"
#include "quicksort.h"

using namespace std;
using namespace libpetey;

namespace libagf {

//returns the difference between the pdfs of two classes plus the gradient vector...
template <class real>
real dgrad(real **mat, dim_ta D, nel_ta ntrain, nel_ta clind, real *vec, real var[2], 
		nel_ta k, real minw, real *grad, agf_diag_param *diag_param) {
  real *knearest;		//distances of k nearest neighbours
  long *ind;			//indexes of k nearest neighbours
  real *weight;		//the current value for the weights
  real tw;			//total of weights

  real dplus, dminus, d;

  real **dwdx;			//gradients of weights
  real varf;			//final filter width

  real **mat2;			//k nearest training samples

  //select out the k nearest:
  knearest=new real[k];
  ind=new long[k];
  kinearest(mat, ntrain, D, vec, k, knearest, ind);

  //select out training samples:
  mat2=new real *[k];
  for (nel_ta i=0; i<k; i++) mat2[i]=mat[ind[i]];

  /*
  for (nel_ta i=0; i<k; i++) printf("%g ", knearest[i]);
  printf("\n");
  */

  //calculate the weights:
  weight=new real[k];
  diag_param->nd=AGF_CALC_W_FUNC(knearest, k, minw, var, weight, varf);

  //we reduce the conditional probabilities to only one number:
  dplus=0;
  dminus=0;
  for (nel_ta i=0; i<k; i++) {
    if (ind[i] < clind) {
      dminus+=weight[i];
    } else {
      dplus+=weight[i];
    }
  }
  tw=dplus+dminus;
  d=(dplus-dminus)/tw;

  //now to calculate the gradient vector:
  //allocate space:
  dwdx=new real *[k];
  dwdx[0]=new real[k*D];
  for (nel_ta i=1; i<k; i++) dwdx[i]=dwdx[0]+i*D;

  agf_grad_w(mat2, D, vec, weight, knearest, k, varf, dwdx);

  for (dim_ta i=0; i<D; i++) {
    //calculate the total change in the weights:
    grad[i]=0;
    for (nel_ta j=0; j<k; j++) {
      if (ind[j] < clind) {
        grad[i]-=dwdx[j][i];
      } else {
        grad[i]+=dwdx[j][i];
      }
    }
    grad[i]=grad[i]/tw;
  }

  //diagnostic parameter:
  diag_param->f=weight[k-1]/weight[0];
  diag_param->W=tw;

  //clean up:
  delete [] dwdx[0];
  delete [] dwdx;

  delete [] knearest;
  delete [] ind;
  delete [] weight;
  delete [] mat2;

  return d;

}

//returns the difference between the pdfs of two classes but no gradient vector...
template <class real>
real dcalc(real **mat, dim_ta D, nel_ta ntrain, nel_ta clind, real *vec, real var[2],
		nel_ta k, real minw, agf_diag_param *diag_param) {
  real *knearest;		//distances of k nearest neighbours
  long *ind;			//indexes of k nearest neighbours
  real *weight;		//the current value for the weights
  real tw;			//total of weights

  real d, dplus, dminus;

  real varf;

  //select out the k nearest:
  knearest=new real[k];
  ind=new long[k];
  kinearest(mat, ntrain, D, vec, k, knearest, ind);

  /*
  for (nel_ta i=0; i<k; i++) printf("%g ", knearest[i]);
  printf("\n");
  */

  //calculate the weights:
  weight=new real[k];
  diag_param->nd=AGF_CALC_W_FUNC(knearest, k, minw, var, weight, varf);

  //we reduce the conditional probabilities to only one number:
  dplus=0;
  dminus=0;
  for (nel_ta i=0; i<k; i++) {
    if (ind[i] < clind) {
      dminus+=weight[i];
    } else {
      dplus+=weight[i];
    }
  }
  tw=dplus+dminus;
  d=(dplus-dminus)/tw;

  //diagnostic parameters:
  diag_param->f=weight[k-1]/weight[0];
  diag_param->W=tw;

  //clean up:
  delete [] knearest;
  delete [] ind;
  delete [] weight;

  return d;

}

//these versions use all the training data:

//returns the difference between the pdfs of two classes plus the gradient vector...
template <class real>
real dgrad(real **mat, dim_ta D, nel_ta ntrain, nel_ta clind, real *vec, real var[2],
		real minw, real *grad, agf_diag_param *diag_param) {
  real *d2;			//the distances (squared)
  real *weight;		//the current value for the weights
  //real dw[k];			//change in the weights
  real tw;			//total of weights

  real d, dplus, dminus;

  real **dwdx;			//gradients of weights
  real varf;			//final filter width

  real **mat2;			//k nearest training samples

  //first we calculate all the distances:
  d2=new real[ntrain];
  for (nel_ta i=0; i<ntrain; i++) {
    d2[i]=metric2(vec, mat[i], D);
  }

  /*
  for (nel_ta i=0; i<k; i++) printf("%g ", knearest[i]);
  printf("\n");
  */

  //calculate the weights:
  weight=new real[ntrain];
  diag_param->nd=AGF_CALC_W_FUNC(d2, ntrain, minw, var, weight, varf);

  //we reduce the conditional probabilities to only one number:
  dplus=0;
  dminus=0;
  for (nel_ta i=0; i<clind; i++) {
    dminus+=weight[i];
  }
  for (nel_ta i=clind; i<ntrain; i++) {
    dplus+=weight[i];
  }
  tw=dplus+dminus;
  d=(dplus-dminus)/tw;

  //now to calculate the gradient vector:
  //allocate space:
  dwdx=new real *[ntrain];
  dwdx[0]=new real[ntrain*D];
  for (nel_ta i=1; i<ntrain; i++) dwdx[i]=dwdx[0]+i*D;

  agf_grad_w(mat, D, vec, weight, d2, ntrain, varf, dwdx);

  for (dim_ta i=0; i<D; i++) {
    //calculate the total change in the weights:
    grad[i]=0;
    for (nel_ta j=0; j<clind; j++) grad[i]-=dwdx[j][i];
    for (nel_ta j=clind; j<ntrain; j++) grad[i]+=dwdx[j][i];
    grad[i]=grad[i]/tw;
  }

  //diagnostic parameters:
  //diag_param->f=0;
  diag_param->W=tw;

  //clean up:
  delete [] dwdx[0];
  delete [] dwdx;

  delete [] d2;
  delete [] weight;

  return d;

}

//returns the difference between the pdfs of two classes but no gradient vector...
template <class real>
real dcalc(real **mat, dim_ta D, nel_ta ntrain, nel_ta clind, real *vec, real var[2],
		real minw, agf_diag_param *diag_param)
		{
  real *d2;			//the distances (squared)
  real *weight;		//the current value for the weights
  real tw;			//total of weights

  real d, dplus, dminus;

  real **dwdx;			//gradients of weights
  real varf;

  //first we calculate all the distances:
  d2=new real[ntrain];
  for (nel_ta i=0; i<ntrain; i++) {
    d2[i]=metric2(vec, mat[i], D);
  }

  /*
  for (nel_ta i=0; i<k; i++) printf("%g ", knearest[i]);
  printf("\n");
  */

  //calculate the weights:
  weight=new real[ntrain];
  diag_param->nd=AGF_CALC_W_FUNC(d2, ntrain, minw, var, weight, varf);

  //we reduce the conditional probabilities to only one number:
  dplus=0;
  dminus=0;
  for (nel_ta i=0; i<clind; i++) {
    dminus+=weight[i];
  }
  for (nel_ta i=clind; i<ntrain; i++) {
    dplus+=weight[i];
  }
  tw=dplus+dminus;
  d=(dplus-dminus)/tw;

  //diagnostic parameters:
  //diag_param->f=0;
  diag_param->W=tw;

  //clean up:
  delete [] d2;
  delete [] weight;

  return d;

}

template double dgrad<double>(double **mat, dim_ta D, nel_ta ntrain, nel_ta clind, double *vec, double var[2], 
		nel_ta k, double minw, double *grad, agf_diag_param *diag_param);

template double dcalc(double **mat, dim_ta D, nel_ta ntrain, nel_ta clind, double *vec, double var[2],
		nel_ta k, double minw, agf_diag_param *diag_param);

template double dgrad<double>(double **mat, dim_ta D, nel_ta ntrain, nel_ta clind, double *vec, double var[2],
		double minw, double *grad, agf_diag_param *diag_param);

template double dcalc<double>(double **mat, dim_ta D, nel_ta ntrain, nel_ta clind, double *vec, double var[2],
		double minw, agf_diag_param *diag_param);

template float dgrad<float>(float **mat, dim_ta D, nel_ta ntrain, nel_ta clind, float *vec, float var[2], 
		nel_ta k, float minw, float *grad, agf_diag_param *diag_param);

template float dcalc(float **mat, dim_ta D, nel_ta ntrain, nel_ta clind, float *vec, float var[2],
		nel_ta k, float minw, agf_diag_param *diag_param);

template float dgrad<float>(float **mat, dim_ta D, nel_ta ntrain, nel_ta clind, float *vec, float var[2],
		float minw, float *grad, agf_diag_param *diag_param);

template float dcalc<float>(float **mat, dim_ta D, nel_ta ntrain, nel_ta clind, float *vec, float var[2],
		float minw, agf_diag_param *diag_param);

}

