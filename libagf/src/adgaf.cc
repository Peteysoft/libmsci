//
// This software is released under the following terms:
//
// 1. No commercial use.
// 2. Copies and derivative works are free to use and modify.
// 3. Attribution must be given to all contributors of both original and derivative works.
//
// (like I care if anyone steals this heap of steaming turds... )
//
// Authors:
//
// 2017-07-16 Peter Mills: added license information 
//

//
// Codes for AGF interpolation ("adgaf").
//

#include <math.h>
#include <stdio.h>

#include "kextreme.h"
#include "quicksort.h"

#include "agf_lib.h"

#include <stdio.h>

using namespace std;
using namespace libpetey;
using namespace libagf;

//interpolates an irregularly gridded set of samples using an
//adaptive Guassian filter:
template <class real>
real libagf::adgaf(real **xmat,       //location of samples (n row by D col)
            dim_ta D,             //number of dimensions
            real *y,           //samples of function
            nel_ta n,             //number of samples
            real *xvec,        //interpolation point
	    real var[2],	//initial filter variance
            nel_ta k,             //number of nearest neighbours
            real Wc,           //objective total weight
	    agf_diag_param *diag_param)		//diagnostic parameters
{

  real *knearest;	//distances of k nearest
  long *ind;		//indices of k nearest
  real *weight;		//the weights
  real tw;		//sum of the weights
  real result;		//final value of the interpolate
  real var_f;		//final filter variance

  //select out the k nearest:
  knearest=new real[k];
  ind=new long[k];
  kinearest(xmat, n, D, xvec, k, knearest, ind);
  
  //calculate the weights:
  weight=new real[k];
  diag_param->nd=AGF_CALC_W_FUNC(knearest, k, Wc, var, weight, var_f);

  //apply weights and normalize by total weight:
  tw=0;
  result=0;
  for (nel_ta i=0; i<k; i++) {
    result+=weight[i]*y[ind[i]];
    tw+=weight[i];
  }

  //set diagnostic parameters:
  diag_param->W=tw;
  diag_param->f=weight[k-1]/weight[0];

  //clean up:
  delete [] knearest;
  delete [] ind;
  delete [] weight;

  return result/tw;
}

template float adgaf<float>(float **xmat, dim_ta D, float *y, nel_ta n,  
            float *xvec, float var[2], nel_ta k, 
            float Wc, agf_diag_param *diag_param);	
template double adgaf<double>(double **xmat, dim_ta D, double *y, nel_ta n,  
            double *xvec, double var[2], nel_ta k, 
            double Wc, agf_diag_param *diag_param);	

//interpolates an irregularly gridded set of samples using an adaptive Guassian filter:
//
//xmat:		location of samples
//m:		number of dimensions
//y:		samples of function
//n:		number of samples
//vec:		location for which interpolates are desired
//var_0:	initial filter width
//k:		number of nearest neighbours to use in the interpolation
//wc:		desired total value of the weights
//err:		error estimate
template <class real>
real libagf::adgaf_err(real **xmat, dim_ta D, real *y, nel_ta n, real *vec, real var[2], 
		nel_ta k, real wc, real &err, agf_diag_param *diag_param) {
  real *knearest;		//distances of k nearest
  long *ind;			//indices of k nearest
  real *weight;		//the weights
  real *grad;			//the gradient vector
  real varf;			//final filter width
  real result;			//final value of the interpolate
  real diff;			//intermediate result;
  real tw;
  real **mat2;			//k nearest training samples

  real **dwdx;			//gradients of weights

  //select out the k nearest:
  knearest=new real[k];
  ind=new long[k];
  kinearest(xmat, n, D, vec, k, knearest, ind);

  //rearrange the training samples:
  mat2=new real * [k];
  for (nel_ta i=0; i<k; i++) mat2[i]=xmat[ind[i]];

  //calculate the weights:
  weight=new real[k];
  diag_param->nd=AGF_CALC_W_FUNC(knearest, k, wc, var, weight, varf);

  //apply them to the result:
  tw=0;
  result=0;
  for (nel_ta i=0; i<k; i++) {
    tw+=weight[i];
    result+=weight[i]*y[ind[i]];
  }
  
  //calculate the gradient vector:
  //printf("Gradient vector: ");
  //allocate space:
  dwdx=new real * [k];
  dwdx[0]=new real[k*D];
  for (nel_ta i=1; i<k; i++) dwdx[i]=dwdx[0]+i*D;
  grad=new real[D];

  agf_grad_w(xmat, D, vec, weight, knearest, k, varf, dwdx);
  for (dim_ta i=0; i<D; i++) {
    grad[i]=0;
    for (nel_ta j=0; j<k; j++) {
      grad[i]+=dwdx[j][i]*(y[ind[j]]-result);
    }
    grad[i]=grad[i]/tw/varf;
    //printf("%g ", grad[i]);
  }
  //printf("\n");
  
  //use the gradient vector to approximate the error:
  err=0;
  for (nel_ta i=0; i<k; i++) {
    diff=0;
    for (dim_ta j=0; j<D; j++) {
      diff+=grad[j]*(mat2[i][j]-vec[j]);
    }
    diff=diff+result-y[ind[i]];
    err+=weight[i]*diff*diff;
  }
  
  err=err/tw;
  result=result/tw;
      
  //set diagnostic parameters:
  diag_param->W=tw;
  diag_param->f=weight[k-1]/weight[0];

//  printf("Interpolated total weight: %g\n", tw);

  //clean up:
  delete [] dwdx[0];
  delete [] dwdx;

  delete [] knearest;
  delete [] ind;
  delete [] weight;
  delete [] grad;
  delete [] mat2;

  return result;

}

template float adgaf_err<float>(float **xmat, dim_ta D, float *y, nel_ta n, 
		float *vec, float var[2], nel_ta k, float wc, 
		float &err, agf_diag_param *diag_param);
template double adgaf_err<double>(double **xmat, dim_ta D, double *y, nel_ta n, 
		double *vec, double var[2], nel_ta k, double wc, 
		double &err, agf_diag_param *diag_param);

//these versions take all the training data:

//interpolates an irregularly gridded set of samples using an
//adaptive Guassian filter:

template <class real>
real libagf::adgaf(real **xmat,       //location of samples (n row by D col)
            dim_ta D,             //number of dimensions
            real *y,           //samples of function
            nel_ta n,             //number of samples
            real *xvec,        //interpolation point
	    real var[2],	//initial filter variance
            real Wc,           //objective total weight
	    agf_diag_param *diag_param)
{

  real *d2;                  //all the distances
  real *weight;              //the weights
  real tw;            		//sum of the weights
  real result;                 //final value of the interpolate
  real var_f;			//final filter variance

  //first we calculate all the distances:
  d2=new real[n];
  for (nel_ta i=0; i<n; i++) d2[i]=metric2(xvec, xmat[i], D);

  //calculate the weights:
  weight=new real[n];
  diag_param->nd=AGF_CALC_W_FUNC(d2, n, Wc, var, weight, var_f);

  //apply weights and normalize by total weight:
  tw=0;
  result=0;
  for (nel_ta i=0; i<n; i++) {
    result+=weight[i]*y[i];
    tw+=weight[i];
  }

  //set diagnostic parameters:
  diag_param->W=tw;
  diag_param->f=0;

  //clean up:
  delete [] d2;
  delete [] weight;

  return result/tw;
}

template float adgaf<float>(float **xmat, dim_ta D, float *y, nel_ta n, float *xvec, 
		float var[2], float Wc, agf_diag_param *diag_param);
template double adgaf<double>(double **xmat, dim_ta D, double *y, nel_ta n, double *xvec, 
		double var[2], double Wc, agf_diag_param *diag_param);

//interpolates an irregularly gridded set of samples using an adaptive Guassian filter:
//
//xmat:		location of samples
//m:		number of dimensions
//y:		samples of function
//n:		number of samples
//vec:		location for which interpolates are desired
//var_0:	initial filter width
//k:		number of nearest neighbours to use in the interpolation
//wc:		desired total value of the weights
//err:		error estimate
template <class real>
real libagf::adgaf_err(real **xmat, dim_ta D, real *y, nel_ta n, real *vec, real var[2], 
		real wc, real &err, agf_diag_param *diag_param) {
  real *d2;			//all the distances
  real *weight;		//the weights
  real *grad;			//the gradient vector
  real varf;			//final filter width
  real result;			//final value of the interpolate
  real diff;			//intermediate result;
  real tw;

  real **dwdx;			//gradients of weights

  //first we calculate all the distances:
  d2=new real[n];
  for (nel_ta i=0; i<n; i++) {
    d2[i]=metric2(vec, xmat[i], D);
  }

  //calculate the weights:
  weight=new real[n];
  diag_param->nd=AGF_CALC_W_FUNC(d2, n, wc, var, weight, varf);

  //apply them to the result:
  tw=0;
  result=0;
  for (nel_ta i=0; i<n; i++) {
    tw+=weight[i];
    result+=weight[i]*y[i];
  }
  
  //calculate the gradient vector:
  //printf("Gradient vector: ");
  //allocate space:
  dwdx=new real * [n];
  dwdx[0]=new real[n*D];
  for (nel_ta i=1; i<n; i++) dwdx[i]=dwdx[0]+i*D;
  grad=new real[D];

  agf_grad_w(xmat, D, vec, weight, d2, n, varf, dwdx);

  //use the gradient vector to approximate the error:
  for (dim_ta i=0; i<D; i++) {
    grad[i]=0;
    for (nel_ta j=0; j<n; j++) {
      grad[i]+=dwdx[j][i]*(y[j]-result);
    }
    grad[i]=grad[i]/tw/varf;
    //printf("%g ", grad[i]);
  }
  //printf("\n");
  
  err=0;
  for (nel_ta i=0; i<n; i++) {
    diff=0;
    for (dim_ta j=0; j<D; j++) {
      diff+=grad[j]*(xmat[i][j]-vec[j]);
    }
    diff=diff+result-y[i];
    err+=weight[i]*diff*diff;
  }
  
  err=err/tw;
  result=result/tw;
      
//  printf("Interpolated total weight: %g\n", tw);

  //clean up:
  delete [] dwdx[0];
  delete [] dwdx;

  //set diagnostic parameters:
  diag_param->W=tw;
  diag_param->f=0;

  delete [] d2;
  delete [] weight;
  delete [] grad;

  return result;

}

template float adgaf_err<float>(float **xmat, dim_ta D, float *y, nel_ta n, float *vec, 
		float var[2], float wc, float &err, agf_diag_param *diag_param);
template double adgaf_err<double>(double **xmat, dim_ta D, double *y, nel_ta n, double *vec, 
		double var[2], double wc, double &err, agf_diag_param *diag_param);
