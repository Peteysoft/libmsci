//
// Codes for AGF density estimation.
//

#include <math.h>

#include "kextreme.h"
#include "quicksort.h"
#include "agf_lib.h"

using namespace std;
using namespace libpetey;

//do this the easy way:
namespace libagf {

//estimates the pdf at a test point based on a set of samples:
//
//mat is a set of n vectors with m elements
//vec is the vector we want to classify
//k is the number of nearest neighbours to use in the classification
//minw is the minimum un-normalized total of the weights
//
//returns the pdf
template <class real>
real agf_calc_pdf(real **mat, dim_ta D, nel_ta n, real *vec, real var[2], 
			nel_ta k, real Wc, agf_diag_param *diag_param) {
  real var_f;			//final value of the filter width (as variance)
  real tw;			//total weight
  real *kn;			//distances of k nearest neighbours
  real *weight;			//the current value for the weights
  real norm;			//normalisation coeff.
  real pdf;			//final calculated value of pdf

  //select out the k nearest:
  kn=new real[k];
  knearest(mat, n, D, vec, k, kn);

  //calculate the weights using the central "engine":
  weight=new real[k];
  diag_param->nd=AGF_CALC_W_FUNC(kn, k, Wc, var, weight, var_f);
  tw=0;
  for (nel_ta i=0; i<k; i++) tw+=weight[i];

  //use the final filter width to normalize the pdf:
  norm=pow(sqrt(var_f*M_PI*2), D);
  
  //norm=pow(var_f*M_PI*2, m/2.);
  //printf("var_f=%g, tw=%g, norm=%g\n", var_f, tw, norm);

  pdf=tw/norm/n;

  //set the diagnostic parameters:
  diag_param->f=weight[k-1]/weight[0];
  diag_param->W=tw;

  delete [] kn;
  delete [] weight;

  return pdf;

}

//estimates the pdf at a test point based on a set of samples:
//this version uses all of the training data...
//
//mat is a set of n vectors with m elements
//vec is the vector we want to classify
//minw is the minimum un-normalized total of the weights
//
//returns the pdf
template <class real>
real agf_calc_pdf(real **mat, dim_ta D, nel_ta n, real *vec, real var[2], 
			real Wc, agf_diag_param *diag_param) {

  real *d2;			//the distances (squared)
  real var_f;			//final value of the filter width (as variance)
  real tw;			//total weight
  real *weight;			//the current value for the weights
  real norm;			//normalisation coeff.
  real pdf;			//final calculated value of pdf

  //first we calculate all the distances:
  d2=new real[n];
  for (nel_ta i=0; i<n; i++) {
    d2[i]=metric2(vec, mat[i], D);
    //printf("%g ", d2[i]);
  }
  //printf("\n");

  //calculate the weights using the central "engine":
  weight=new real[n];
  diag_param->nd=AGF_CALC_W_FUNC(d2, n, Wc, var, weight, var_f);
  tw=0;
  for (nel_ta i=0; i<n; i++) tw+=weight[i];

  //use the final filter width to normalize the pdf:
  norm=pow(sqrt(var_f*M_PI*2), D);
  
  //norm=pow(var_f*M_PI*2, m/2.);
  //printf("var_f=%g, tw=%g, norm=%g\n", var_f, tw, norm);

  pdf=tw/norm/n;

  //set the diagnostic parameters:
  //diag_param->f=0;
  diag_param->W=tw;

  delete [] d2;
  delete [] weight;

  return pdf;

}

template float agf_calc_pdf<float>(float **mat, dim_ta D, nel_ta n, float *vec, 
		float var[2], nel_ta k, float Wc, agf_diag_param *diag_param);
template double agf_calc_pdf<double>(double **mat, dim_ta D, nel_ta n, double *vec, 
		double var[2], nel_ta k, double Wc, agf_diag_param *diag_param);

template float agf_calc_pdf<float>(float **mat, dim_ta D, nel_ta n, float *vec, 
		float var[2], float Wc, agf_diag_param *diag_param);
template double agf_calc_pdf<double>(double **mat, dim_ta D, nel_ta n, 
		double *vec, double var[2], double Wc, agf_diag_param *diag_param);

}
