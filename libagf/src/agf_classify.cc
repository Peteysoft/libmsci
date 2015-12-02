
//Copyright (C) 2007 Peter Mills.  All rights reserved.

#include <math.h>
#include <stdio.h>

#include "kextreme.h"
#include "quicksort.h"
#include "agf_lib.h"

#include <stdio.h>

using namespace std;
using namespace libpetey;

namespace libagf {

template <class real, class cls_t>
cls_t agf_classify(real **mat, 		//matrix of training data
		dim_ta D, 		//number of dimensions
		cls_t *cl,		//class labels
		nel_ta n, 		//number of samples
		cls_t ncls,		//number of class types
		real *vec, 		//test point
		real var[2],		//initial filter variance
		nel_ta k, 		//number of nearest neighbours to use
		real Wc, 		//objective total weight
		real *pdf,		//returned confidence cond. prob.
		agf_diag_param *diag_param,	//returned diagnostics
		flag_a joint)		//joint or conditional probabilities?
{

  real *d2;
  long *ind;			//indices for k nearest
  real *knearest;		//k nearest distances
  real *weight;		//the weights
  real tw;
  cls_t cls;			//final class
  real var_f;			//final filter variance
  //if we want joint probabilities:
  real norm;

  //first we calculate all the distances:
  d2=new real[n];
  for (nel_ta i=0; i<n; i++) {
    d2[i]=metric2(vec, mat[i], D);
    //printf("%g ", d2[i]);
  }
  //printf("\n");

  //select out the k nearest:
  knearest=new real[k];
  ind=new long[k];
  //sorter.decompose(knearest, ind, k);
  KLEAST_FUNC(d2, n, k, knearest, ind);

  //calculate the weights using the central "engine":
  weight=new real[k];
  diag_param->nd=AGF_CALC_W_FUNC(knearest, k, Wc, var, weight, var_f);
  //printf("f= %f\n", weight[k-1]/weight[0]);

  //calculate the non-normalized conditional probabilities:
  for (cls_t i=0; i<ncls; i++) pdf[i]=0;
  for (nel_ta i=0; i<k; i++) {
    if (cl[ind[i]]==-1) continue;		//exlude class values of -1
    pdf[cl[ind[i]]]+=weight[i];
  }

  //use these to select the class:
  cls=choose_class(pdf, ncls);

  //total of weights to normalize probabilities:
  tw=pdf[0];
  for (cls_t i=1; i<ncls; i++) tw+=pdf[i];

  //normalize the pdfs:
  if (joint) {
    //use the final filter width to normalize the pdf:
    norm=n*pow(sqrt(var_f*M_PI*2), D);
  } else {
    norm=tw;
  }
  for (cls_t i=0; i<ncls; i++) pdf[i]/=norm;

  //diagnostic parameter:
  diag_param->f=weight[k-1]/weight[0];
  diag_param->W=tw;

  //clean up:
  delete [] d2;
  delete [] knearest;
  delete [] ind;
  delete [] weight;

  return cls;

}

//this version uses all of the training data:
template <class real, class cls_t>
cls_t agf_classify(real **mat, 
		dim_ta D, 
		cls_t *cl,
		nel_ta n, 
		cls_t ncls,
		real *vec, 
		real var[2],
		real Wc, 
		real *pdf,
		agf_diag_param *diag_param, 
		flag_a joint)		
{

  real *d2;
  real *weight;		//the weights
  real tw;
  cls_t cls;			//final class
  real var_f;			//final filter variance

  real minw, maxw;

  //if we want joint probabilities:
  real norm;

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

  //calculate the non-normalized conditional probabilities:
  minw=1;
  maxw=0;
  for (cls_t i=0; i<ncls; i++) pdf[i]=0;
  for (nel_ta i=0; i<n; i++) {
    if (cl[i]==-1) continue;		//exclude class values of -1
    pdf[cl[i]]+=weight[i];
    //if (weight[i]<minw) minw=weight[i]; 
    //		else if (weight[i]>maxw) maxw=weight[i];
  }

  //use these to select the class:
  cls=choose_class(pdf, ncls);

  //total of weights to normalize probabilities:
  tw=pdf[0];
  for (cls_t i=1; i<ncls; i++) tw+=pdf[i];

  //normalize the pdfs:
  if (joint) {
    //use the final filter width to normalize the pdf:
    norm=n*pow(sqrt(var_f*M_PI*2), D);
  } else {
    norm=tw;
  }
  for (cls_t i=0; i<ncls; i++) pdf[i]/=norm;

  //diagnostics:
  //diag_param->f=minw/maxw;	//should normally be ~0, but we'll return it anyway...
  //diag_param->f=0;		//slows it down too much--for a useless diagnostic parm.
  diag_param->W=tw;

  //clean up:
  delete [] d2;
  delete [] weight;

  return cls;

}

template cls_ta agf_classify<float,cls_ta>(float **mat, dim_ta D, cls_ta *cl, 
		nel_ta n, cls_ta ncls, float *vec, float var[2], nel_ta k,
		float Wc, float *pdf, agf_diag_param *diag_param, flag_a joint);

template cls_ta agf_classify<double,cls_ta>(double **mat, dim_ta D, cls_ta *cl, 
		nel_ta n, cls_ta ncls, double *vec, double var[2], nel_ta k,
		double Wc, double *pdf, agf_diag_param *diag_param, flag_a joint);

template cls_ta agf_classify<float,cls_ta>(float **mat, dim_ta D, cls_ta *cl, 
		nel_ta n, cls_ta ncls, float *vec, float var[2], 
		float Wc, float *pdf, agf_diag_param *diag_param, flag_a joint);

template cls_ta agf_classify<double,cls_ta>(double **mat, dim_ta D, cls_ta *cl, 
		nel_ta n, cls_ta ncls, double *vec, double var[2], 
		double Wc, double *pdf, agf_diag_param *diag_param, flag_a joint);

}

