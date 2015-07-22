#ifndef SAMPLE_CLASS_BORDERS__H
#define SAMPLE_CLASS_BORDERS__H

#include <stdio.h>

#include <gsl/gsl_rng.h>

#include "tree_lg.h"
#include "agf_defs.h"

namespace libagf {

  //- parameter structure for general borders training when this borders training
  //    is done with discretely sampled training data
  //- to be passed to "rfunc", and "sample" 
  //- use rparam for different versions of "rfunc"
  //    which calculates the difference in conditional probabilities
  template <class real>
  struct bordparam {
    void *rparam;		//pass to rfunc

    //training data:
    real **train;		//sample locations, sorted by class
    nel_ta n;			//number of samples
    dim_ta D;			//number of dimensions
    nel_ta ind;			//class 0/class 1

    //keep track of samples to avoid duplicates (large number of samples):
    libpetey::tree_lg<int64_t> *sind1;

    //for sampling small datasets:
    long s;			//current sample number
    int64_t *sind2;		//random indices

    //maximum possible number of samples:
    int64_t maxsample;

    //generate random numbers:
    gsl_rng *rann;
  };

  //stick the training data into the parameter structure:
  template <class real>
  int bordparam_init(bordparam<real> *param,
		real **train, 			//training samples
		dim_ta D, 			//number of dimensions
		nel_ta n, 			//number of samples
		nel_ta ind,			//index of class label change
		int smallflag=0);		//small dataset?

  template <class real>
  void bordparam_clean(bordparam<real> *param);

  //returns a sample point from each class:
  //negative return values indicates error (no more samples left...)
  template <class real>
  int oppositesample(void *param, 		//parameters: training samples etc.
		real *x1, 			//sample from first class
		real *x2);			//sample from second class

  //for small datasets: 
  //(returns every combination of samples)
  template <class real>
  int oppositesample_small(void *param, 
		real *x1, 
		real *x2);

  //given a function returning the difference in conditional probabilities
  //and a function to sample opposite classes, samples the borders between
  //a pair of classes; returns number of border samples found
  template <class real>
  nel_ta sample_class_borders(real (*rfunc) (real *, void *, real *), 	
			//returns difference in conditional prob. plus derivs
		int (*sample) (void *, real *, real *),	//returns a random sample from each class
		void *param,			//these are just along for the ride
		nel_ta n,			//number of times to sample
		dim_ta D,			//number of dimensions
		real tol,			//desired tolerance
		iter_ta maxit,			//maximum number of iterations
		real **border,			//returned border samples
		real **gradient,		//returned border gradients
		real rthresh=0);		//location of Bayesian border

  //now for classification with borders & gradients:
  //classify a test point based on a set of border samples:
  //this version return the index of the nearest sample and the disance
  //also, it doesn't take the hyperbolic tangent of the result
  template <class real>
  real border_classify0(real **brd,	//border samples
                real **grd,		//gradient vectors
                dim_ta D,		//number of dimensions
                nel_ta n,		//number of border samples
                real *x, 		//test point
		nel_ta &k,		//index of nearest sample
		real &d);		//distance

  //works as usual: given a set of border sample and accompanying
  //gradients, returns the approximated difference in conditional probabilities:
  template <class real>
  inline real border_classify(real **brd,	//border samples
                real **grd,		//gradient vectors
                dim_ta D,		//number of dimensions
                nel_ta n,		//number of border samples
                real *x) { 		//test point
    real d;		//both are throw-aways
    nel_ta k;
    return tanh(border_classify0(brd, grd, D, n, x, k, d));
  }

  //test codes:
  template <class real>
  int test_oppositesample(nel_ta n1, nel_ta n2);
  template <class real>
  int test_oppositesample_small(nel_ta n);

} //end namespace libagf

#endif

