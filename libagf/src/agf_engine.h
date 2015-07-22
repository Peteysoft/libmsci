
//Copyright (C) 2007 Peter Mills.  All rights reserved.

#ifndef AGF_ENGINE_H_INCLUDED
#define AGF_ENGINE_H_INCLUDED

#include "agf_defs.h"

namespace libagf {

  //global variables (parameters):
  extern long agf_global_weights_maxiter; 
  extern real_a agf_global_weights_tol; 

  //given a set of distances,
  //calculates the weights of an "adaptive Gaussian filter"
  //--uses repeated squaring of the weights
  template <class real>
  iter_ta agf_calc_w(real *d2,		//distances squared
              nel_ta k,			//number of distances
              real Wc,			//objective total weight
	      real var[2],		//initial filter width
              real *weight,		//returned weights
              real &var_f);		//returned final filter variance

  //given a set of distances,
  //calculates the weights of an "adaptive Gaussian filter"
  //--uses supernewton root-finding algorithm
  template <class real>
  iter_ta agf_calc_w2(real *d2,		//distances squared
              nel_ta k,			//number of distances
              real Wc,			//objective total weight
	      real var[2],		//filter width brackets
              real *weight,		//returned weights
              real &var_f);		//returned final filter variance

  //calculate the gradient of the weights:
  template <class real>
  void agf_grad_w(real **x,		//matrix of samples
		dim_ta n,		//# of dimensions
		real *xvec,		//test point
		real *w,		//weights
		real *d2,		//squared distances
		nel_ta k,		//number of weights and indices
		real var_f,		//filter width
		real **dwdx);		//returned gradient vectors
 
}

#endif
