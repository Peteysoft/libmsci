//
// Routines for AGF interpolation.
//

#ifndef AGFGAF_H_INCLUDED
#define AGFGAF_H_INCLUDED 1

#include "agf_defs.h"

namespace libagf {

  //interpolation routines:
  template <class real>
  real adgaf(real **xmat,		//location of samples (n row by D col)
            dim_ta D,			//number of dimensions
            real *y,			//samples of function
            nel_ta n,			//number of samples
            real *xvec,			//interpolation point
            real var[2],		//filter variance bracket
            nel_ta k,			//number of nearest neighbours
            real Wc,			//objective total weight
            agf_diag_param *diag_param);//diagnostic parameters

  //uses all of the data (no k parameter):
  template <class real>
  real adgaf(real **mat, 
		dim_ta D, 
		real *y, 
		nel_ta n, 
		real *vec, 
		real var[2], 
		real wc, 
		agf_diag_param* diag_param);

  //returns RMS error estimates:
  template <class real>
  real adgaf_err(real **mat, 
		dim_ta D, 
		real *y, 
		nel_ta n, 
		real *vec, 
		real var[2], 
		nel_ta k, 
		real wc, 
		real &err, 		//returned error estimates
		agf_diag_param *diag_param);

  //uses all of the data and returns RMS error estimates:
  template <class real>
  real adgaf_err(real **mat, 
		dim_ta D, 
		real *y, 
		nel_ta n, 
		real *vec, 
		real var[2], 
		real wc, 
		real &err, 
		agf_diag_param *diag_param);
}

#endif

