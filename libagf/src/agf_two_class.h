#ifndef AGF_TWO_CLASS_H_DEF
#define AGF_TWO_CLASS_H_DEF 1

#include "agf_defs.h"

namespace libagf {

  //returns the difference between the pdfs of two classes plus the gradient vector...
  template <class real>
  real dgrad(real **mat, 		//coordinate data sorted by class
		dim_ta D, 		//dimension of coordinate data
		nel_ta ntrain, 		//number of training samples
		nel_ta clind, 		//index of transition from class 0 to class 1
		real *vec, 		//test vector
		real var[2],		//filter variance brackets
		nel_ta k, 		//number of nearest neighbours to use 
					// in calculation
		real wc, 		//total of filter weights
		real *grad, 		//gradient of estimate
		agf_diag_param *diag_param); //diagnostic parameters

  //returns the difference between the pdfs of two classes but no gradient vector...
  template <class real>
  real dcalc(real **mat, 
		dim_ta D, 
		nel_ta ntrain, 
		nel_ta clind, 
		real *vec, 
		real var[2],
		nel_ta k, 
		real wc, 
		agf_diag_param *diag_param);

  //these versions use all the training data:

  //returns the difference between the pdfs of two classes plus the gradient vector...
  template <class real>
  real dgrad(real **mat, 
		dim_ta D, 
		nel_ta ntrain, 
		nel_ta clind, 
		real *vec, 
		real var[2],
		real wc, 
		real *grad, 
		agf_diag_param *diag_param);

  //returns the difference between the pdfs of two classes but no gradient vector...
  template <class real>
  real dcalc(real **mat, 
		dim_ta D, 
		nel_ta ntrain, 
		nel_ta clind, 
		real *vec, 
		real var[2],
		real wc, 
		agf_diag_param *diag_param);
}

#endif

