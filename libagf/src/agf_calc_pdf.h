//
// AGF density estimation.
//

#ifndef AGF_CALC_PDF_H_INCLUDED
#define AGF_CALC_PDF_H_INCLUDED 1

#include "agf_defs.h"

namespace libagf {

  //for PDF calculation:
  //restricts calculations to a set of k-nearest neighbours:
  template <class real>
  real agf_calc_pdf(real **mat, 		//coordinate samples
		dim_ta D, 			//dimension of samples
		nel_ta n, 			//number of samples
		real *vec, 			//test point
		real var[2], 			//filter variance brackets
		nel_ta k, 			//number of nearest neighbours
		real wc, 			//total of filter weights
		agf_diag_param *diag_param);	//diagnostic params

  //uses all of the training data:
  template <class real>
  real agf_calc_pdf(real **mat, 
		dim_ta D, 
		nel_ta n, 
		real *vec, 
		real var[2], 
		real wc, 
		agf_diag_param *diag_param);

  //AGF calculation with optimal filter variance (experimental):
  template <class real>
  real agf_calc_pdf_opt(real **mat, 
		dim_ta D, 
		nel_ta n, 
		real *vec, 
		real *parm, 			//parameters (see main .cc file)
		nel_ta ntest, 
		real &err, 
		int opts=0);			//options (see main .cc file)

  //AGF calculation with optimal filter variance (experimental):
  template <class real>
  real agf_calc_pdf_opt2(real **mat, 
		dim_ta D, 
		nel_ta n, 
		real *vec, 
		double *parm, 			//parameters (see main .cc file)
		int ntrial,
		real &err); 

}

#endif

