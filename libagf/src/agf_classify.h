//
// This software is released under the following terms:
//
// 1. No commercial use.
// 2. Copies and derivative works are free to use and modify.
// 3. Attribution must be given to all contributors of both original and derivative works.
//
// Authors:
//
// 2017-07-16 Peter Mills: added license information 
//

//
// Direct AGF classification.
//

#ifndef AGF_CLASSIFY_H_INCLUDED
#define AGF_CLASSIFY_H_INCLUDED 1

#include "agf_defs.h"

namespace libagf {

  //for straight classification:
  template <class real, class cls_t>
  cls_t agf_classify(real **mat,	//matrix of training data
                dim_ta D,		//number of dimensions
                cls_t *cl,		//class labels
                nel_ta n,		//number of samples
                cls_t ncls,		//number of class types
                real *vec,		//test point
                real var[2],		//filter variance bracket
                nel_ta k,		//number of nearest neighbours to use
                real Wc,		//objective total weight
                real *pdf,		//returned joint/ cond. prob.
                agf_diag_param *diag_param,	//returned diagnostics
                flag_a joint=0);	//joint or conditional probabilities?

  //uses all of the data:
  template <class real, class cls_t>
  cls_t agf_classify(real **mat, 
		dim_ta D, 
		cls_t *cl, 
		nel_ta n, 
		cls_t ncl, 
		real *vec, 
		real var[2], 
		real wc, real *pdf, 
		agf_diag_param *diag_param, 
		flag_a joint=0);

}

#endif

