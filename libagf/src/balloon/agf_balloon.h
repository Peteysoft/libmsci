#ifndef _AGF_BALLOON__H_
#define _AGF_BALLOON__H_

#include "balloon_param.h"

namespace libagf {
  //classify a test point using a variable bandwidth kernel density
  //estimator with constant weight:
  template <typename real, typename cls_t, typename index_t>
  cls_t balloon_classify(
		  balloon_param *wt_calc,	//classifier parameters
		  index_t n,			//number of training samples
		  real **x,			//coordinates
		  cls_t ncls,			//number of classes
		  cls_t *cls,			//discrete ordinates [0-ncls]
		  real *test,			//test point
		  real *p);			//returned cond. prob.

  template <typename real, typename cls_t, typename index_t>
  cls_t balloon_classify_joint(balloon_param *wt_calc,
		  index_t n,
		  real **x,
		  cls_t *cls,
		  cls_t ncls,
		  real *test,
		  real *p);		//returned joint probabilities

  //variable bandwidth kernel density estimator with constant weight:
  template <typename real, typename index_t>
  real balloon_density(balloon_param *wt_calc,
		  index_t n,
		  real **x,
		  real *test);

  template <typename real, typename index_t>
  real balloon_interpol(balloon_param *wt_calc,
		  index_t n,
		  real **x,
		  real *y,
		  real *test,
		  real *e2);

} //end namespace libagf

