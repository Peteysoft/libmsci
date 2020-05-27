#ifndef _AGF_BALLOON__H_
#define _AGF_BALLOON__H_

#include "agf_defs.h"

#include "balloon_param.h"

namespace libagf {
  //classify a test point using a variable bandwidth kernel density
  //estimator with constant weight:
  template <typename real, typename cls_t>
  cls_t balloon_classify(
		  balloon_param<real, long> *wt_calc,	//classifier parameters
		  nel_ta n,			//number of training samples
		  real **x,			//coordinates
		  cls_t ncls,			//number of classes
		  cls_t *cls,			//discrete ordinates [0-ncls]
		  real *test,			//test point
		  real *p);			//returned cond. prob.

  template <typename real, typename cls_t>
  cls_t balloon_classify_joint(balloon_param<real, long> *wt_calc,
		  nel_ta n,
		  real **x,
		  cls_t ncls,
		  cls_t *cls,
		  real *test,
		  real *p);		//returned joint probabilities

  //variable bandwidth kernel density estimator with constant weight:
  template <typename real>
  real balloon_density(balloon_param<real, long> *wt_calc,
		  nel_ta n,
		  real **x,
		  real *test);

  template <typename real>
  real balloon_interpol(balloon_param<real, long> *wt_calc,
		  nel_ta n,
		  real **x,
		  real *y,
		  real *test,
		  real *e2);

} //end namespace libagf

#endif
