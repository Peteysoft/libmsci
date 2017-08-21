#ifndef SOLVE_COND_PROB_H_INCLUDED
#define SOLVE_COND_PROB_H_INCLUDED

#include <gsl/gsl_linalg.h>

using namespace libpetey;

namespace libagf {

  int solve_cond_prob2(gsl_matrix *a,		//decision matrix
		gsl_vector *r,			//original "raw" probabilities
		gsl_vector *p,			//returned cond. prob.
		int *gind,			//current active columns
		int &ng,			//number of active columns
		int iter);			//iteration

  int solve_cond_prob(gsl_matrix *a,		//decision matrix
		gsl_vector *r,			//original "raw" probabilities
		gsl_vector *p);			//returned cond. prob.

  //nothing fancy: 1. apply constraint sum_i p_i=1 by removing one variable
  //		and modifying decision matrix and solution vector appropriately
  //		2. if p_i dips below zero, remove ith column
  //		3. if p_i rises above 1, return p_i=1, all others 0
 

  //naive method works--produces results within the constraints 
  //doesn't improve classification accuracy
  //but does improve accuracy of probability estimates!
  template <class real, class cls_t>
  int derive_probability_map(real **r,		//binary probabilities
		cls_t *truth,			//true classes
		cls_t nmodel,			//number of binary classifiers
		cls_t ncls,			//number of classes
		nel_ta n,			//number of samples
		real (*skill)(cls_t *, cls_t *, nel_ta),
		real **map
		);

} //end namespace libagf

#endif

