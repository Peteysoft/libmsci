#ifndef _BALLOON_PARAM__H_
#define _BALLOON_PARAM__H_

#include <math.h>

namespace libagf {

  template <typename real>
  real gaussian_kernel(real r2, real *dfdr);

  template <typename real, typename index_t>
  class balloon_param {
    private:
      //kernel function:
      real (* kernel)(real, real *);		//kernel function
      index_t k;				//number of weight to calc.
      real var[2];				//variance brackets
      index_t maxiter;
      real tol;
    protected:
      index_t D;				//number of features
      real W;					//total weight
    public:
      balloon_param();
      balloon_param(
		      real (*K) (real, real*),	//kernel function
		      index_t nvar,		//number of features
		      real total_weight, 	//total of weights
		      index_t k2,		//number to use in calculation
		      real v1,			//lower variance bracket
		      real v2,			//upper variance bracket
		      real tolerance,		//tolerance (should have default)
		      index_t maxN);		//maximum iterations
      ~balloon_param();

      virtual void maxn(index_t *);	//max size needed for weights
      virtual index_t dim();		//number of features

      //compute the weights, not including derivs:
      virtual index_t operator () (real *d2, index_t n, index_t *ind, real *w,
		      real *norm=NULL);
  };

} //end namespace libagf

#endif
