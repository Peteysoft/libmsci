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
      real (* kernel) (real, real *);		//kernel function
      index_t k;				//number of weight to calc.
      real var[2];				//variance brackets
      index_t maxiter;
      real tol;
    protected:
      index_t D;				//number of features
      real W;					//total weight
    public:
      balloon_param(real (*) (real), index_t, real, index_t);
      ~balloon_param();

      virtual void n(index_t *);	//max size needed for weights
      virtual index_t dim();		//number of features

      //compute the weights, not including derivs:
      virtual index_t operator () (real *d2, index_t n, index_t *ind, real *w,
		      real *norm=NULL);
  };

} //end namespace libagf

#endif
