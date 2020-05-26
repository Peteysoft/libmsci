#ifndef _KNN_PARAM__H_
#define _KNN_PARAM__H_

#include "balloon_param.h"

namespace libagf {

  template <typename real, typename index_t>
  class knn_param:public balloon_param<real, index_t> {
    private:
      index_t k;				//number of weight to calc.
    public:
      knn_param(real (*) (real), index_t, real, index_t);
      ~knn_param();

      virtual void n(index_t *);	//max size needed for weights

      //compute the weights, not including derivs:
      virtual index_t operator () (real *d2, index_t n, index_t *ind, real *w,
		      real *norm=NULL);
  };

} //end namespace libagf

#endif
