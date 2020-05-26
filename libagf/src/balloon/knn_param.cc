#include <stdint.h>
#include <math.h>
#include <gsl/gsl_sf_gamma.h>

#include "quicksort.h"

#include "knn_param.h"

using namespace libpetey;

namespace libagf {

  template <typename real, typename index_t>
  knn_param<real, index_t>::knn_param(index_t kay, index_t nvar) {
    k=kay;
    this->W=kay;
    this->D=nvar;
  }

  template <typename real, typename index_t>
  knn_param<real, index_t>::~knn_param() {
  }

  template <typename real, typename index_t>
  void knn_param<real, index_t>::maxn(index_t *mn) {
    *mn=k+1;
  }

  template <typename real, typename index_t>
  index_t knn_param<real, index_t>::operator () (
		real *d2, 				//distances
		index_t n,				//number of distances
		index_t *ind,				//returned indices
		real *w,				//returned weights
		real *norm) {				//returned norm. coef.
    real kleast[k+1];
    index_t D=this->D;
    real V;
    index_t extra=0;

    if (norm != NULL) extra=1;

    kleast_quick(d2, n, k+extra, kleast, ind);
    for (index_t i=0; i<k+extra; i++) w[i]=1;
    if (norm != NULL) {
      real sqrtpir=sqrt(M_PI)*(kleast[k-1]+kleast[k])/2;
      V=1.
      for (index_t i=0; i<D; i++) V*=sqrtpir;
      *norm=V/gsl_sf_gamma(D/2.+1);
    }
    return k;
  }

  template class knn_param<float, int32_t>
  template class knn_param<double, int32_t>
} //end namespace libagf

