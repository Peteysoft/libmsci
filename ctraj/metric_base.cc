#include "metric_base.h"

namespace ctraj {

  template <class real>
  Cart_t<real>::Cart_t(int nd) {
    n=nd;
  }

  template <class real>
  real Cart_t<real>::ds2(real *x1, real *x2) {
    real result;
    real diff;
    result=0;
    for (int i=0; i<n; i++) {
      diff=x2[i]-x1[i];
      result+=diff*diff;
    }

    return result;
  }

  template <class real>
  void Cart_t<real>::mcoef2(real *x, real *c) {
    for (int i=0; i<n; i++) c[i]=1;
  }

  template <class real>
  real Cart_t<real>::gcurv(real *x) {
    return 0;
  }

  template <class real>
  int Cart_t<real>::ndim() {
    return n;
  }

  template class Cart_t<float>;
  template class Cart_t<double>;

} //end namespace ctraj


