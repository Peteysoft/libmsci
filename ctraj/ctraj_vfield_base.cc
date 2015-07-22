#include <stdio.h>

#include "ctraj_vfield_base.h"

namespace ctraj {

  //since not every class will have a working one of these:
  template <class real>
  int ctraj_vfield_base<real>::jmat(int32_t d, double t, real *x, real *j) {
  }

  //linker is looking for this for some reason:
  template <class real>
  ctraj_vfield_base<real>::~ctraj_vfield_base() {
  }

  template <class real>
  ctraj_vfield_single_domain<real>::ctraj_vfield_single_domain() {
  }

  template <class real>
  int32_t ctraj_vfield_single_domain<real>::fix(int32_t domain, double tind, real *x) {
    return 0;
  }

  template <class real>
  int32_t ctraj_vfield_single_domain<real>::absolute(int32_t domain, real *x) {
    if (domain < 0) return 0; 
    return -1;
  }

  template <class real>
  int32_t ctraj_vfield_single_domain<real>::reference(int32_t domain, real *x) {
    if (domain < 0) return 0;
    
    return 0;
  }

  template <class real>
  double ctraj_vfield_single_domain<real>::get_tind(char *date) {
    double tind;

    return sscanf(date, "%lg", &tind);
  }

  template <class real>
  int ctraj_vfield_single_domain<real>::get_t(double tind, char *date) {
    return sprintf(date, "%lg", tind);
  }

  template class ctraj_vfield_base<float>;
  template class ctraj_vfield_single_domain<float>;

} //end namespace ctraj

