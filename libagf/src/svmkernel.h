#ifndef _LIBAGF__SVMKERNEL__H_
#define _LIBAGF__SVMKERNEL__H_

#include "agf_defs.h"

namespace libagf {
  //basis functions ("kernels"):
  template <class real>
  real linear_basis(real *x, real *y, dim_ta n, void *param);

  template <class real>
  real polynomial_basis(real *x, real *y, dim_ta n, void *param);

  template <class real>
  real radial_basis(real *x, real *y, dim_ta n, void *param);

  template <class real>
  real sigmoid_basis(real *x, real *y, dim_ta n, void *param);

  template <class real>
  real linear_basis_deriv(real *x, real *y, dim_ta n, void *param, real *deriv);

  template <class real>
  real polynomial_basis_deriv(real *x, real *y, dim_ta n, void *param, real *deriv);

  template <class real>
  void radial_basis_deriv_num(real *x, real *y, dim_ta n, real dx, void *param, real *deriv);

  template <class real>
  real radial_basis_deriv(real *x, real *y, dim_ta n, void *param, real *deriv);

  template <class real>
  real sigmoid_basis_deriv(real *x, real *y, dim_ta n, void *param, real *deriv);

}

#endif

