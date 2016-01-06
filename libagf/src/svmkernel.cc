#include <math.h>
#include <string.h>

#include "agf_defs.h"
#include "agf_metric2.h"

namespace libagf {
  //basis functions ("kernels"):
  template <class real>
  real linear_basis(real *x, real *y, dim_ta n, void *param) {
    real dot=0;
    for (dim_ta i=0; i<n; i++) dot+=x[i]*y[i];
    return dot;
  }

  template <class real>
  real polynomial_basis(real *x, real *y, dim_ta n, void *param) {
    real *p2=(real *) param;
    real gamma=p2[0];
    real coef0=p2[1];
    real degree=p2[2];
    real dot=linear_basis(x, y, n, NULL);
    return pow(gamma*dot+coef0, degree);
  }

  template <class real>
  real radial_basis(real *x, real *y, dim_ta n, void *param) {
    real gamma=((real *) param)[0];
    real d2=metric2(x, y, n);
    return exp(-gamma*d2);
  }

  template <class real>
  real sigmoid_basis(real *x, real *y, dim_ta n, void *param) {
    real *p2=(real *) param;
    real gamma=p2[0];
    real coef0=p2[1];
    real dot=linear_basis(x, y, n, NULL);
    return tanh(gamma*dot+coef0);
  }

  //basis functions ("kernels"):
  template <class real>
  real linear_basis_deriv(real *x, real *y, dim_ta n, void *param, real *deriv) {
    real dot=0;
    for (dim_ta i=0; i<n; i++) {
      dot+=x[i]*y[i];
      deriv[i]=y[i];
    }
    return dot;
  }

  template <class real>
  real polynomial_basis_deriv(real *x, real *y, dim_ta n, void *param, real *deriv) {
    real *p2=(real *) param;
    real gamma=p2[0];
    real coef0=p2[1];
    real degree=p2[2];
    real dot=linear_basis_deriv(x, y, n, NULL, deriv);
    for (dim_ta i=0; i<n; i++) {
      deriv[i]*=pow(gamma*dot+coef0, degree-1)*degree*gamma;
    }
    return pow(gamma*dot+coef0, degree);
  }

  template <class real>
  void radial_basis_deriv_num(real *x, real *y, dim_ta n, real dx, void *param, real *deriv) {
    real x1[n];
    real x2[n];
    real r1, r2;
    for (dim_ta i=0; i<n; i++) {
      x1[i]=x[i];
      x2[i]=x[i];
    }
    for (dim_ta i=0; i<n; i++) {
      x1[i]-=dx;
      x2[i]+=dx;
      r1=radial_basis(x1, y, n, param);
      r2=radial_basis(x2, y, n, param);
      deriv[i]=(r2-r1)/dx/2;
      x1[i]=x[i];
      x2[i]=x[i];
    }
  }

  template <class real>
  real radial_basis_deriv(real *x, real *y, dim_ta n, void *param, real *deriv) {
    real deriv1[n];
    real gamma=((real *) param)[0];
    real t1;
    real d2=0;
    for (dim_ta i=0; i<n; i++) {
      deriv[i]=x[i]-y[i];
      d2+=deriv[i]*deriv[i];
    }
    t1=exp(-gamma*d2);
    for (dim_ta i=0; i<n; i++) deriv[i]*=-2*t1*gamma;

    return t1;

    //sanity check (checks out...):
    radial_basis_deriv_num(x, y, n, real(0.001), param, deriv1);
    for (dim_ta i=0; i<n; i++) printf(" %g", deriv[i]);
    printf("\n");
    for (dim_ta i=0; i<n; i++) printf(" %g", deriv1[i]);
    printf("\n");
  }

  template <class real>
  real sigmoid_basis_deriv(real *x, real *y, dim_ta n, void *param, real *deriv) {
    real *p2=(real *) param;
    real gamma=p2[0];
    real coef0=p2[1];
    real dot=linear_basis_deriv(x, y, n, NULL, deriv);
    real t1=tanh(gamma*dot+coef0);
    for (dim_ta i=0; i<n; i++) deriv[i]*=gamma*(1-t1*t1);
    return t1;
  }

  template float linear_basis<float>(float *, float *, dim_ta, void *);
  template float polynomial_basis<float>(float *, float *, dim_ta, void *);
  template float radial_basis<float>(float *, float *, dim_ta, void *);
  template float sigmoid_basis<float>(float *, float *, dim_ta, void *);

  template double linear_basis<double>(double *, double *, dim_ta, void *);
  template double polynomial_basis<double>(double *, double *, dim_ta, void *);
  template double radial_basis<double>(double *, double *, dim_ta, void *);
  template double sigmoid_basis<double>(double *, double *, dim_ta, void *);

  template float linear_basis_deriv<float>(float *, float *, dim_ta, void *, float *);
  template float polynomial_basis_deriv<float>(float *, float *, dim_ta, void *, float *);
  template float radial_basis_deriv<float>(float *, float *, dim_ta, void *, float *);
  template float sigmoid_basis_deriv<float>(float *, float *, dim_ta, void *, float *);

  template double linear_basis_deriv<double>(double *, double *, dim_ta, void *, double *);
  template double polynomial_basis_deriv<double>(double *, double *, dim_ta, void *, double *);
  template double radial_basis_deriv<double>(double *, double *, dim_ta, void *, double *);
  template double sigmoid_basis_deriv<double>(double *, double *, dim_ta, void *, double *);

}
