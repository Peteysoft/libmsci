#include <assert.h>
#include <complex>

#include "error_codes.h"

#include "matrix_base.h"
#include "sparse_solve.h"

using namespace std;

namespace libpetey {
namespace libsparse {

inline double conj(double x) {return x;}

//copied more or less verbatim from wikipedia:
template <class index_t, class value, class matrix>
value conj_grad(matrix &a, value *b, value *x, value tol, index_t maxiter) {
  index_t m, n;
  value *r;
  value *p;
  value *ap;
  value alpha;
  float rsold, rsnew;
  index_t i;

  tol=tol*tol;

  a.dimensions(m, n);
  assert(m == n);

  r=new value[n];
  p=new value[n];
  ap=new value[n];

  a.vect_mult(x, r);

  rsold=0;
  for (index_t j=0; j<n; j++) {
    r[j]=b[j]-r[j];
    p[j]=r[j];
    rsold+=conj(r[j])*r[j];
  }

  for (i=1; i<maxiter; i++) {
    a.vect_mult(p, ap);
    alpha=0;
    for (index_t j=0; j<n; j++) {
      alpha+=p[j]*ap[j];
    }
    alpha=rsold/alpha;
    rsnew=0;
    for (index_t j=0; j<n; j++) {
      x[j]=x[j]+alpha*p[j];
      r[j]=r[j]-alpha*ap[j];
      rsnew+=conj(r[j])*r[j];
    }
    //printf("conj_grad: err=%g\n", rsnew);
    if (rsnew < tol) break;
    for (index_t j=0; j<n; j++) p[j]=r[j]+rsnew*p[j]/rsold;
    rsold=rsnew;
  }

  //printf("conj_grad: tol=%g; err=%g; maxiter=%d; iter=%d\n", tol, rsnew, maxiter, i);

  delete [] r;
  delete [] p;
  delete [] ap;

  return rsnew;
}

//copied more or less verbatim from wikipedia:
template <class index_t, class value, class matrix>
value biconj_grad(matrix &A, matrix &P, value *b, value *x, value tol, index_t maxiter) {
  index_t m, n;
  index_t m1, n1;

  value *r, *ra;
  value *xa, *ba;
  value *p, *pa;

  value *rn, *rna;
  value *swap;

  value *Pr, *Ap;
  value *PTra, *ATpa;

  value alpha, beta;
  value alpha_c;
  value raPr_old, raPr_new;

  value err;

  A.dimensions(m, n);
  P.dimensions(n1, m1);

  if (m!=m1 || n!=n1) {
    fprintf(stderr, "biconj_grad: dimensions of multiplicand [%dx%d] and conditioner [%dx%d] must match\n", m, n, m1, n1);
    return -1;
  }

  tol=tol*tol;

  r=new value[m];
  ra=new value[n];

  xa=new value[m];
  ba=new value[n];

  p=new value[n];
  pa=new value[m];

  Pr=new value[n];
  Ap=new value[m];
  PTra=new value[m];
  ATpa=new value[n];

  for (index_t j=0; j<m; j++) xa[j]=b[j];
  for (index_t j=0; j<n; j++) ba[j]=x[j];

  A.vect_mult(x, r);
  for (index_t j=0; j<m; j++) r[j]=b[j]-r[j];
  A.left_mult(xa, ra);
  for (index_t j=0; j<n; j++) ra[j]=ba[j]-ra[j];

  P.vect_mult(r, Pr);
  P.left_mult(ra, PTra);
  A.vect_mult(p, Ap);
  A.left_mult(pa, ATpa);
  raPr_old=0;
  for (index_t j=0; j<n; j++) {
    p[j]=Pr[j];
    pa[j]=PTra[j];
    raPr_old+=ra[j]*Pr[j];
  }

  for (index_t i=0; i<maxiter; i++) {

    alpha=0;
    for (index_t j=0; j<m; j++) alpha+=pa[j]*Ap[j];
    alpha=raPr_old/alpha;

    alpha_c=conj(alpha);
    for (index_t j=0; j<n; j++) {
      x[j]=x[j]+p[j]*alpha;
      ra[j]=ra[j]-alpha_c*ATpa[j];
    }
    err=0;
    for (index_t j=0; j<m; j++) {
      //xa[j]=xa[j]+pa[j]*alpha_c;
      r[j]=r[j]-alpha*Ap[j];
      err+=conj(r[j])*r[j];
    }

    printf("biconj_grad: err=%g\n", err);
    if (err < tol) break;

    A.vect_mult(p, Ap);
    P.vect_mult(r, Pr);
    A.left_mult(pa, ATpa);
    P.left_mult(ra, PTra);

    raPr_new=0;
    for (index_t j=0; j<n; j++) raPr_new+=ra[j]*Pr[j];
    beta=raPr_new/raPr_old;

    for (index_t j=0; j<n; j++) p[j]=Pr[j]+beta*p[j];
    for (index_t j=0; j<m; j++) pa[j]=PTra[j]+conj(beta)*pa[j];

    raPr_old=raPr_new;
  }

  
  delete [] r;
  delete [] ra;

  delete [] xa;
  delete [] ba;

  delete [] p;
  delete [] pa;

  delete [] Pr;
  delete [] Ap;
  delete [] PTra;
  delete [] ATpa;

  return err;

}

template <class index_t, class value, class matrix>
int conj_grad(matrix &a, value *b, value *x, void *parm) {
  sparse_solver_parm<matrix> *p1;
  p1=(sparse_solver_parm<matrix> *) parm;

  p1->err=(double) conj_grad(a, b, x, (value) p1->tol, p1->maxiter);
  return (p1->err>p1->tol);
}

template <class index_t, class value, class matrix>
int biconj_grad(matrix &A, value *b, value *x, void *parm) {
  sparse_solver_parm<matrix> *p1;
  p1=(sparse_solver_parm<matrix> *) parm;
  int err;

  p1->err=(double) biconj_grad(A, *p1->cond, b, x, (value) p1->tol, p1->maxiter);

  if (p1->err < 0) err=DIMENSION_MISMATCH;
  if (p1->err > p1->tol) err=NUMERICAL_ERROR;

  return err;
}

template int conj_grad<int32_t, float, matrix_base<int32_t, float> >
		(matrix_base<int32_t, float> &, 
		float *, float *, void *);

template int biconj_grad<int32_t, float, matrix_base<int32_t, float> >
		(matrix_base<int32_t, float> &, 
		float *, float *, void *);

} } //end namespace libpetey::libsparse

