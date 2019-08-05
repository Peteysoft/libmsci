#include <stdint.h>
#include <math.h>
#include <assert.h>

#include "full_util.h"
#include "matrix_base.h"
#include "sparse_solve.h"

//SVD decomposition and eigenvalue decomposition (for Hermitian matrices)

namespace libpetey {
namespace libsparse {

template <class index_t, class real, class sparse_t>
int sparse_svd_solver(sparse_t &A, real *b, real *x, int32_t nev, int32_t ncv, 
		int32_t maxiter, real tol) {
  index_t m, n;
  sparse_t *AT;
  sparse_t *A2;
  //right singular vectors:
  real **v;
  real *s;		//singular values
  real *eval1;
  real *xp, *y;		//intermediate vectors
  index_t k1, k2;

  A.dimensions(m, n);

  eval1=new real[nev];
  s=new real[nev];

  //get the right singular vectors:
  AT=A.clone();
  AT->transpose();
  A2=AT->mat_mult(&A);

  v=cc_arsvd(n, nev, ncv, eval1, A2, "LM", maxiter, tol);

  // x = V S^-2 V^T A^T b
  xp=A.vect_mult(b);
  y=vector_mult(v, xp, nev, n);
  for (index_t i=0; i<nev; i++) {
    printf("%g\n", eval1[i]);
    y[i]/=eval1[i];
  }
  left_vec_mult(y, v, nev, n);

  delete [] xp;
  delete [] y;
  delete [] v[0];
  delete [] v;
  delete [] eval1;
  delete [] s;
  delete AT;
  delete A2;

  return 0;

}

//this only works for Hermitian matrices:
template <class index_t, class real, class sparse_t>
int sparse_evd_solver(sparse_t &A, real *b, real *x, int32_t nev, int32_t ncv, 
		int32_t maxiter, real tol) {
  using namespace std;

  index_t m, n;
  complex<real> **u;
  complex<real> **v;
  complex<real> *eval;
  complex<real> sum;

  A.dimensions(m, n);
  //assert(m==n);

  //v=new real[m*nev];

  eval=new complex<real>[nev];

  //FORTRAN_FUNC(ev_sparse_array)(&n, &nev, &ncv, v, eval, (void *) &A);
  u=cc_arevd(m, nev, ncv, eval, &A, "LM", maxiter, tol);
  A.transpose();
  v=cc_arevd(n, nev, ncv, eval, &A, "LM", maxiter, tol);

  //this works exactly the same as in the above routine:
  //apparently, sneupd can return an orthonormal basis in which case the
  //following will work
  //presumeably we can also do something similar with left and right 
  //eigenvectors
  //in either case, I don't know how to extract the imaginary parts...
  //perhaps simplest just to use SVD and not bother with this at all
  //(except we still need it for solving linear systems of ODEs...)
  for (index_t i=0; i<m; i++) {
    x[i]=0;
    for (index_t j=0; j<nev; j++) {
      sum=0;
      for (index_t k=0; k<n; k++) {
        sum+=v[j][k]*b[k];
      }
      //result should be real (at least in theory): 
      x[i]=x[i]+std::real(sum*u[j][i]/eval[j]);
    }
  }

  delete [] v[0];
  delete [] v;
  delete [] u[0];
  delete [] u;
  delete [] eval;

  //doesn't always succeed, but still, need to return sthng...

  return 0;

}

template <class index_t, class real, class sparse_t>
int sparse_svd_solver(sparse_t &A, real *b, real *x, void *parm) {
  sparse_solver_parm<sparse_t> *p1;
  p1=(sparse_solver_parm<sparse_t> *) parm;
  return sparse_svd_solver<index_t, real, sparse_t>
		(A, b, x, (index_t) p1->nev, (index_t) p1->ncv, 
		p1->maxiter, p1->tol);
}

template <class index_t, class real, class sparse_t>
int sparse_evd_solver(sparse_t &A, real *b, real *x, void *parm) {
  sparse_solver_parm<sparse_t> *p1;
  p1=(sparse_solver_parm<sparse_t> *) parm;
  return sparse_evd_solver<index_t, real, sparse_t>
		(A, b, x, (index_t) p1->nev, (index_t) p1->ncv,
		p1->maxiter, p1->tol);
}

template int sparse_svd_solver<int32_t, float, matrix_base<int32_t, float> >
	(matrix_base<int32_t, float> &, float *, float *, void *);

template int sparse_evd_solver<int32_t, float, matrix_base<int32_t, float> >
	(matrix_base<int32_t, float> &, float *, float *, void *);

//SVD decomposition and eigenvalue decomposition (for Hermitian matrices)

template <class index_t, class real, class sparse_t>
real ** sparse_svd_invert(sparse_t &A, int32_t nev, int32_t ncv,
		int32_t maxiter, real tol) {
  index_t m, n;
  sparse_t *AT;
  sparse_t *A2;
  real **u, **v;
  real *s;		//singular values
  real *eval1, *eval2;
  real sum;
  index_t k1, k2;
  real **result;

  A.dimensions(m, n);

  eval1=new real[nev];
  eval2=new real[nev];
  s=new real[nev];

  //get the right singular vectors:
  AT=A.clone();
  AT->transpose();
  A2=AT->mat_mult(&A);

  v=cc_arsvd(n, nev, ncv, eval1, A2, "LM", maxiter, tol);

  //get the right singular vectors:
  //(is this faster or better than just multiplying through with A?)
  //delete A2;
  //A2=AT->mat_mult(&A);
  //v=cc_arsvd(n, nev, ncv, eval2, A2);

  //multiply through to get the left sinular vectors:
  //get the singular values by averaging the eigenvalues:
  //(are they in the same order??)
  u=allocate_matrix<real, index_t>(nev, m);
  for (index_t i=0; i<nev; i++) {
    printf("%g\n", eval1[i]);
    A.vect_mult(v[i], u[i]);
    //s[i]=sqrt(eval1[i]);
    for (index_t j=0; j<m; j++) u[i][j]/=eval1[i];
  }

  result=allocate_matrix<real, index_t>(n, m);

  for (index_t i=0; i<n; i++) {
    for (index_t j=0; j<m; j++) {
      result[i][j]=0;
      for (index_t k=0; k<nev; k++) {
        result[i][j]+=u[k][j]*v[k][i];
      }
    }
  }

  delete [] u[0];
  delete [] u;
  delete [] v[0];
  delete [] v;
  delete [] eval1;
  delete [] eval2;
  delete [] s;
  delete AT;
  delete A2;

  return result;

}

//this only works for Hermitian matrices: (?? --probably doesn't work at all...)
template <class index_t, class real, class sparse_t>
real ** sparse_evd_invert(sparse_t &A, int32_t nev, int32_t ncv, 
		int32_t maxiter, real tol) {
  using namespace std;

  index_t m, n;
  complex<real> **u;
  complex<real> **v;
  complex<real> *eval;
  complex<real> sum;

  real **result;

  A.dimensions(m, n);
  //assert(m==n);

  //v=new real[m*nev];

  eval=new complex<real>[nev];

  //FORTRAN_FUNC(ev_sparse_array)(&n, &nev, &ncv, v, eval, (void *) &A);
  u=cc_arevd(m, nev, ncv, eval, &A, "LM", maxiter, tol);
  A.transpose();
  v=cc_arevd(n, nev, ncv, eval, &A, "LM", maxiter, tol);

  //this works exactly the same as in the above routine:
  //apparently, sneupd can return an orthonormal basis in which case the
  //following will work
  //presumeably we can also do something similar with left and right 
  //eigenvectors
  //in either case, I don't know how to extract the imaginary parts...
  //perhaps simplest just to use SVD and not bother with this at all
  //(except we still need it for solving linear systems of ODEs...)
  result=allocate_matrix<real, index_t>(n, m);
  for (index_t i=0; i<m; i++) {
    for (index_t k=0; k<n; k++) {
      result[k][i]=0;
      for (index_t j=0; j<nev; j++) {
        result[k][i]=result[k][i]+std::real(v[j][k]*u[j][i]/eval[j]);
      }
      //result should be real (at least in theory): 
    }
  }

  delete [] v[0];
  delete [] v;
  delete [] u[0];
  delete [] u;
  delete [] eval;

  return result;

}

template <class index_t, class real, class sparse_t>
real ** sparse_svd_invert(sparse_t &A, void *parm) {
  sparse_solver_parm<sparse_t> *p1;
  p1=(sparse_solver_parm<sparse_t> *) parm;
  return sparse_svd_invert<index_t, real, sparse_t>
		(A, (index_t) p1->nev, (index_t) p1->ncv,
		p1->maxiter, p1->tol);
}

template <class index_t, class real, class sparse_t>
real ** sparse_evd_invert(sparse_t &A, void *parm) {
  sparse_solver_parm<sparse_t> *p1;
  p1=(sparse_solver_parm<sparse_t> *) parm;
  return sparse_evd_invert<index_t, real, sparse_t>
		(A, (index_t) p1->nev, (index_t) p1->ncv,
		p1->maxiter, p1->tol);
}

template float ** sparse_svd_invert<int32_t, float, matrix_base<int32_t, float> >
	(matrix_base<int32_t, float> &, void *);

template float ** sparse_evd_invert<int32_t, float, matrix_base<int32_t, float> >
	(matrix_base<int32_t, float> &, void *);

} //end namespace libsparse
} //end namespace libpetey


