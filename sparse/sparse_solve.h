#ifndef SPARSE_SOLVE__H
#define SPARSE_SOLVE__H 1

#include <stdint.h>

#include "av.h"

namespace libpetey {
  namespace libsparse {

    template <class matrix>
    struct sparse_solver_parm {
      int maxiter;			//maximum number of iterations
      double tol;			//desired tolerance
      int nev;			//number of eigenvectors/values
      int ncv;			//number of Arnoldi vectors
      matrix *cond;			//conditioning matrix
      double err;			//returned error
    };

    template <class index_t, class value, class matrix>
    value conj_grad(matrix &a, value *b, value *x, value tol, index_t maxiter);

    template <class index_t, class value, class matrix>
    value biconj_grad(matrix &A, matrix &P, value *b, value *x, value tol, index_t maxiter);

    template <class index_t, class value, class matrix>
    int conj_grad(matrix &a, value *b, value *x, void *parm);

    template <class index_t, class value, class matrix>
    int biconj_grad(matrix &A, value *b, value *x, void *parm);

    template <class index_t, class real, class sparse_t>
    int sparse_svd_solver(sparse_t &A, real *b, real *x, int32_t nev, int32_t ncv,
		int32_t maxiter=DEF_MAXNITER, real tol=DEF_TOL);

    template <class index_t, class real, class sparse_t>
    int sparse_evd_solver(sparse_t &A, real *b, real *x, int32_t nev, int32_t ncv,
		int32_t maxiter=DEF_MAXNITER, real tol=DEF_TOL);

    template <class index_t, class real, class sparse_t>
    int sparse_svd_solver(sparse_t &A, real *b, real *x, void *parm);

    template <class index_t, class real, class sparse_t>
    int sparse_evd_solver(sparse_t &A, real *b, real *x, void *parm);

    template <class index_t, class real, class sparse_t>
    real ** sparse_svd_invert(sparse_t &A, int32_t nev, int32_t ncv);

    template <class index_t, class real, class sparse_t>
    real ** sparse_evd_invert(sparse_t &A, int32_t nev, int32_t ncv);

    template <class index_t, class real, class sparse_t>
    real ** sparse_svd_invert(sparse_t &A, void *parm);

    template <class index_t, class real, class sparse_t>
    real ** sparse_evd_invert(sparse_t &A, void *parm);

    //normal equation and general inverse

    template <class index_t, class real, class sparse_t>
    int sparse_normal(sparse_t &a, real *b, real *x, void *parm,
		int (* solver) (sparse_t &, real *, real *, void *)
		);

    //find the inverse through simple, repeated application of a solver;
    //returns transpose of inverse
    template <class index_t, class real, class sparse_t>
    int sparse_invert(sparse_t &a, real **ait, void *parm,
		int (* solver) (sparse_t &, real *, real *, void *)
		);
  } //end namespace libsparse
} //end namespace libpetey

#endif

