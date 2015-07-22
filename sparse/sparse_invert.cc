#include "matrix_base.h"
#include "sparse_solve.h"

//normal equation and general inverse

namespace libpetey {
namespace libsparse {

template <class index_t, class real, class sparse_t>
int sparse_normal(sparse_t &a, real *b, real *x, void *parm,
		int (* solver) (sparse_t &, real *, real *, void *)
		) {

  index_t m, n;
  sparse_t *at;
  sparse_t *a2;
  real *atb;
  int err;

  a.dimensions(m, n);

  at=a.clone();
  at->transpose();

  a2=at->mat_mult(&a);

  atb=at->vect_mult(b);

  err=(*solver) (*a2, atb, x, parm);

  delete at;
  delete a2;
  delete [] atb;

  return err;
}

//find the inverse through simple, repeated application of a solver;
//returns transpose of inverse
template <class index_t, class real, class sparse_t>
int sparse_invert(sparse_t &a, real **ait, void *parm,
		int (* solver) (sparse_t &, real *, real *, void *)
		) {
  index_t m, n;
  index_t dmin;
  real *b;
  int err;

  a.dimensions(m, n);

  b=new real[m];

  for (index_t i=0; i<m; i++) {
    for (index_t j=0; j<m; j++) b[j]=0;
    b[i]=1;
    err=(*solver) (a, b, ait[i], parm);
  }

  delete [] b;

}

template int sparse_normal<int32_t, float, matrix_base<int32_t, float> >
		(matrix_base<int32_t, float> &, float *, float *, void *,
		int (* solver) (matrix_base<int32_t, float> &, 
		float *, float *, void *));

template int sparse_invert<int32_t, float, matrix_base<int32_t, float> >
		(matrix_base<int32_t, float> &, float **, void *,
		int (* solver) (matrix_base<int32_t, float> &, float *, float *, void *));

} //end namespace libsparse
} //end namespace libpetey

