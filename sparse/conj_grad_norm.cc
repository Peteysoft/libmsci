#include <assert.h>
#include "sparse.h"

namespace libpetey {
namespace libsparse {

//solves the normal equation for a sparse matrix using the conjugate gradient method:
void conj_grad_norm(sparse_matrix *a, float *x, float *b, float tol, long maxiter) {
  ind_type m, n;
  sparse_matrix a_t, a1;
  ind_type nx, nb;
  float *b1, *b2;
  float b3, v;
  float divisor, dividend;
  float err, old_err;
  float tol1;
  ind_type this_i, this_j;
  float val;

  a->mat_mult_t(*a, a1);
  a->dimensions(nx, nb);	//a is the transpose of the matrix we want to solve
  a1.dimensions(m, n);
  assert(m==n);			//a^2 should be square
  assert(nx==m);		//outer dimension of a^2 should be same as

  tol1=tol*tol*nb;

  b1=new float[nb];		//result vector in the normal equation
  b2=new float[nb];		//result vector for minimizing the conjugate equation in one dimension

  a->vect_mult(b, b1);

  for (long i=0; i<maxiter; i++) {
    a1.vect_mult(x, b2);
    old_err=err;
    err=0;
    for (long k=0; k<n; k++) {
      b3=b1[k]-b2[k];
      err+=b3*b3;		//while we're at it, we might as well calculate the error
    }

    if (tol1 > err || (i>0 && err >= old_err)) {
      break;
    }

    this_i=a1.matrix[0].i;
    b3=b1[this_i];
    v=0;
    for (long j=0; j<a1.nel; j++) {
      //check for change in row index:
      if (a1.matrix[j].i != this_i) {
        //printf("%d: v: %g, b3: %g\n", this_i, v, b3);
        //update ith value of x:
	if (v!=0) {
          x[this_i]+=b3/v;
	} else {
	  x[this_i]=0;
	}
        this_i=a1.matrix[j].i;
        //reset other parameters:
        b3=b1[this_i];
        v=0;
      }
      this_j=a1.matrix[j].j;
      val=a1.matrix[j].value;
      b3-=x[this_j]*val;
      if (this_j==this_i) v=val;
    }

    printf("%5d %10.5f\n", i, err);
  }
}

}} //end namespace libpetey::libsparse

