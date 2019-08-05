#include <assert.h>

#include "matrix_base.h"

#include "av.h"

using namespace std;
using namespace libpetey;
using namespace libsparse;

void FORTRAN_FUNC(av)(int32_t *nx, float * v, float * w, void *parm) {
  matrix_base<int32_t, float> *A;
  int32_t m, n;

  A=(matrix_base<int32_t, float> *) parm;

  A->dimensions(m, n);
  //printf("matrix dimensions: %d %d\n", m, *nx);
  assert(*nx==n);
  assert(m==n);

  A->vect_mult(v, w);
  //printf("*\n");

}

void FORTRAN_FUNC(dav)(int32_t *nx, double * v, double * w, void *parm) {
  matrix_base<int32_t, double> *A;
  int32_t m, n;

  A=(matrix_base<int32_t, double> *) parm;

  A->dimensions(m, n);
  //printf("matrix dimensions: %d %d\n", m, *nx);
  assert(*nx==n);
  assert(m==n);

  A->vect_mult(v, w);
  //printf("*\n");

}

namespace libpetey {
namespace libsparse {

//C++ glue functions to glue the Fortran wrapper functions...
float **                         //returned eigenvectors
        cc_arsvd(int32_t n,             //matrix side length
                int32_t nev,                    //number of eigenvalues
                int32_t ncv,                    //number of Arnoldi vectors
                float *eval,                    //returned eigenvalues
                matrix_base<int32_t, float> *parm, //matrix
		const char which[2],
		int32_t maxitr,
		float tol)
		
{
  float ev[nev];
  float **v;
  v=new float*[ncv];
  v[0]=new float[ncv*n];
  for (int32_t i=1; i<ncv; i++) v[i]=v[0]+i*n;
  FORTRAN_FUNC(sarsvd)(&n, &nev, &ncv, v[0], ev, &maxitr, &tol, which, 
		(void *) parm);
  for (int32_t i=0; i<nev; i++) eval[i]=ev[i];
  return v;
}

double ** cc_arsvd(int32_t n,
        int32_t nev, int32_t ncv, double *eval,
        matrix_base<int32_t, double> *parm,
	const char which[2], int32_t maxiter, double tol) {
  double ev[nev];
  double **v;
  v=new double*[ncv];
  v[0]=new double[ncv*n];
  for (int32_t i=1; i<ncv; i++) v[i]=v[0]+i*n;
  FORTRAN_FUNC(darsvd)(&n, &nev, &ncv, v[0], ev, &maxiter, &tol, which,
		(void *) parm);
  for (int32_t i=0; i<nev; i++) eval[i]=ev[i];
  return v;
}

complex<float> ** cc_arevd(int32_t n,
        int32_t nev, int32_t ncv, complex<float> *eval,
        matrix_base<int32_t, float> *parm,
	const char which[2], int32_t maxitr, float tol) {
  float ev0[nev*2];
  float v0[ncv*n];
  complex<float> **v;
  int32_t k1, k2;

  FORTRAN_FUNC(sarevd)(&n, &nev, &ncv, v0, ev0, &maxitr, &tol, which,
		(void *) parm);

  v=new complex<float>*[nev];
  v[0]=new complex<float>[nev*n];
  for (int32_t i=0; i<nev; i++) {
    v[i]=v[0]+i*n;
    eval[i]=complex<float>(ev0[i], ev0[i+nev]);
  }

  //lets hope this works...
  for (int32_t i=0; i<nev; i++) {
    if (imag(eval[i])!=0) {
      for (int32_t j=0; j<n; j++) {
        k1=i*n+j;
        k2=(i+1)*n+j;
        v[i][j]=complex<float>(v0[k1], v0[k2]);
        v[i+1][j]=conj(v[i][j]);
      }
      i++;
    } else {
      for (int32_t j=0; j<n; j++) {
        k1=i*n+j;
        v[i][j]=v0[k1];
      }
    }
  }
  return v;
}

complex<double> ** cc_arevd(int32_t n,
        int32_t nev, int32_t ncv, complex<double> *eval,
        matrix_base<int32_t, double> *parm,
	const char which[2], int32_t maxiter, double tol) {
  double ev0[nev*2];
  double v0[ncv*n];
  complex<double> **v;
  int32_t k1, k2;

  FORTRAN_FUNC(darevd)(&n, &nev, &ncv, v0, ev0, &maxiter, &tol, which, 
		(void *) parm);

  v=new complex<double>*[nev];
  v[0]=new complex<double>[nev*n];
  for (int32_t i=0; i<nev; i++) {
    v[i]=v[0]+i*n;
    eval[i]=complex<float>(ev0[i], ev0[i+nev]);
  }

  //lets hope this works...
  for (int32_t i=0; i<nev; i++) {
    if (imag(eval[i])!=0) {
      for (int32_t j=0; j<n; j++) {
        k1=i*n+j;
        k2=(i+1)*n+j;
        v[i][j]=complex<double>(v0[k1], v0[k2]);
        v[i+1][j]=conj(v[i][j]);
      }
      i++;
    } else {
      for (int32_t j=0; j<n; j++) {
        k1=i*n+j;
        v[i][j]=v0[k1];
      }
    }
  }
  return v;
}

} //end namespace libsparse
} //end namespace libpetey
