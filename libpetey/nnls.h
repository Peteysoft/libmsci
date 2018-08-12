#ifndef _NNLS_H_INCLUDED_
#define _NNLS_H_INCLUDED_

#include <stdint.h>
#include "av.h"

extern "C" {
	void FORTRAN_FUNC(ldp)(double **g, 
			int32_t *mdg, 
			int32_t *m, 
			double *h,
			double *x,
			double *xnorm,
			double *w,
			int32_t *index,
			int32_t *mode);

	void FORTRAN_FUNC(nnls)(double **a,
			int32_t *mda,
			int32_t *m,
			int32_t *n,
			double *b,
			double *x,
			double *rnorm,
			double *w,
			double *zz,
			int32_t *mode);

	void FORTRAN_FUNC(svdrs)(double **a, 
			int32_t *mda, 
			int32_t *M1, 
			int32_t *N1, 
			double **B, 
			int32_t *MDB, 
			int32_t *NB, 
			double *S, 
			double **WORK);

}

#endif

