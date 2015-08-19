#ifndef __LIBPETEY__CONSTRAINED__H
#define __LIBPETEY__CONSTRAINED__H

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

namespace libpetey {

  int constrained(gsl_matrix *a,	//matrix to solve
		gsl_vector *b,		//solution vector
		gsl_matrix *v,		//constraint normals
		gsl_vector *c,		//constraint thresholds
		gsl_vector *x);		//result
}

#endif

