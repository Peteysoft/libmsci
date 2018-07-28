#ifndef __LIBPETEY__CONSTRAINED__H
#define __LIBPETEY__CONSTRAINED__H

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

namespace libpetey {

  //solve under-determined least-squares problem using SVD:
  int SV_solve_transpose(gsl_matrix *vt,
		  gsl_matrix *u,
		  gsl_vector *s,
		  gsl_vector *b,
		  gsl_vector *x);

  //solve linear least-squares problem using SVD:
  int gsl_lsq_solver(gsl_matrix *a,
		gsl_vector *b,
		gsl_vector *x);

  int find_interior(gsl_matrix *v,	//constraint normals
		gsl_vector *c,		//constraint thresholds
    		gsl_vector *pp,		//interior point
		double offset,		//how far away from constraint borders
  					//to make interior point
		int (* solver1) (gsl_matrix *,
				gsl_vector *,
				gsl_vector *)=&gsl_lsq_solver);

  void normalize_constraints(gsl_matrix *v, gsl_vector *c);

  int check_constraints(gsl_matrix *v,		//matrix of constraint normals
		  gsl_vector *c,		//vector of constraint thresholds
		  gsl_vector *x,		//point to test
		  int *ind);			//list of violated constraints

  int test_interior(int n);

  int constrain_row(gsl_vector *row,		//row to transform
		  gsl_vector *constraint, 	//constraint normal
		  int index,			//index of variable to eliminate
		  gsl_vector *row2);		//transformed row

  int apply_constraint(gsl_matrix *v,		//constraint normals
		gsl_vector *c,		//constraint thresholds
		int cindex,		//index of constraint
		int vindex,		//index of variable to eliminate
		gsl_matrix *v2,
		gsl_vector *c2);
 
  int constrained(gsl_matrix *v,	//transformed constraint normals
		gsl_vector *c,		//transformed constraint thresholds
		gsl_vector *p,		//interior point
		gsl_vector *x);		//solution

  int constrained(gsl_matrix *a,	//matrix to solve
		gsl_vector *b,		//solution vector
		gsl_matrix *v,		//constraint normals
		gsl_vector *c,		//constraint thresholds
		gsl_vector *p,		//interior point
		gsl_vector *x,		//solution vector
		int (* solver1) (gsl_matrix *,
				gsl_vector *,
				gsl_vector *)=&gsl_lsq_solver);

  int constrained(gsl_matrix *a,	//matrix to solve
		gsl_vector *b,		//solution vector
		gsl_matrix *v,		//constraint normals
		gsl_vector *c,		//constraint thresholds
		gsl_vector *x,		//result
		int (* solver1) (gsl_matrix *,
				gsl_vector *,
				gsl_vector *)=&gsl_lsq_solver);
}

#endif

