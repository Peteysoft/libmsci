#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

#include "error_codes.h"

namespace libpetey {

  int solver(gsl_matrix *a,
		gsl_vector *b,
		gsl_vector *x) {
    int m, n;
    gsl_matrix *vt;
    gsl_vector *s;
    gsl_vector *work;

    if (a->size1 != b->size ||
		    a->size2 != x->size) {
      fprintf(stderr, "solver: dimension mismatch\n");
      return DIMENSION_MISMATCH;
    }
    if (a->size1 < a->size2) {
      m=a->size2;
      n=a->size1;
      gsl_matrix *u=gsl_matrix_alloc(m, n);
      gsl_matrix_transpose_memcpy(u, a);
    } else {
      m=a->size1;
      n=a->size2;
      gsl_matrix *u=gsl_matrix_alloc(a->size1, a->size2);
      gsl_matrix_memcpy(u, a);
    }

    vt=gsl_matrix_alloc(n, n);
    s=gsl_matrix_alloc(n);
    work=gsl_matrix_alloc(n);

    gsl_linalg_SV_decomp(u, vt, s, work);
    gsl_vector_free(work);

    if (a->size1 < a->size2) {
      //under-determined case,
      //just pretend it's a square problem:
      u->size1=n;
      gsl_linalg_SV_solve(vt, u, s, b, x);
      gsl_matrix_free(up);
    } else {
      //over-determined or square case:
      gsl_linalg_SV_solve(u, vt, s, b, x);
    }

    gsl_matrix_free(u);
    gsl_matrix_free(vt);
    gsl_vector_free(s);

    return 0;

  }

  int find_interior(gsl_matrix *a,	//matrix to solve
		gsl_vector *b,		//solution vector
		gsl_matrix *v,		//constraint normals
		gsl_vector *c,		//constraint thresholds
    		gsl_vector *interior,	//interior point
		double offset) {	//offset (if applicable)
    int err;
    gsl_vector *xp;		//interior point in transformed coords

    if (a->size1 != b->size ||
		    a->size2 != x->size ||
		    v->size2 != a->size2 ||
		    v->size1 > a->size2+1 ||
		    c->size != v->size1 ||
		    interior->size!= a->size2) {
      fprintf(stderr, "find_interior: dimension mismatch\n");
      return DIMENSION_MISMATCH;
    }

    //all the constraint thresholds must be positive:
    for (int i=0; i<c->size; i++) {
      double c_el=gsl_vector_get(c, i);
      if (c_el < 0) {
        gsl_vector_set(c, i, -c_el);
	for (int j=0; j<v->size2; j++) {
          gsl_matrix_set(v, i, j, -gsl_matrix_get(v, i, j));
	}
      }
    }

    xp=gsl_vector_alloc(c->size);

    if (c->size==a->size2+1) {
      gsl_vector *ubound=gsl_vector_alloc(c->size);
      solver(v, c, ubound);
      for (int i=0; i<c->size; i++) {
        double ub_el=gsl_vector_get(ubound, i);
	double c_el=gsl_vector_get(c, i);
	if (ub_el < c_el) {
          fprintf(stderr, "get_interior: contraints not consistent\n");
	  gsl_vector_free(ubound);
	  gsl_vector_free(xp);
	  return PARAMETER_OUT_OF_RANGE;
	}
        gsl_vector_set(xp, i, (c_el+ub_el)/2);
      }
      gsl_vector_free(ubound);
    } else {
      for (int i=0; i<c->size; i++) {
	double c_el=gsl_vector_get(c, i);
        gsl_vector_set(xp, i, c_el+offset);
      }
    }

    err=gsl_blas_dgemv(CblasNoTrans, 1., v, xp, 0., interior);

    return err;

  }

  void check_constraints(gsl_matrix *v,
		  gsl_vector *c,
		  gsl_vector *x,
		  int *flag) {
    if (v->size1 != c->size || v->size2 != x->size) {
      fprintf(stderr, "check constraints: dimension mismatch\n");
      exit(DIMENSION_MISMATCH);
    }
    for (int i=0; i<v->size2; i++) flag[i]=0;
    for (int i=0; i<v->size1; i++) {
      double *vrow;
      double dot;
      vrow=gsl_matrix_ptr(v, i, 0);
      dot=cblas_ddot(v->size2, vrow, 1, x->data, x->stride);
      if (dot<gsl_vector_get(c, i)) flag[i]=1;
    }
  }

  //returns the index of the eliminated variable:
  int apply_constraint(gsl_matrix *a, 
		gsl_vector *b,		//solution vector
		gsl_matrix *v,		//constraint normals
		gsl_vector *c,		//constraint thresholds
		int cindex,		//index of constraint
		gsl_matrix *a2,
		gsl_matrix *b2,
		gsl_matrix *v2,
		gsl_matrix *c2) {
    int vindex;		//index of eliminated variable
    for (int i=0; i<a->size2; i++) {
      if (gsl_matrix_get(v, i)!=0) {
        vindex=i;
	break;
      }
    }
    for (int i=0; i<vindex; i++) {
  }

  int constrained(gsl_matrix *a,	//matrix to solve
		gsl_vector *b,		//solution vector
		gsl_matrix *v,		//constraint normals
		gsl_vector *c,		//constraint thresholds
    		gsl_vector *interior,	//interior point
		gsl_vector *x){		//result

    gsl_vector *xtrial;		//unconstrained solution
    gsl_vector *interior;	//interior point
    gsl_matrix *

    if (a->size1 != b->size ||
		    a->size2 != x->size ||
		    v->size2 != a->size2 ||
		    v->size1 > a->size2+1 ||
		    c->size != v->size1 ||
		    interior->size!= a->size2) {
      fprintf(stderr, "constrained: dimension mismatch\n");
      return DIMENSION_MISMATCH;
    }

    //first, we find the unconstrained solution:


