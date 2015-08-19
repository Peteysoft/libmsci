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

  //
  int find_interior(
		gsl_matrix *v,		//constraint normals
		gsl_vector *c,		//constraint thresholds
    		gsl_vector *interior,	//interior point
		double offset) {	//offset (if applicable)
    int err;
    gsl_vector *pp;		//interior point in transformed coords

    if (v->size1 != c->size ||
		    interior->size!= a->size2) {
      fprintf(stderr, "find_interior: dimension mismatch\n");
      return DIMENSION_MISMATCH;
    }

    pp=gsl_vector_alloc(v->size2);

    if (c->size==v->size2+1) {
      gsl_vector *xp=gsl_vector_alloc(a->size2);
      double dot;		//dot product--must be less than c_n
      c->size--;
      v->size1--;
      solver(v, c, xp);
      dot=cblas_ddot(c->size, c, c->stride, xp->data, xp->stride);
      c->size++;
      v->size1++;
      if (dot > gsl_vector_get(c, c->size-1)) {
        fprintf(stderr, "get_interior: contraints not consistent\n");
        exit(PARAMETER_OUT_OF_RANGE);
      }
      for (int i=0; i<v->size2; i++) {
        double pp_el=gsl_vector_get(pp, i);
        gsl_vector_set(pp, pp_el+dot/pp_el/a->size2);
      }
      gsl_vector_free(xp);
    } else {
      fprintf(stderr, "get_interior: case for less than n+1 contraints not implemented yet, sorry.\n");
      exit(PARAMETER_OUT_OF_RANGE);
    }

    err=gsl_blas_dgemv(CblasNoTrans, 1., v, pp, 0., interior);

    return err;

  }

  int check_constraints(gsl_matrix *v,
		  gsl_vector *c,
		  gsl_vector *x,
		  int *ind) {
    int n=0;
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
      if (dot>gsl_vector_get(c, i)) {
        ind[n]=i;
	n++;
      }
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
    double v_cv;	//constraint coeff. at cindex, vindex
    double c_c;		//constraint limit at cindex
    //we eliminate the first variable with a non-zero coefficient
    //in the desired constraint:
    for (int i=0; i<a->size2; i++) {
      if (gsl_matrix_get(v, i)!=0) {
        vindex=i;
	break;
      }
    }
    v_cv=gsl_matrix_get(v, cindex, vindex);
    c_c=gsl_vector_get(c, cindex);
    //apply constraint to the constraints themselves:
    for (int i=0; i<cindex; i++) {
      double c_el=gsl_vector_get(c, i);
      double v_iv=gsl_matrix_get(v, i, vindex);
      gsl_vector_set(c2, i, gsl_matrix_get(c, i)-c_c);
      for (int j=0; j<vindex; j++) {
        double v_el=gsl_matrix_get(v, i, j);
	double v_cj=gsl_matrix_get(v, cindex, j);
	gsl_matrix_set(v2, i, j, v_el-v_cj/v_iv);
      }
      for (int j=vindex+1; j<v->size2; j++) {
        double v_el=gsl_matrix_get(v, i, j);
	double v_cj=gsl_matrix_get(v, cindex, j);
	gsl_matrix_set(v2, i, j-1, v_el-v_cj/v_iv);
      }
    }
    for (int i=cindex+1; i<v->size1; i++) {
      double c_el=gsl_vector_get(c, i);
      double v_iv=gsl_matrix_get(v, i, vindex);
      gsl_vector_set(c2, i-1, gsl_matrix_get(c, i)-c_c);
      for (int j=0; j<vindex; j++) {
        double v_el=gsl_matrix_get(v, i, j);
	double v_cj=gsl_matrix_get(v, cindex, j);
	gsl_matrix_set(v2, i-1, j, v_el-v_cj/v_iv);
      }
      for (int j=vindex+1; j<v->size2; j++) {
        double v_el=gsl_matrix_get(v, i, j);
	double v_cj=gsl_matrix_get(v, cindex, j);
	gsl_matrix_set(v2, i-1, j-1, v_el-v_cj/v_iv);
      }
    }
    //apply constraints to the matrix and solution vector:
    for (int i=0; i<cindex; i++) {
      double b_el=gsl_vector_get(b, i);
      double v_iv=gsl_matrix_get(v, i, vindex);
      gsl_vector_set(b2, i, gsl_matrix_get(b, i)-c_c);
      for (int j=0; j<vindex; j++) {
        double a_el=gsl_matrix_get(v, i, j);
	double v_cj=gsl_matrix_get(v, cindex, j);
	gsl_matrix_set(v2, i, j, a_el-v_cj/v_iv);
      }
      for (int j=vindex+1; j<v->size2; j++) {
        double a_el=gsl_matrix_get(a, i, j);
	double v_cj=gsl_matrix_get(v, cindex, j);
	gsl_matrix_set(a2, i, j-1, a_el-v_cj/v_iv);
      }
    }
  }

  int constrained(gsl_matrix *a,	//matrix to solve
		gsl_vector *b,		//solution vector
		gsl_matrix *v,		//constraint normals
		gsl_vector *c,		//constraint thresholds
    		gsl_vector *interior,	//interior point
		gsl_vector *x){		//result

    gsl_vector *xtrial;		//unconstrained solution
    gsl_vector *interior;	//interior point
    int bind[c->size];		//which constraints are broken
    int nb;			//number of broken constraints

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
    solver(a, b, xtrial);
    nb=check_constraints(v, c, xtrail, bind);

