#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

#include "error_codes.h"

#define FUCK(X) X-1
#define TM 1

namespace libpetey {

  int SV_solve_transpose(gsl_matrix *vt,
		  gsl_matrix *u,
		  gsl_vector *s,
		  gsl_vector *b,
		  gsl_vector *x) {
    if (vt->size1 != s->size ||
		    vt->size2 != b->size ||
		    s->size != b->size ||
		    u->size1 != x->size ||
		    u->size2 != s->size) {
      fprintf(stderr, "SV_solve_transpose: dimension mismatch\n");
      throw DIMENSION_MISMATCH;
    }
    for (int i=0; i<u->size2; i++) {
      double xi=0;
      gsl_vector_set(x, i, 0);
      for (int j=0; j<s->size; j++) {
        double vtb=0;
        for (int k=0; k<vt->size2; k++) {
          vtb+=gsl_matrix_get(vt, j, k)*gsl_vector_get(b, k);
	}
	xi+=gsl_matrix_get(u, i, j)*vtb/gsl_vector_get(s, j);
      }
      gsl_vector_set(x, i, xi);
    }
    return 0;
  }

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
      throw DIMENSION_MISMATCH;
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
      SV_solve_transpose(vt, u, s, b, x);
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
		    interior->size!= v->size2) {
      fprintf(stderr, "find_interior: dimension mismatch\n");
      throw DIMENSION_MISMATCH;
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
        throw PARAMETER_OUT_OF_RANGE;
      }
      for (int i=0; i<v->size2; i++) {
        double pp_el=gsl_vector_get(pp, i);
        gsl_vector_set(pp, pp_el+dot/pp_el/a->size2);
      }
      gsl_vector_free(xp);
    } else {
      fprintf(stderr, "get_interior: case for less than n+1 contraints not implemented yet, sorry.\n");
      throw PARAMETER_OUT_OF_RANGE;
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
      throw DIMENSION_MISMATCH;
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

  //***NOTE***: matrix coefficients do not transform the same as vectors
  //themselves!!
  //apply constraint to a single vector:
  //if index==-1 then we use the first non-zero value in the constraint
  //normal, the value for this index is returned
  int constrain_row(gsl_vector *row,		//row to transform
		  gsl_vector *constraint, 	//constraint normal
		  int index,			//index of variable to eliminate
		  gsl_vector *row2) {		//transformed row
    assert(row->size == constraint->size ||
		    row->size == row2->size+1);
    if (index == -1) {
      for (int i=0; i<constraint->size; i++) {
        if (gsl_matrix_get(constraint, i)!=0) {
          index=i;
          break;
        }
      }
    }
    assert(index<x->size && index>=0);

    double v_i=gsl_vector_get(constraint, index);
    for (int j=0; j<index; j++) {
      double x_el=gsl_vector_get(row, j);
      double v_j=gsl_vector_get(constraint, j);
      gsl_vector_set(row2, j, x_el-v_cj/v_i);
    }
    for (int j=index+1; j<row->size; j++) {
      double x_el=gsl_vector_get(row, j);
      double v_j=gsl_vector_get(constraint, j);
      gsl_vector_set(row2, j-1, x_el-v_j/v_i);
    }

    return index;
  }
    
  //returns the index of the eliminated variable:
  int apply_constraint(gsl_matrix *a, 
		gsl_vector *b,		//solution vector
		gsl_matrix *v,		//constraint normals
		gsl_vector *c,		//constraint thresholds
		int cindex,		//index of constraint
		int vindex,		//index of variable to eliminate
		gsl_matrix *a2,
		gsl_matrix *b2,
		gsl_matrix *v2,
		gsl_matrix *c2) {
    double c_c=gsl_vector_get(c, cindex);	//constraint limit at cindex
    //pull out the desired constraint:
    gsl_vector_view applied_constraint=gsl_matrix_row(v, cindex);
    double v_cv;		//constraint coeff.

    //apply constraint to the constraints themselves:
    for (int i=0; i<cindex; i++) {
      //pull out row views for original matrices:
      gsl_vector_view vrow=gsl_matrix_row(v, i);
      //pull out row views for transformed matrices:
      gsl_vector_view v2row=gsl_matrix_row(v2, i);

      //then pass them to apply_constraint for vector only:
      vindex=constrain_row(&vrow.vector, &applied_constraint.vector, 
		      vindex, &v2row.vector);
      v_cv=gsl_matrix_get(v, cindex, vindex);
      gsl_vector_set(c2, i, gsl_matrix_get(c, i)-c_c/v_cv);

    }
    for (int i=cindex+1; i<v->size1; i++) {
      //pull out row views for original matrices:
      gsl_vector_view vrow=gsl_matrix_row(v, i);
      //pull out row views for transformed matrices:
      gsl_vector_view v2row=gsl_matrix_row(v2, i-1);

      //then pass them to apply_constraint for vector only:
      constrain_row(&vrow.vector, &applied_constraint.vector, 
		      vindex, &v2row.vector);
      gsl_vector_set(c2, i-1, gsl_matrix_get(c, i)-c_c/v_cv);
    }
    //apply constraints to the matrix and solution vector:
    for (int i=0; i<a->size1; i++) {
      //set-up:
      gsl_vector_view arow=gsl_matrix_row(a, i);
      gsl_vector_view a2row=gsl_matrix_row(a2, i);
      //actual work:
      constrain_row(&arow.vector, &applied_constraint.vector, 
		      vindex, &a2row.vector);
      gsl_vector_set(b2, i, gsl_matrix_get(b, i)-c_c/v_cv);
    }

    return vindex;
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
      throw DIMENSION_MISMATCH;
    }

    //first, we find the unconstrained solution:
    solver(a, b, xtrial);

    //see how many constraints it violates:
    nb=check_constraints(v, c, xtrail, bind);
    //find distance between the interior point and each of the violated
    //constraints along the line to the unconstrained solution:
    //(in other words, we always keep both a solution and a point that satisfies
    //all the constraints, so in a sense "bracketing" the constrained solution)
    if (nb>0) {
      gsl_vector *int2;		//new interior point
      gsl_matrix *a2;		//reduced problem
      gsl_vector *b2;		//reduced, transformed solution vector
      gsl_vector *v2;		//reduced, transformed constraint normals
      gsl_vector *c2;		//reduced, transformed constraint thresholds

      double dir[x->size];		//direction vector
      double s[nb];			//line parameter for each broken constraint
      double smin=0;		//min. line parameter
      double bimin=-1;		//index for constraint closest to interior point in the direction of the unconstrained solution
      for (int i=0; i<nb; i++) {
        double vdotx0=0;	//constraint dotted with unconstr. soln.
	double vdotint=0;	//constraint normal dotted with interior pt.
	for (int j=0; j<x->size; j++) {
          double vind=gsl_matrix_get(v, bind[i], j);
          vdotx0+=vel*gsl_vector_get(xtrial, j);
	  vdotint+=vel*gsl_vector_get(interior, j);
	}
	s[i]=(gsl_vector_get(c, bind[i])-vdotint-vdotx0)/vdotint;
	if (s[i]<smin) {
          smin=s[i];
	  bimin=i;
	}
	if (bimin==-1) {
          fprintf(stderr, "constrained: something went wrong\n");
	  throw INTERNAL_ERROR;
	}

	//allocate new variables for solving reduced problem:
	a2=gsl_matrix_alloc(a->size1, a->size2-1);
	b2=gsl_matrix_alloc(a->size2);
	v2=gsl_matrix_alloc(v->size1-1, v->size2-1);
	c2=gsl_matrix_alloc(v->size1-1);

	//apply the constraint:
	int vind=-1;
	vind=apply_constraint(a, b, v, c, bind[bimin], vind, a2, b2, v2, c2);

	//find the new interior point which is the location along the line
	//between the old interior point and the unconstrained solution
	//at the eliminated constraint:
	gsl_vector *interior2=gsl_vector_alloc(interior->size-1);
	for (int i=0; i<vind; i++) {
          gsl_vector_set(interior2, i, (1-smin)*gsl_vector_get(interior, i) +
			  smin*gsl_vector_get(xtrial, i));
	}
	for (int i=vind+1; i<interior->size; i++) {
          gsl_vector_set(interior2, i+1, (1-smin)*gsl_vector_get(interior, i) +
			  smin*gsl_vector_get(xtrial, i));
	}

	//same problem, now slightly reduced in scale:
	gsl_vector *x2=gsl_vector_alloc(x->size-1);
	constrained(a2, b2, v2, c2, interior2t, x2);

	//reconstitute eliminated variable:
	double x_v=gsl_vector_get(c, bind[bimin]);
	for (int i=0; i<vind; i++) {
          x2_i=gsl_vector_get(x2, i);
          gsl_vector_set(x, i, x2_i);
          x_v-=x2_i*gsl_matrix_get(v, bind[bimin], i);
	}
	for (int i=vind+1; i<interior->size; i++) {
          x2_i=gsl_vector_get(x2, i-1);
          gsl_vector_set(x, i, x2_i);
          x_v-=x2_i*gsl_matrix_get(v, bind[bimin], i);
	}
	gsl_vector_set(x, vind, x_v/gsl_matrix_get(v, bind[bimin], vind);

	//delete a shit-load of variables:
	gsl_matrix_free(a2);
	gsl_vector_free(b2);
	gsl_matrix_free(v2);
	gsl_vector_free(c2);
	gsl_vector_free(interior2);
	gsl_vector_free(x2);
      }

      //if all constraints are satisfied, then we're done! yay!
      
      return FUCK(TM);

    }

