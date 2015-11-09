#include <math.h>
#include <assert.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

#include "gsl_util.h"
#include "error_codes.h"
#include "randomize.h"
#include "constrained.h"

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
    gsl_matrix *u;
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
      u=gsl_matrix_alloc(m, n);
      gsl_matrix_transpose_memcpy(u, a);
    } else {
      m=a->size1;
      n=a->size2;
      u=gsl_matrix_alloc(a->size1, a->size2);
      gsl_matrix_memcpy(u, a);
    }

    vt=gsl_matrix_alloc(n, n);
    s=gsl_vector_alloc(n);
    work=gsl_vector_alloc(n);

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

  //given a set of constraints, finds and interior point:
  int find_interior(
		gsl_matrix *v,		//constraint normals
		gsl_vector *c,		//constraint thresholds
    		gsl_vector *pp,		//interior point
		double offset,		//offset (if applicable)
	  	int (* solver1) (gsl_matrix *, gsl_vector *, gsl_vector *)) {
    int err=0;
    int nb;
    int bind[pp->size];

    if (v->size1 != c->size ||
		    pp->size!= v->size2) {
      fprintf(stderr, "find_interior: dimension mismatch\n");
      throw DIMENSION_MISMATCH;
    }

    //d=v_k*V_k^-1*z_0 - c_k
    //k=n
    if (c->size==v->size2+1) {
      gsl_vector_view v_k = gsl_matrix_row(v, v->size2);
      gsl_vector *xp=gsl_vector_alloc(v->size2);
      double dot;		//dot product--must be less than c_n
      double sum_v_k=0;		//sum of elements of nth constraint
      double c_k=gsl_vector_get(c, v->size2);
      c->size--;		//z_0i=c_i; k ommitted
      v->size1--;
      err=(*solver1)(v, c, xp);
      if (err!=0) {
        fprintf(stderr, "find_interior: solver returned error, %d\n", err);
        throw err;
      }
      c->size++;
      v->size1++;
      //dot=cblas_ddot(c->size, c->data, c->stride, xp->data, xp->stride);
      err=gsl_blas_ddot(&v_k.vector, xp, &dot);
      if (err!=0) {
        fprintf(stderr, "find_interior: gsl_blas_ddot returned error, %d\n", err);
        throw err;
      }
      if (dot < c_k) {
        fprintf(stderr, "get_interior: contraints not consistent\n");
	gsl_vector_free(xp);
        throw PARAMETER_OUT_OF_RANGE;
      }
      //p=V_k^-1*z_0 + f*d*1/v_k*1
      for (int i=0; i<v->size2; i++) sum_v_k+=gsl_vector_get(&v_k.vector, i);
      for (int i=0; i<v->size2; i++) {
        double xp_el=gsl_vector_get(xp, i);
	gsl_vector_set(pp, i, xp_el+offset*(c_k-dot)/sum_v_k);
      }
      gsl_vector_free(xp);
    } else {
      fprintf(stderr, "get_interior: case for less than n+1 contraints not implemented yet, sorry.\n");
      throw PARAMETER_OUT_OF_RANGE;
    }

    nb=check_constraints(v, c, pp, bind);
    assert(nb==0);

    //err=gsl_blas_dgemv(CblasNoTrans, 1., v, pp, 0., interior);

    return err;

  }

  void normalize_constraints(gsl_matrix *v,
		  gsl_vector *c) {
    assert(v->size1==c->size);
    
    for (int i=0; i<v->size1; i++) {
      double d=0;
      for (int j=0; j<v->size2; j++) {
        double v_j=gsl_matrix_get(v, i, j);
	d+=v_j*v_j;
      }
      d=sqrt(d);
      for (int j=0; j<v->size2; j++) {
	gsl_matrix_set(v, i, j, gsl_matrix_get(v, i, j)/d);
      }
      gsl_vector_set(c, i, gsl_vector_get(c, i)/d);
    }
  }

  //tests each of the constraint and returns the indices of each of the
  //violated constraints:
  //returns the number of violated constraints
  int check_constraints(gsl_matrix *v,		//matrix of constraint normals
		  gsl_vector *c,		//vector of constraint thresholds
		  gsl_vector *x,		//point to test
		  int *ind) {			//list of violated constraints
    int n=0;
    double eps=1e-15;

    if (v->size1 != c->size || v->size2 != x->size) {
      fprintf(stderr, "check constraints: dimension mismatch\n");
      throw DIMENSION_MISMATCH;
    }
    for (int i=0; i<v->size1; i++) {
      double *vrow;
      double dot;
      double c_i=gsl_vector_get(c, i);
      vrow=gsl_matrix_ptr(v, i, 0);
      dot=cblas_ddot(v->size2, vrow, 1, x->data, x->stride);
      //we have to fudge this a little bit otherwise it screws up:
      //(although I wonder if the errors might not compound after a while??)
      if (dot<c_i) {
      //if (c_i-dot > eps) {
        ind[n]=i;
	n++;
      }
    }
    return n;
  }

  int test_interior(int n) {		//dimension of problem
    gsl_vector *x0=gsl_vector_alloc(n);	//random point in interior of constraint hyper-pyramid
    gsl_matrix *v=gsl_matrix_alloc(n+1, n);	//set of random constraint normals
    gsl_vector *c=gsl_vector_alloc(n+1);	//constraint constants
    gsl_vector *p=gsl_vector_alloc(n);
    int err;
    int ind[n];

    for (int i=0; i<n; i++) gsl_vector_set(x0, i, rang());

    do {
      //create a set of random constraint normals:
      for (int i=0; i<n+1; i++) {
        double c_i=0;
        for (int j=0; j<n; j++) {
          double val=rang();
	  double v_ij=gsl_vector_get(x0, j)-val;
          gsl_matrix_set(v, i, j, v_ij);
	  c_i+=v_ij*val;
	}
      }
      //constraints may not enclose a volume but find_interior should catch
      //this:
      err=find_interior(v, c, p, 0.5);
    } while (err!=0);

    err=check_constraints(v, c, p, ind);

    if (err!=0) {
      fprintf(stderr, "test_interior: unit test for find_interior and check_constraints failed!\n");
      return -1;
    }
    return 0;
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
        if (gsl_vector_get(constraint, i)!=0) {
          index=i;
          break;
        }
      }
    }
    assert(index < row->size && index >= 0);

    double v_i=gsl_vector_get(constraint, index);
    double x_i=gsl_vector_get(row, index);
    for (int j=0; j<index; j++) {
      double x_el=gsl_vector_get(row, j);
      double v_j=gsl_vector_get(constraint, j);
      gsl_vector_set(row2, j, x_el-x_i*v_j/v_i);
    }
    for (int j=index+1; j<row->size; j++) {
      double x_el=gsl_vector_get(row, j);
      double v_j=gsl_vector_get(constraint, j);
      gsl_vector_set(row2, j-1, x_el-x_i*v_j/v_i);
    }

    return index;
  }
    
  //returns the index of the eliminated variable:
  int apply_constraint( 
		gsl_matrix *v,		//constraint normals
		gsl_vector *c,		//constraint thresholds
		int cindex,		//index of constraint
		int vindex,		//index of variable to eliminate
		gsl_matrix *v2,
		gsl_vector *c2) {
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
      gsl_vector_set(c2, i, gsl_vector_get(c, i)-
		      gsl_matrix_get(v, i, vindex)*c_c/v_cv);

    }
    for (int i=cindex+1; i<v->size1; i++) {
      //pull out row views for original matrices:
      gsl_vector_view vrow=gsl_matrix_row(v, i);
      //pull out row views for transformed matrices:
      gsl_vector_view v2row=gsl_matrix_row(v2, i-1);

      //then pass them to apply_constraint for vector only:
      vindex=constrain_row(&vrow.vector, &applied_constraint.vector, 
		      vindex, &v2row.vector);
      v_cv=gsl_matrix_get(v, cindex, vindex);
      gsl_vector_set(c2, i-1, gsl_vector_get(c, i)-
		      gsl_matrix_get(v, i, vindex)*c_c/v_cv);
    }

    return vindex;
  }

  //returns the index of the eliminated variable:
  int apply_constraint(gsl_matrix *a, 
		gsl_vector *b,		//solution vector
		gsl_matrix *v,		//constraint normals
		gsl_vector *c,		//constraint thresholds
		int cindex,		//index of constraint
		int vindex,		//index of variable to eliminate
		gsl_matrix *a2,
		gsl_vector *b2,
		gsl_matrix *v2,
		gsl_vector *c2) {
    double c_c=gsl_vector_get(c, cindex);	//constraint limit at cindex
    //pull out the desired constraint:
    gsl_vector_view applied_constraint=gsl_matrix_row(v, cindex);
    double v_cv;		//constraint coeff.

    assert(a->size2==v->size2 && a->size1==b->size && v->size1==c->size
		    && a2->size2==a->size2-1 && a2->size1==a->size1 
		    && v2->size1==v->size1-1 && v2->size2==a2->size2
		    && c2->size==v2->size1 && b2->size==a2->size1);

    //apply constraint to the constraints themselves:
    for (int i=0; i<cindex; i++) {
      //pull out row views for original matrices:
      gsl_vector_view vrow=gsl_matrix_row(v, i);
      //pull out row views for transformed matrices:
      gsl_vector_view v2row=gsl_matrix_row(v2, i);

      //then pass them to apply_constraint for vector only:
      vindex=constrain_row(&vrow.vector, &applied_constraint.vector, 
		      vindex, &v2row.vector);
      assert(vindex >=0 && vindex < a->size2);
      v_cv=gsl_matrix_get(v, cindex, vindex);
      gsl_vector_set(c2, i, gsl_vector_get(c, i)-
		      gsl_matrix_get(v, i, vindex)*c_c/v_cv);

    }
    for (int i=cindex+1; i<v->size1; i++) {
      //pull out row views for original matrices:
      gsl_vector_view vrow=gsl_matrix_row(v, i);
      //pull out row views for transformed matrices:
      gsl_vector_view v2row=gsl_matrix_row(v2, i-1);

      //then pass them to apply_constraint for vector only:
      vindex=constrain_row(&vrow.vector, &applied_constraint.vector, 
		      vindex, &v2row.vector);
      v_cv=gsl_matrix_get(v, cindex, vindex);
      gsl_vector_set(c2, i-1, gsl_vector_get(c, i)-
		      gsl_matrix_get(v, i, vindex)*c_c/v_cv);
    }
    //apply constraints to the matrix and solution vector:
    for (int i=0; i<a->size1; i++) {
      //set-up:
      gsl_vector_view arow=gsl_matrix_row(a, i);
      gsl_vector_view a2row=gsl_matrix_row(a2, i);
      //actual work:
      constrain_row(&arow.vector, &applied_constraint.vector, 
		      vindex, &a2row.vector);
      gsl_vector_set(b2, i, gsl_vector_get(b, i)-
		      gsl_matrix_get(a, i, vindex)*c_c/v_cv);
    }
    //printf("apply_constraint: vindex=%d\n", vindex);

    return vindex;
  }

  int constrained(gsl_matrix *v,	//transformed constraint normals
		  gsl_vector *c,	//transformed constraint thresholds
		  gsl_vector *p,	//interior point
		  gsl_vector *x) {
    int err=0;
    int nb;		//number of broken constraints
    int bind[c->size];		//list of broken constraints
    int n=x->size;

    assert(p->size==x->size);
    assert(x->size==v->size2);
    assert(v->size2==v->size1-1);
    assert(c->size==v->size1);

    //check for violated constraints:
    nb=check_constraints(v, c, x, bind);
    //if a constraint has been violated, find the least distance to the
    //constraint hyperplane:
    if (nb > 0) {
      double smin=1;		//distance to nearest constraint
      int minind=-1;		//index of nearest constraint
      int vind=-1;		//index of next eliminated variable
      gsl_vector *pn=gsl_vector_alloc(n);
      gsl_vector *xn=gsl_vector_alloc(n);
      gsl_vector *pp=gsl_vector_alloc(n-1);
      gsl_vector *xp=gsl_vector_alloc(n-1);
      gsl_matrix *vp=gsl_matrix_alloc(n, n-1);
      gsl_vector *cp=gsl_vector_alloc(n);

      for (int i=0; i<nb; i++) {
        double s;
	double vdotp, vdotx;
        gsl_vector_view v_ind=gsl_matrix_row(v, bind[i]);

	gsl_blas_ddot(&v_ind.vector, p, &vdotp);
	gsl_blas_ddot(&v_ind.vector, x, &vdotx);

	s=(gsl_vector_get(c, bind[i])-vdotp)/(vdotx-vdotp);
	assert(s>=0 && s<=1);
	if (s<smin) {
          smin=s;
	  minind=bind[i];
	}
      }
      assert(minind>=0);

      //transform the other constraints according to the closest violated
      //constraint:
      vind=apply_constraint(v, c, minind, vind, vp, cp);

      //find the location of the both the new interior point and the solution
      //on the constraint hyperplane as projected onto the new variables:
      gsl_vector_memcpy(pn, p);
      gsl_vector_scale(pn, -1);
      gsl_vector_add(pn, x);
      gsl_vector_scale(pn, smin);
      gsl_vector_add(pn, p);

      double vidotx;
      double vdotv;
      gsl_vector_view v_i=gsl_matrix_row(v, minind);

      gsl_blas_ddot(&v_i.vector, x, &vidotx);
      gsl_blas_ddot(&v_i.vector, &v_i.vector, &vdotv);

      gsl_vector_memcpy(xn, &v_i.vector);
      gsl_vector_scale(xn, (gsl_vector_get(c, minind)-vidotx)/vdotv);
      gsl_vector_add(xn, x);

      //exclude chosen variables:
      for (int i=0; i<vind; i++) {
        gsl_vector_set(xp, i, gsl_vector_get(xn, i));
        gsl_vector_set(pp, i, gsl_vector_get(pn, i));
      }
      for (int i=vind; i<n-1; i++) {
        gsl_vector_set(xp, i, gsl_vector_get(xn, i+1));
        gsl_vector_set(pp, i, gsl_vector_get(pn, i+1));
      }

      //restart the solution:
      constrained(vp, cp, pp, xp);

      //reconstitute the missing variable before exiting:
      double xvind=gsl_vector_get(c, minind);

      for (int i=0; i<vind; i++) {
	double xp_i=gsl_vector_get(xp, i);
	gsl_vector_set(x, i, xp_i);
        xvind-=gsl_matrix_get(v, minind, i)*xp_i;
      }
      for (int i=vind; i<n-1; i++) {
	double xp_i=gsl_vector_get(xp, i);
	gsl_vector_set(x, i+1, xp_i);
        xvind-=gsl_matrix_get(v, minind, i+1)*xp_i;
      }
      gsl_vector_set(x, vind, xvind/gsl_matrix_get(v, minind, vind));

      //clean up:
      gsl_vector_free(xn);
      gsl_vector_free(pn);
      gsl_vector_free(xp);
      gsl_vector_free(pp);
      gsl_matrix_free(vp);
      gsl_vector_free(cp);
    }

    return err;

  }


  int constrained(gsl_matrix *a,	//matrix to solve
		gsl_vector *b,		//solution vector
		gsl_matrix *v,		//constraint normals
		gsl_vector *c,		//constraint thresholds
    		gsl_vector *interior,	//interior point
		gsl_vector *x,		//result
	  	int (* solver1) (gsl_matrix *, gsl_vector *, gsl_vector *)) {

    gsl_vector *xtrial;		//unconstrained solution
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
    //xtrial=gsl_vector_alloc(x->size);
    (*solver1)(a, b, x);
/*
    printf("constrained: solve:\n");
    printf("a=\n");
    print_gsl_matrix(stdout, a);
    printf("b=");
    gsl_vector_fprintf(stdout, b, "%g ");
    printf("\n");
    printf("subject to:\n");
    printf("v=\n");
    print_gsl_matrix(stdout, v);
    printf("c=");
    gsl_vector_fprintf(stdout, c, "%g ");
    printf("trial solution=\n");
    gsl_vector_fprintf(stdout, x, "%g ");
    printf("interior point=\n");
    gsl_vector_fprintf(stdout, interior, "%g ");
*/
    //see how many constraints it violates:
    nb=check_constraints(v, c, x, bind);
    //find distance between the interior point and each of the violated
    //constraints along the line to the unconstrained solution:
    //(in other words, we always keep both a solution and a point that satisfies
    //all the constraints, so in a sense "bracketing" the constrained solution)
    if (nb>0) {
      //if there is only one variable left, we solve based on the violated
      //constraint: (vx=c)
      /*
      printf("constraints broken: ");
      for (int i=0; i<nb; i++) printf("%d ", bind[i]);
      printf("\n");
      */
      if (x->size==1) {
        assert(nb==1);	//if problem has been specified properly, only
	  			//one constraint can be violated
        gsl_vector_set(x, 0, gsl_vector_get(c, bind[0])/gsl_matrix_get(v, bind[0], 0));
        return 0;
      }

      gsl_vector *int2;		//new interior point
      gsl_matrix *a2;		//reduced problem
      gsl_vector *b2;		//reduced, transformed solution vector
      gsl_matrix *v2;		//reduced, transformed constraint normals
      gsl_vector *c2;		//reduced, transformed constraint thresholds

      double dir[x->size];		//direction vector
      double s[nb];			//line parameter for each broken constraint
      double smin=1;		//min. line parameter
      int bimin=0;		//index for constraint closest to interior point in the direction of the unconstrained solution
      for (int i=0; i<nb; i++) {
        double vdotx0=0;	//constraint dotted with unconstr. soln.
	double vdotint=0;	//constraint normal dotted with interior pt.
	for (int j=0; j<x->size; j++) {
          double vel=gsl_matrix_get(v, bind[i], j);
          vdotx0+=vel*gsl_vector_get(x, j);
	  vdotint+=vel*gsl_vector_get(interior, j);
	}
//***this will not work in all cases!
//must be redone!
	s[i]=(gsl_vector_get(c, bind[i])-vdotint)/(vdotx0-vdotint);
	//this is the source of our problem; solution: just get rid of it...
	//assert(s[i]>=0 && s[i]<=1);
	if (s[i]<smin) {
          smin=s[i];
	  bimin=i;
	}
      }
      if (bimin==-1) {
        fprintf(stderr, "constrained: something went wrong\n");
        throw INTERNAL_ERROR;
      }
      /*
      printf("distances to each constraint: ");
      for (int i=0; i<nb; i++) printf("%g ", s[i]);
      printf("\n");
      */
      //allocate new variables for solving reduced problem:
      a2=gsl_matrix_alloc(a->size1, a->size2-1);
      b2=gsl_vector_alloc(a->size1);
      v2=gsl_matrix_alloc(v->size1-1, v->size2-1);
      c2=gsl_vector_alloc(v->size1-1);

      //apply the constraint:
      int vind=-1;
      vind=apply_constraint(a, b, v, c, bind[bimin], vind, a2, b2, v2, c2);

      //find the new interior point which is the location along the line
      //between the old interior point and the unconstrained solution
      //at the eliminated constraint:
      int2=gsl_vector_alloc(interior->size-1);
      for (int i=0; i<vind; i++) {
        gsl_vector_set(int2, i, (1-smin)*gsl_vector_get(interior, i) +
			  smin*gsl_vector_get(x, i));
      }
      for (int i=vind+1; i<interior->size; i++) {
        gsl_vector_set(int2, i-1, (1-smin)*gsl_vector_get(interior, i) +
			  smin*gsl_vector_get(x, i));
      }

      //same problem, now slightly reduced in scale:
      gsl_vector *x2=gsl_vector_alloc(x->size-1);
      constrained(a2, b2, v2, c2, int2, x2);

      //reconstitute eliminated variable:
      double x_v=gsl_vector_get(c, bind[bimin]);
      for (int i=0; i<vind; i++) {
        double x2_i=gsl_vector_get(x2, i);
        gsl_vector_set(x, i, x2_i);
        x_v-=x2_i*gsl_matrix_get(v, bind[bimin], i);
      }
      for (int i=vind+1; i<interior->size; i++) {
        double x2_i=gsl_vector_get(x2, i-1);
        gsl_vector_set(x, i, x2_i);
        x_v-=x2_i*gsl_matrix_get(v, bind[bimin], i);
      }
      gsl_vector_set(x, vind, x_v/gsl_matrix_get(v, bind[bimin], vind));

      //delete a shit-load of variables:
      gsl_matrix_free(a2);
      gsl_vector_free(b2);
      gsl_matrix_free(v2);
      gsl_vector_free(c2);
      gsl_vector_free(int2);
      gsl_vector_free(x2);
    }

    //if all constraints are satisfied, then we're done! yay!
      
    return FUCK(TM);

  }

  int constrained(gsl_matrix *a,	//matrix to solve
		gsl_vector *b,		//solution vector
		gsl_matrix *v,		//constraint normals
		gsl_vector *c,		//constraint thresholds
		gsl_vector *x,		//result
	  	int (* solver1) (gsl_matrix *, gsl_vector *, gsl_vector *)) {

    gsl_vector *p;		//interior point
    int err=0;

    if (a->size1 != b->size ||
		    a->size2 != x->size ||
		    v->size2 != a->size2 ||
		    v->size1 > a->size2+1 ||
		    c->size != v->size1) {
      fprintf(stderr, "constrained: dimension mismatch\n");
      throw DIMENSION_MISMATCH;
    }

    gsl_error_handler_t *old_handler=gsl_set_error_handler(&gsl_throw_handler);

    p=gsl_vector_alloc(x->size);
    find_interior(v, c, p, 0.5);

    //if the matrix is square, we can use the more efficient method:
    //if (a->size1 == a->size2) {
    //there is actually no advantage to doing it this way:
    if (0) {
      gsl_matrix *vt=gsl_matrix_alloc(v->size1, v->size2);	//transformed constraint normals
      gsl_vector *ct=gsl_vector_alloc(c->size);			//transformed constraint thresholds
      gsl_vector *xt=gsl_vector_alloc(b->size);
      gsl_vector_view v_i;
      gsl_vector_view vt_i;
      double ct_i;
      //procedure is pretty basic: transform the constraints, then pass 
      //everything to the procedure that minimizes just |x|, then transform
      //the results back...
      gsl_matrix_transpose(a);
      for (int i=0; i<v->size1; i++) {
        v_i=gsl_matrix_row(v, i);
	vt_i=gsl_matrix_row(vt, i);
	(*solver1)(a, &v_i.vector, &vt_i.vector);
        gsl_blas_ddot(&vt_i.vector, b, &ct_i);
	gsl_vector_set(ct, i, gsl_vector_get(c, i)+ct_i);
      }

      p=gsl_vector_alloc(b->size);
      find_interior(vt, ct, p, 0.5);

      constrained(vt, ct, p, xt);

      //transform result back:
      gsl_matrix_transpose(a);
      gsl_vector_sub(xt, b);
      (*solver1)(a, xt, x);

      gsl_matrix_free(vt);
      gsl_vector_free(ct);
      gsl_vector_free(xt);
    } else {
      p=gsl_vector_alloc(x->size);
      find_interior(v, c, p, 0.5);
      constrained(a, b, v, c, p, x);

      //we've solved it using the flawed algorithm, lets see if the criteria
      //for the algorithm being always valid are satisfied:
      //calculate a^t*a:
      gsl_matrix *ata=gsl_matrix_alloc(a->size2, a->size2);
      gsl_matrix *at=gsl_matrix_alloc(a->size2, a->size1);
      gsl_vector *t1=gsl_vector_alloc(a->size1);
      gsl_vector *t2=gsl_vector_alloc(a->size1);
      gsl_matrix_transpose_memcpy(at, a);
      gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1., a, a, 0., ata);
      print_gsl_matrix(stdout, a);
      print_gsl_matrix(stdout, ata);
      double vpmag[a->size2];
      //calculate v_i (A^T A)^{-1} v_j/sqrt(v_i (A^T A) v_i)/sqrt(v_j (A^T A)^{-1} v_j) for each i and j:
      for (int i=0; i<v->size1; i++) {
	gsl_vector_view v_i=gsl_matrix_row(v, i);
        //(*solver1) (ata, &v_i.vector, t1);
        (*solver1) (at, &v_i.vector, t1);
	//gsl_blas_ddot(t1, &v_i.vector, vpmag+i);
	gsl_blas_ddot(t1, t1, vpmag+i);
	vpmag[i]=sqrt(vpmag[i]);
	for (int j=0; j<i; j++) {
	  gsl_vector_view v_j=gsl_matrix_row(v, j);
          double vidotvj;
          (*solver1) (at, &v_j.vector, t2);
	  //gsl_blas_ddot(t1, &v_j.vector, &vidotvj);
	  gsl_blas_ddot(t1, t2, &vidotvj);
	  //printf("%12.5g ", vidotvj/sqrt(vpmag[i]*vpmag[j]));
	  printf("%12.5g ", vidotvj);
	  printf("%12.5g ", sqrt(vpmag[i]*vpmag[j]));
	}
	printf("\n");
      }
      gsl_vector_free(t1);
      gsl_vector_free(t2);
      gsl_matrix_free(ata);
      gsl_matrix_free(at);
    }

    gsl_vector_free(p);
    gsl_set_error_handler(old_handler);

    return err;
  }

} //end namespace libpetey

