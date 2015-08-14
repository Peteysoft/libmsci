#include <gsl/gsl_linalg.h>
#include <gsl/gsl_multifit.h>

#include "full_util.h"
#include "randomize.h"
#include "gsl_util.h"

#include "agf_lib.h"


using namespace libpetey;

namespace libagf {

//strictly for profiler so we don't have to recompile the gsl library:
void gsl_lsq_solver(gsl_matrix *a, gsl_vector *r, gsl_vector *p) {
  //gsl_matrix *u;
  gsl_matrix *vt;
  gsl_vector *work;
  gsl_vector *s;

  int m=a->size1;
  int n=a->size2;

  vt=gsl_matrix_alloc(n, n);
  //u=gsl_matrix_alloc(m, n);
  s=gsl_vector_alloc(n);
  work=gsl_vector_alloc(n);

  //gsl_matrix_memcpy(u, a);
  gsl_linalg_SV_decomp(a, vt, s, work);
  gsl_linalg_SV_solve(a, vt, s, r, p);

  gsl_vector_free(s);
  gsl_vector_free(work);
  gsl_matrix_free(vt);
  //gsl_matrix_free(u);

}

int solve_cond_prob2(gsl_matrix *a,		//decision matrix
		gsl_vector *r,			//original "raw" probabilities
		gsl_vector *p,			//returned cond. prob.
		int *gind,			//current active columns
		int &ng,			//number of active columns
		int iter) {			//iteration
  gsl_matrix *u;
  gsl_matrix *vt;
  gsl_vector *work;
  gsl_vector *s;
  gsl_matrix *a1;
  gsl_vector *p1;
  int m=a->size1;
  int n=a->size2;
  int bind[ng];		//columns to potentially discard
  int n1;

  a1=gsl_matrix_alloc(m, ng);

  //printf("solve_cond_prob: indices:\n");
  //for (int j=0; j<ng; j++) printf(" %d", gind[j]);
  //printf("\n");

  for (int i=0; i<m; i++) {
    for (int j=0; j<ng; j++) {
      gsl_matrix_set(a1, i, j, gsl_matrix_get(a, i, gind[j]));
    }
  }

  vt=gsl_matrix_alloc(ng, ng);
  s=gsl_vector_alloc(ng);
  work=gsl_vector_alloc(ng);

  /*
  printf("\n");
  print_gsl_matrix(stdout, a);
  printf("\n");
  */

  //for (int i=0; i<iter; i++) printf("  ");
  //printf("iteration %d:", iter);
  //for (int i=0; i<n; i++) printf(" %g", gsl_vector_get(p, i));

  p1=gsl_vector_alloc(ng);

  gsl_linalg_SV_decomp(a1, vt, s, work);
  gsl_linalg_SV_solve(a1, vt, s, r, p1);

  gsl_vector_free(s);
  gsl_vector_free(work);
  gsl_matrix_free(vt);
  gsl_matrix_free(a1);

  //gsl_lsq_solver(u, r, p);

  //find probabilities less than 0 and pull them out:
  n1=0;
  for (int i=0; i<ng; i++) {
    double p1_i=gsl_vector_get(p1, i);
    //printf(" %d:%g", gind[i], p1_i);
    if (p1_i>=0) {
      gind[n1]=gind[i];
      n1++;
    } else {
      bind[i-n1]=gind[i];
    }
  }
  //printf("\n");

  if (n1==ng) {
    double pmax=0;
    int ind;
    for (int i=0; i<ng; i++) {
      double p1_i=gsl_vector_get(p1, i);
      if (p1_i>pmax) {
        if (p1_i>1) {
          if (pmax>1) {
            fprintf(stderr, "solve_cond_prob2: warning, found two p's greater than 1!\n");
            printf("solve_cond_prob2: warning, found two p's greater than 1!\n");
          }
          pmax=p1_i;
	  ind=i;
	}
      }
    }
    if (pmax>1) {
      for (int i=0; i<ng; i++) gsl_vector_set(p, gind[i], 0);
      gsl_vector_set(p, gind[ind], 1);
    } else {
      for (int i=0; i<ng; i++) gsl_vector_set(p, gind[i], gsl_vector_get(p1, i));
    }
  } else if (n1>0) {
    //if there is more than one value less than one, we have to pull out 
    //each of them in turn and find the solution with the fewest zeroes:
    if (n1<ng-1) {
      int gind1[ng];	//new columns that are still good
      int ng1;
      int gind2[ng];
      int maxng=n1;
      for (int j=0; j<ng; j++) gind2[j]=gind[j];
      //this is terribly inefficient because solutions are repeated:
      for (int i=0; i<ng-n1; i++) {
        for (int j=0; j<n1; j++) gind1[j]=gind[j];
        for (int j=0; j<i; j++) {
	  gind1[n1+j]=bind[j];
	}
        for (int j=i+1; j<ng-n1; j++) {
	  gind1[n1+j-1]=bind[j];
	}
	ng1=ng-1;
        solve_cond_prob2(a, r, p, gind1, ng1, iter+1);
	if (ng1>maxng) {
          maxng=ng1;
	  for (int j=0; j<ng1; j++) gind2[j]=gind1[j];
	}
      }
      if (maxng>ng1) {
        for (int i=0; i<maxng; i++) gind[i]=gind2[i];
	ng=maxng;
      } else {
        for (int i=0; i<n1; i++) gind[i]=gind1[i];
	ng=ng1;
      }
    } else {
      ng=n1;
      solve_cond_prob2(a, r, p, gind, ng, iter+1);
    }
  }

  gsl_vector_free(p1);

  return 0;
}


int solve_cond_prob(gsl_matrix *a,		//decision matrix
		gsl_vector *r,			//original "raw" probabilities
		gsl_vector *p) {		//returned cond. prob.

  //nothing fancy: 1. apply constraint sum_i p_i=1 by removing one variable
  //		and modifying decision matrix and solution vector appropriately
  //		2. if p_i dips below zero, remove ith column
  //		3. if p_i rises above 1, return p_i=1, all others 0
 
  gsl_matrix *a1;
  gsl_vector *r1;

  int m=a->size1-1;
  int n=a->size2;

  int ind;		//use this variable to enforce normalization constraint
  int *gind;		//indices of variable still included
  int ng;		//number of variables still included
  double p_n_1;		//variable left out to enforce norm. constraint

  //apply first constraint:
  gind=new int[n];

  //to avoid any biases produced by using the same variable each time
  ind=ranu()*n;
  //printf("solv_cond_prob: ind=%d\n", ind);

  for (int j=0; j<ind; j++) gind[j]=j;
  for (int j=ind; j<n-1; j++) gind[j]=j+1;
  ng=n-1;

  a1=gsl_matrix_alloc(m, n);
  r1=gsl_vector_alloc(m);

  do {
    for (int i=0; i<n; i++) gsl_vector_set(p, i, 0);

    //apply normalization constraint:
    for (int i=0; i<m; i++) {
      double aind=gsl_matrix_get(a, i, ind);
      gsl_vector_set(r1, i, gsl_vector_get(r, i)-aind);
      for (int j=0; j<ind; j++) {
        gsl_matrix_set(a1, i, j, gsl_matrix_get(a, i, j)-aind);
      }
      for (int j=ind+1; j<n; j++) {
        gsl_matrix_set(a1, i, j, gsl_matrix_get(a, i, j)-aind);
      }
    }

    solve_cond_prob2(a1, r1, p, gind, ng, 1);

    p_n_1=1;
    for (int i=0; i<ng; i++) {
      p_n_1-=gsl_vector_get(p, gind[i]);
    }

    if (p_n_1<0) {
      int ind2=ranu()*ng;
      ind=gind[ind2];
      for (int i=ind2; i<ng-1; i++) gind[i]=gind[i+1];
      ng--;
    } else {
      gsl_vector_set(p, ind, p_n_1);
    }
  } while (p_n_1<0 && ng>0);
      
  delete [] gind;

  gsl_matrix_free(a1);
  gsl_vector_free(r1);

  return 0;
}

//naive method works--produces results within the constraints 
//doesn't improve classification accuracy
//but does improve accuracy of probability estimates!

template <class real, class cls_t>
int derive_probability_map(real **r,		//binary probabilities
		cls_t *truth,			//true classes
		cls_t nmodel,			//number of binary classifiers
		cls_t ncls,			//number of classes
		nel_ta n,			//number of samples
		real (*skill)(cls_t *, cls_t *, nel_ta),
		real **map
		)
{
  cls_t ret[n];		//retrieved classes
  long **sind;		//sorted conditional probabilities
  real **p;		//derived probabilities
  real **dmap;		//minimum change in the mapping to produce a change
  			//in the classes
  real **dsk;		//resultant change in the skill
  real base;		//base skill score

  nel_ta ind;		//index of class that changes
  cls_t newcls;		//class it changes to
  cls_t oldcls;		//the old class

  p=allocate_matrix<real, int32_t>(n, ncls);

  for (nel_ta i=0; i<n; i++) {
    p[i]=vector_mult(map, r[i], ncls, nmodel);
    sind[i]=heapsort(p[i], ncls);
    ret[i]=sind[i][ncls-1];
  }
  base=(*skill)(truth, ret, n);

  dmap=allocate_matrix<real, int32_t>(ncls, nmodel);
  dsk=allocate_matrix<real, int32_t>(ncls, nmodel);

  //calculate "gradient" matrix:
  for (cls_t i=0; i<ncls; i++) {
    for (cls_t j=0; j<nmodel; j++) {
      ind=0;
      for (nel_ta k=0; k<n; k++) {
	real change;		//change in map
        if (ret[k]==i) {
          //how much do we have to change map element before largest prob.
	  //equals second largest prob.?
          change=(p[k][sind[k][ncls-2]]-p[k][i])/r[i][j];
	} else {
          //how much do we have to change map element before prob.
	  //equals largest prob.?
          change=(p[k][i]-p[k][ret[k]])/r[i][j];
	}
	if (fabs(change)<fabs(dmap[i][j]) || k==0) {
          dmap[i][j]=change;
	  ind=k;
	}
      }
      //calculate skill score:
      if (ret[ind]==i) {
        ret[ind]=sind[ind][ncls-2];
      } else {
        ret[ind]=i;
      }
      dsk[i][j]=base-(*skill)(truth, ret, n);
      ret[ind]=sind[ind][ncls-1];
    }
  }
}

} //end namespace libagf
