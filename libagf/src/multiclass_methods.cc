#ifndef MULTICLASS_METHODS_H
#define MULTICLASS_METHODS_H

#include <assert.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>

#include "gsl_util.h"
#include "constrained.h"
#include "randomize.h"
#include "full_util.h"

#include "agf_lib.h"

using namespace libpetey;

//methods for solving for multi-class
namespace libagf {
  template <class real>
  void prep_nonstrict(real **a, int m, int n, real *r,
		  gsl_matrix *q, gsl_vector *b) {
    //transform coding matrix to a linear system: 
    for (int i=0; i<m; i++) {
      for (int j=0; j<n; j++) {
	if (a[i][j]==0) {
          //gsl_matrix_set(q, i, j, 1);
          gsl_matrix_set(q, i, j, r[i]);
	} else {
          //gsl_matrix_set(q, i, j, map_el/r[i]);
          gsl_matrix_set(q, i, j, a[i][j]);
	}
      }
      gsl_vector_set(b, i, r[i]);
    }
  }

  //various lame-ass methods of re-normalizing conditional probabilities to
  //remove values less than 0 and ensure they sum to 1:
  //
  //zero negative values, sum and normalize:
  template <class real>
  void p_renorm(real *p, int ncls) {
    real pt=0;
    for (int i=0; i<ncls; i++) p[i]/=pt;
  }

  //zero negative values; additive normalization constant; repeat until
  //no negatives are left:
  template <class real>
  void p_constrain_renorm1(real *tly, int ncls) {
    real k;
    real pt=0;
    int ind[ncls];
    int ng=0, ng2;
    for (int i=0; i<ncls; i++) {
      if (tly[i]<0) {
        tly[i]=0;
      } else {
        pt+=tly[i];
        ind[ng]=i;
	ng++;
      }
    }
    do {
      k=(1-pt)/ng;
      ng2=ng;
      ng=0;
      pt=0;
      for (int j=0; j<ng2; j++) {
        tly[ind[j]]=tly[ind[j]]+k;
        if (tly[ind[j]]<0) {
          tly[ind[j]]=0;
        } else {
          pt+=tly[ind[j]];
          ind[ng]=ind[j];
          ng++;
        }
      }
    } while (ng!=ng2);
  }

  //another version of the same (didn't realize I'd written it already...):
  template <class real>
  void p_constrain_renorm1a(real *p, int n) {
    int ind[n];			//indices of non-zero probabilities
    int nleft=n;		//number non-zero probabilities
    int nnew0;			//number of new zeros
    int done;

    for (int i=0; i<n; i++) ind[i]=i;
    do {
      real pt;
      //first we normalize the probabilities:
      pt=0;
      for (int i=0; i<nleft; i++) pt+=p[ind[i]];
      for (int i=0; i<nleft; i++) p[ind[i]]+=(1-pt)/n;
      done=1;
      nnew0=0;
      for (int i=0; i<nleft-nnew0; i++) {
	ind[i]=ind[i+nnew0];
        if (p[ind[i]]<0) {
          p[ind[i]]=0;
          nnew0++;
	  ind[i]=ind[i+nnew0];
          done=0;
	}
      }
    } while (done==0);
  }

  //do it recursively:
  template <class real>
  void p_constrain_renorm1b(real **p, int n) {
    real pt=0;
    int nbad=0;
    //re-normalize:
    for (int i=0; i<n; i++) {
      pt+=*p[i];
    }
    for (int i=0; i<n; i++) {
      *p[i]+=(1-pt)/n;
    }
    //check for out-of-range (<0) values:
    for (int i=0; i<n; i++) {
      p[i-nbad]=p[i];
      if (*p[i]<0) {
        *p[i]=0;
        nbad++;
      }
    }
    if (nbad>0) p_constrain_renorm1b(p, n-nbad);
  }

  template <class real>
  void p_constrain_renorm1b(real *p, int n) {
    real *p2[n];
    for (int i=0; i<n; i++) p2[i]=p+i;
    p_constrain_renorm1b(p2, n);
  }

  template <class real>
  void solve_class_scratch(real **a0, int m, int n, real *r, real *p) {
    gsl_matrix *u1;
    gsl_vector *b1;
    gsl_vector *p1;

    //prepare the linear system:
    u1=gsl_matrix_alloc(m+1, n);
    b1=gsl_vector_alloc(m+1);

    prep_nonstrict(a0, m, n, r, u1, b1);

    //if problem is under-determined or otherwise unstable, this will
    //help, otherwise it makes no difference:
    for (int i=0; i<n; i++) gsl_matrix_set(u1, m, i, 1);
    gsl_vector_set(b1, m, 1);

    //print_gsl_matrix(stdout, u1);
    //gsl_vector_fprintf(stdout, b1, "%g");

    //solve the linear system:
    p1=gsl_vector_alloc(n);
    gsl_lsq_solver(u1, b1, p1);
    //gsl_vector_fprintf(stdout, p1, "%g");

    for (int i=0; i<n; i++) p[i]=gsl_vector_get(p1, i);

    //clean up:
    gsl_vector_free(p1);
    gsl_matrix_free(u1);
    gsl_vector_free(b1);

  }

  template <class real>
  void solve_class_vote(real **a, int m, int n, real *r, real *p) {
    for (int i=0; i<n; i++) p[i]=0;
    for (int i=0; i<m; i++) {
      for (int j=0; j<n; j++) {
        p[j]+=a[i][j]*r[i]/fabs(r[i]);
      }
    }
  }

  template <class real>
  void solve_class_vote_pdf(real **a, int m, int n, real *r, real *p) {
    for (int i=0; i<n; i++) p[i]=0;
    for (int i=0; i<m; i++) {
      for (int j=0; j<n; j++) {
        p[j]+=a[i][j]*r[i];
      }
    }
  }

  //here we add the constraint in the form of a "slack" variable after first
  //forming the normal equations:
  template <class real>
  void solve_class_norm2(real **a0, int m, int n, real *r, real *p) {
    real **a;
    real **at;
    real **ata;

    //since normalization constraint has been added, we can use a linear
    //system in which the result is 0:
    a=copy_matrix(a0, m, n);
    for (int i=0; i<m; i++) {
      for (int j=0; j<n; j++) {
	if (a[i][j]!=0) a[i][j]-=r[i];
      }
    }
    at=matrix_transpose(a, m, n);
    ata=matrix_mult(at, a, n, m, n);

    gsl_matrix *q=gsl_matrix_alloc(n+1, n+1);
    gsl_vector *b=gsl_vector_alloc(n+1);
    
    for (int i=0; i<n; i++) {
      for (int j=0; j<n; j++) {
        gsl_matrix_set(q, i, j, ata[i][j]);
      }
      gsl_vector_set(b, i, 0);
      gsl_matrix_set(q, i, n, 1);
      gsl_matrix_set(q, n, i, 1);
    }
    gsl_matrix_set(q, n, n, 0);
    gsl_vector_set(b, n, 1);

    gsl_vector *p1=gsl_vector_alloc(n+1);
    gsl_lsq_solver(q, b, p1);

    for (int i=0; i<n; i++) p[i]=gsl_vector_get(p1, i);

    delete_matrix(a);
    delete_matrix(at);
    delete_matrix(ata);

    gsl_matrix_free(q);
    gsl_vector_free(b);
    gsl_vector_free(p1);

  }

  //here we add the constraint in the form of a "slack" variable before
  //solving the least-squares problem:
  template <class real>
  void solve_class_norm1(real **a, int m, int n, real *r, real *p) {
    gsl_matrix *q;
    gsl_vector *b;
    //solution:
    gsl_vector *p1;

    //prepare linear system:
    q=gsl_matrix_alloc(m+1, n+1);
    gsl_matrix_set_zero(q);
    b=gsl_vector_alloc(m+1);

    //since normalization constraint has been added, we can use a linear
    //system in which the result is 0:
    for (int i=0; i<m; i++) {
      for (int j=0; j<n; j++) {
	if (a[i][j]!=0) gsl_matrix_set(q, i, j, a[i][j]-r[i]);
      }
      gsl_vector_set(b, i, 0);
    }

    //set normalization constraints:
    for (int i=0; i<m; i++) gsl_matrix_set(q, i, n, 1);
    for (int i=0; i<n; i++) gsl_matrix_set(q, m, i, 1);
    gsl_matrix_set(q, m, n, 0);
    gsl_vector_set(b, m, 1);

    //solve least squares problem:
    p1=gsl_vector_alloc(n+1);
    gsl_lsq_solver(q, b, p1);

    //extract solution:
    for (int i=0; i<n; i++) p[i]=gsl_vector_get(p1, i);

    //clean up:
    gsl_vector_free(p1);
    gsl_matrix_free(q);
    gsl_vector_free(b);

  }

  //old method:
  template <class real>
  void solve_class_constrained1(real **a, int m, int n, real *r, real *p) {
    gsl_matrix *q;
    gsl_vector *b;
    gsl_vector *p1;
    
    //prepare linear system:
    q=gsl_matrix_alloc(m+1, n);
    b=gsl_vector_alloc(m+1);
    prep_nonstrict(a, m, n, r, q, b);

    //print_gsl_matrix(stdout, q);
    //gsl_vector_fprintf(stdout, b, "%g");

    for (int i=0; i<n; i++) gsl_matrix_set(q, m, i, 1);
    gsl_vector_set(b, m, 1);

    p1=gsl_vector_alloc(n);

    solve_cond_prob(q, b, p1);

    for (int i=0; i<n; i++) p[i]=gsl_vector_get(p1, i);

    gsl_vector_free(p1);
    gsl_vector_free(b);
    gsl_matrix_free(q);

  }

  //new method:
  template <class real>
  void solve_class_constrained2(real **a, int m, int n, real *r, real *p) {
    gsl_vector *bt=gsl_vector_alloc(m);
    gsl_vector *p2=gsl_vector_alloc(n-1);
    gsl_matrix *at=gsl_matrix_alloc(m, n-1);
    //for normalization constraints:
    gsl_matrix *cnorm=gsl_matrix_alloc(n, n-1);
    gsl_vector *cthresh=gsl_vector_alloc(n);
    gsl_vector_view lastrow;

    long *ind0=randomize(n);

    //initialize the constraints:
    gsl_matrix_set_identity(cnorm);
    gsl_vector_set_zero(cthresh);
    lastrow=gsl_matrix_row(cnorm, n-1);
    gsl_vector_set_all(&lastrow.vector, -1);
    gsl_vector_set(cthresh, n-1, -1);

    for (int k=0; k<n; k++) {
      int ind=ind0[k];
      //int ind=2;
      for (int i=0; i<m; i++) {
        //double aind=gsl_matrix_get(map1, i, ind);
	double aind=a[i][ind];
	if (aind==0) aind=r[i];
        gsl_vector_set(bt, i, r[i]-aind);
        for (int j=0; j<ind; j++) {
          //gsl_matrix_set(at, i, j, gsl_matrix_get(map1, i, j)-aind);
	  if (a[i][j]==0) {
            gsl_matrix_set(at, i, j, r[i]-aind);
	  } else {
            gsl_matrix_set(at, i, j, a[i][j]-aind);
	  }
        }
        for (int j=ind+1; j<n; j++) {
          //gsl_matrix_set(at, i, j-1, gsl_matrix_get(map1, i, j)-aind);
	  if (a[i][j]==0) {
            gsl_matrix_set(at, i, j-1, r[i]-aind);
	  } else {
            gsl_matrix_set(at, i, j-1, a[i][j]-aind);
	  }
        }
      }

      //print_matrix(stdout, a, m, n);
      //printf("\n");

      //print_gsl_matrix(stdout, at);
      //printf("\n");

      constrained(at, bt, cnorm, cthresh, p2);

      //reconstitute missing variable and extract the rest:
      p[ind]=1;
      for (int j=0; j<ind; j++) {
        p[j]=gsl_vector_get(p2, j);
        p[ind]-=p[j];
      }
      for (int j=ind+1; j<n; j++) {
        p[j]=gsl_vector_get(p2, j-1);
        p[ind]-=p[j];
      }
      if (p[ind] >= 0) break;
    }

    delete [] ind0;
    gsl_vector_free(p2);
    gsl_vector_free(bt);
    gsl_matrix_free(at);
    gsl_matrix_free(cnorm);
    gsl_vector_free(cthresh);
  
  }

  template <class real>
  void solve_class_renorm(real **a, int m, int n, real *r, real *p) {
    solve_class_scratch(a, m, n, r, p);
    p_constrain_renorm1(p, n);
  }
      
  template <class real>
  void solve_class_vote_pdf2(real **a, int m, int n, real *r, real *p) {
    solve_class_vote_pdf(a, m, n, r, p);
    p_constrain_renorm1b(p, n);
  }
      
  template <class real>
  void solve_class_1vR(real **a, int m, int n, real *r, real *p) {
    printf("%d %d\n", m, n);
    assert(m<=n);
    for (int i=0; i<m; i++) p[i]=(1+r[i])/2;
  }

  template void p_renorm<float>(float *, int);
  template void p_renorm<double>(double *, int);

  template void p_constrain_renorm1<float>(float *, int);
  template void p_constrain_renorm1<double>(double *, int);

  template void p_constrain_renorm1a<float>(float *, int);
  template void p_constrain_renorm1a<double>(double *, int);

  template void solve_class_scratch<float>(float **, int, int, float *, float *);
  template void solve_class_scratch<double>(double **, int, int, double *, double *);

  template void solve_class_vote<float>(float **, int, int, float *, float *);
  template void solve_class_vote<double>(double **, int, int, double *, double *);

  template void solve_class_vote_pdf<float>(float **, int, int, float *, float *);
  template void solve_class_vote_pdf<double>(double **, int, int, double *, double *);

  template void solve_class_norm1<float>(float **, int, int, float *, float *);
  template void solve_class_norm1<double>(double **, int, int, double *, double *);

  template void solve_class_norm2<float>(float **, int, int, float *, float *);
  template void solve_class_norm2<double>(double **, int, int, double *, double *);

  template void solve_class_constrained1<float>(float **, int, int, float *, float *);
  template void solve_class_constrained1<double>(double **, int, int, double *, double *);

  template void solve_class_constrained2<float>(float **, int, int, float *, float *);
  template void solve_class_constrained2<double>(double **, int, int, double *, double *);

  template void solve_class_renorm<float>(float **, int, int, float *, float *);
  template void solve_class_renorm<double>(double **, int, int, double *, double *);

  template void solve_class_vote_pdf2<float>(float **, int, int, float *, float *);
  template void solve_class_vote_pdf2<double>(double **, int, int, double *, double *);

  template void solve_class_1vR<float>(float **, int, int, float *, float *);
  template void solve_class_1vR<double>(double **, int, int, double *, double *);

}

#endif

