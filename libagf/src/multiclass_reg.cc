#include <stdio.h>
#include <getopt.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <assert.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>

#include "full_util.h"
#include "peteys_tmpl_lib.h"
#include "gsl_util.h"

#include "agf_lib.h"

using namespace std;
using namespace libpetey;

//regularized solutions for non-hierarchical multiclass class:

namespace libagf {

#define CONSTRAINT_WEIGHT 1.
#define CONSTRAINT_WEIGHT2 10.
#define CONSTRAINT_WEIGHT3 2.

  template <class real, class cls_t, class twoclass_t>
  multiclass<real, cls_t, twoclass_t>::init2(char *fname, int **part, int npart) {
    gsl_vector *work;

    cls_ta count;

    this->ncls=0;
    for (int i=0; i<nmodel; i++) {
      for (int j=0; part[i*2][j]>=0; j++) {
        if (part[i*2][j]>=this->ncls) this->ncls=part[i*2][j]+1;
      }
      for (int j=0; part[i*2+1][j]>=0; j++) {
        if (part[i*2+1][j]>=this->ncls) this->ncls=part[i*2+1][j]+1;
      }
    }
    nmodel=npart;

    twoclass=new binaryclassifier<real, cls_t> *[nmodel];
    for (cls_t i=0; i<nmodel; i++) {
      //twoclass[i]=new agf2class<real, cls_t>(fname[i]);
      twoclass[i]=new twoclass_t(fname[i]);
    }

    //create the mapping:
    map=gsl_matrix_alloc(nmodel, this->ncls-1);
    gsl_matrix_set_zero(map);

    //sum of the conditional probabilities on one side of the partition
    //is equal to the conditional probability returned from the 2-class
    //classification result:
    for (int i=0; i<nmodel*2; i++) {
      double sgn=2*(i%2)-1;
      for (int j=0; part[i][j]>=0; j++) gsl_matrix_set(map, i/2, part[i][j], sgn);
    }

    /*
    printf("multiclass: pdf inverse mapping:\n");
    for (int i=0; i<nmodel+1; i++) {
      for (int j=0; j<this->ncls; j++) printf("%5.2lf ", gsl_matrix_get(map, i, j));
      printf("\n");
    }
    */
    //gsl_matrix_fprintf(stdout, map, "%5.2f");

    //now we find the inverse of this matrix:
    u=gsl_matrix_alloc(nmodel, this->ncls-1);
    gsl_matrix_memcpy(u, map);

    for (int i=0; i<nmodel; i++) {
      real val=-gsl_matrix_get(map, i, this->ncls-1);
      for (int j=0; j<this->ncls-1; j++) {
        val+=gsl_matrix_get(map, i, j);
        gsl_matrix_set(u, i, j, val);
      }
    }

    //print_gsl_matrix(stdout, u);
    vt=gsl_matrix_alloc(this->ncls-1, this->ncls-1);
    s=gsl_vector_alloc(this->ncls-1);
    work=gsl_vector_alloc(this->ncls-1);
    gsl_linalg_SV_decomp(u, vt, s, work);
    //gsl_linalg_SV_decomp_jacobi(u, vt, s);

    /*
    printf("U:\n");
    print_gsl_matrix(stdout, u);
    printf("S:\n");
    for (int i=0; i<s->size; i++) printf("%10.5g ", gsl_vector_get(s, i));
    printf("\nV^T:\n");
    print_gsl_matrix(stdout, vt);
    */

    gsl_vector_free(work);

    for (int i=0; i<nmodel; i++) delete [] fname[i];
    for (int i=0; i<nmodel*2; i++) delete [] part[i];
  
  }

  template <class real, class cls_t, class twoclass_t>
  cls_t multiclass<real, cls_t, twoclass_t>::classify_reg(real *x, real *pdf) {
    gsl_vector *b;
    real *b2;
    gsl_vector *p;
    cls_t cls1, cls2;
    real pt, pmin;
    real *tally;

    b=gsl_vector_alloc(nmodel+2*this->ncls+1);
    b2=new real[nmodel];
    p=gsl_vector_alloc(this->ncls);

    //printf("multiclass raw pdfs: ");
    for (int i=0; i<nmodel; i++) {
      real R=twoclass[i]->R(x);
      //gsl_vector_set(b, i, R-gsl_matrix_get(map, i, this->ncls-1));
      gsl_vector_set(b, i, R);
      b2[i]=R;
    }
    //printf("\n");
    gsl_vector_set(b, nmodel, CONSTRAINT_WEIGHT);

    //gsl_linalg_SV_solve(u, vt, s, b, p);

    tally=new real[this->ncls];
    //vote_label(b2, pdf);
    vote_pdf(b2, tally);

    //use the voting solution to regularize the matrix:
    gsl_matrix *a_p=gsl_matrix_alloc(nmodel+2*this->ncls+1, this->ncls);
    gsl_matrix *v_p=gsl_matrix_alloc(this->ncls, this->ncls);
    gsl_vector *wk_p=gsl_vector_alloc(this->ncls);
    gsl_vector *s_p=gsl_vector_alloc(this->ncls);
   
    //apply regularization:
    gsl_matrix_set_zero(a_p);
    for (int i=0; i<nmodel+1; i++) {
      for (int j=0; j<this->ncls; j++) {
        gsl_matrix_set(a_p, i, j, gsl_matrix_get(map, i, j));
      }
    } 
    for (int i=0; i<this->ncls; i++) {
      real tlysqr=tally[i]*tally[i];
      gsl_matrix_set(a_p, i+nmodel+1, i, tlysqr/CONSTRAINT_WEIGHT2);
      gsl_matrix_set(a_p, i+nmodel+1+this->ncls, i, 1/tlysqr/CONSTRAINT_WEIGHT3);
      gsl_vector_set(b, i+nmodel+1, tlysqr/CONSTRAINT_WEIGHT2);
      gsl_vector_set(b, i+nmodel+1+this->ncls, 0.);
    }
    gsl_linalg_SV_decomp(a_p, v_p, s_p, wk_p);
    gsl_linalg_SV_solve(a_p, v_p, s_p, b, p);

    gsl_matrix_free(a_p);
    gsl_matrix_free(v_p);
    gsl_vector_free(wk_p);
    gsl_vector_free(s_p);

    //reconstitute last cond. prob.:
    //pdf[this->ncls-1]=1;
    for (int i=0; i<this->ncls; i++) {
      pdf[i]=gsl_vector_get(p, i);
      //pdf[this->ncls-1]-=pdf[i];
    }

    cls1=0;
    cls2=0;
    for (int i=0; i<this->ncls; i++) {
      //pdf[i]=0;
      if (pdf[i]>pdf[cls1]) cls1=i;
      if (tally[i]>tally[cls2]) cls2=i;
    }

    //if voting is different from matrix inversion, use the voting result:
/*    if (cls1!=cls2) {
      cls1=cls2;
      for (int i=0; i<this->ncls; i++) pdf[i]=tally[i];
    }
*/

    //renormalize the conditional prob.:
    pt=0;
    pmin=0;
    for (int i=0; i<this->ncls; i++) {
      //if (pdf[i]<pmin) {
      if (pdf[i]<0) {
        //pmin=pdf[i];
        pdf[i]=0;
      }
      pt+=pdf[i];
    }
    printf("pt=%g\n", pt);

    //correct the resultant conditional probabilities:
    pt-=this->ncls*pmin;
    for (int i=0; i<this->ncls; i++) pdf[i]=(pdf[i]-pmin)/pt;

    //for (int i=0; i<this->ncls; i++) printf("%g ", pdf[i]);
    //printf("\n");

    gsl_vector_free(b);
    delete [] b2;
    gsl_vector_free(p);
    delete [] tally;

    return cls1;
  }

}

