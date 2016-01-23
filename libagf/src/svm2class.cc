#include <math.h>
#include <string.h>

#include "read_ascii_all.h"
#include "error_codes.h"
#include "full_util.h"

#include "agf_lib.h"

using namespace libpetey;

namespace libagf {

  template <class real, class cls_t>
  svm2class<real, cls_t>::svm2class() {
    classifier=NULL;
    ttype=0;
    dflag=0;
  }

  template <class real, class cls_t>
  svm2class<real, cls_t>::svm2class(char *modfile, int tc) {
    ttype=tc;
    cls_t *clist;
    classifier=new svm_multi<real, cls_t>(modfile);
    this->ncls=classifier->n_class();
    if (this->ncls != 2) {
      fprintf(stderr, "svm2class: only binary classifiers accepted (file, %s)\n", modfile);
      throw PARAMETER_OUT_OF_RANGE;
    }
    dflag=1;
    ind1=0;
    ind2=1;
    clist=new cls_t[this->ncls];
    classifier->class_list(clist);
    label1=clist[0];
    label2=clist[1];
    delete [] clist;
  }

  template <class real, class cls_t>
  svm2class<real, cls_t>::svm2class(svm_multi<real, cls_t> *other, cls_t i, cls_t j, int cflag) {
    cls_t *clist;
    cls_t oncls=other->n_class();
    if (cflag) {
      classifier=new svm_multi<real, cls_t>(other);
      dflag=1;
    } else {
      classifier=other;
      dflag=0;
    }
    this->ncls=2;
    if (i==j || i>=oncls || j>=oncls) throw PARAMETER_OUT_OF_RANGE;
    ind1=i;
    ind2=j;
    clist=new cls_t[this->ncls];
    classifier->class_list(clist);
    label1=clist[i];
    label2=clist[j];
    delete [] clist;
  }

  template <class real, class cls_t>
  svm2class<real, cls_t>::~svm2class() {
    if (dflag) delete classifier;
  }

  //some book-keeping:
  template <class real, class cls_t>
  cls_t svm2class<real, cls_t>::class_list(cls_t *list) {
    list[0]=label1;
    list[1]=label2;
    return this->ncls;
  }

  /*
  template <class real, class cls_t>
  void svm2class<real, cls_t>::R_deriv_num(real *x, real dx, real *drdx) {
    real x1[this->D];
    real x2[this->D];
    real r1, r2;
    for (dim_ta i=0; i<this->D; i++) {
      x1[i]=x[i];
      x2[i]=x[i];
    }
    for (dim_ta i=0; i<this->D; i++) {
      x1[i]-=dx;
      x2[i]+=dx;
      r1=R(x1);
      r2=R(x2);
      //printf("r1=%g; r2=%g\n", r1, r2);
      //printf("x1[%d]=%g; x2[%d]=%g\n", i, x1[i], i, x2[i]);
      drdx[i]=(r2-r1)/dx/2;
      //printf("%g ", drdx[i]);
      x1[i]=x[i];
      x2[i]=x[i];
    }
    printf("\n");
    printf("R_d_num: x=");
    for (dim_ta i=0; i<this->D; i++) printf(" %g", x[i]);
    printf("\n");
  }
  */

  template <class real, class cls_t>
  cls_t svm2class<real, cls_t>::classify(real *x, real *p, real *praw) {
    real r;
    r=R(x, praw);
    p[0]=(1-r)/2;
    p[1]=1-p[0];
    if (r<0) return label1; else return label2;
  }

  template <class real, class cls_t>
  cls_t svm2class<real, cls_t>::classify(real *x, real &p, real *praw) {
    real r;
    r=R(x, praw);
    if (r<0) {
      p=(1-r)/2;
      return label1;
    } else {
      p=(1+r)/2;
      return label2;
    }
  }

  template <class real, class cls_t>
  int svm2class<real, cls_t>::ltran_model(real **mat1, real *b1, dim_ta d1, dim_ta d2) {
    return classifier->ltran_model(mat1, b1, d1, d2);
  }

  template <class real, class cls_t>
  real svmrfunc(real *x, void *param, real *deriv) {
    bordparam<real> *p1=(bordparam<real> *) param;
    svm2class<real, cls_t> *p2=(svm2class<real, cls_t> *) p1->rparam;
    return p2->R_deriv(x, deriv);
  }

  template class svm2class<real_a, cls_ta>;

  template real_a svmrfunc<real_a, cls_ta>(real_a *, void *, real_a *);

}
