#ifndef __LIBAGF_SVM_MULTI_H__DEFINED__
#define __LIBAGF_SVM_MULTI_H__DEFINED__

#include "classifier_obj.h"

namespace libagf {
  template <class real>
  int solve_cond_prob_1v1(real **r, int ncls, real *p);

  //we can unify this at some later time with "svm2class" binary classifier
  template <class real, class cls_t>
  class svm_multi:public classifier_obj<real, cls_t> {
    protected:
      real **sv;		//support vectors
      real **coef;		//coefficients
      nel_ta nsv_total;		//total number of support vectors
      nel_ta *nsv;		//number of support vectors
      real *rho;		//constant terms

          //kernel function:
      real (* kernel) (real *, real *, dim_ta, void *);
      //parameters for kernel function:
      real *param;
      //kernel function with derivatives:
      real (* kernel_deriv) (real *, real *, dim_ta, void *, real *);

      //for calculating probabilities:
      real *probA;
      real *probB;

      //class labels:
      cls_t *label;

      real ** classify_raw(real *x);

    public:
      svm_multi();
      svm_multi(char *file);
      virtual ~svm_multi();

      virtual cls_t classify(real *x, real *p, real *praw=NULL);
      virtual cls_t class_list(cls_t *cls);

      virtual void print(FILE *fs, char *fbase=NULL, int depth=0);
      virtual int commands(multi_train_param &param, cls_t **clist, char *fbase);

  };
}

#endif

