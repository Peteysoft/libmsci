#ifndef _LIBAGF__SVM2CLASS__H_
#define _LIBAGF__SVM2CLASS__H_

#include "agf_defs.h"
#include "binary_classifier.h"

namespace libagf {
  template <class real, class cls_t>
  class svm2class:public binaryclassifier<real, cls_t> {
    protected:
      real **sv;		//support vectors
      real *coef;		//coefficients
      nel_ta nsv;		//number of support vectors
      real rho;			//constant term

      //kernel function:
      real (kernel *) (real *, real *, dim_ta n, void *param);
      //parameters for kernel function:
      real *param;

      //for calculating probabilities:
      real probA;
      real probB;
    public:
      svm2class(char *modfile);
      ~svm2class();
      int init(char *modfile);
      virtual real R(real *x, real *praw=NULL);
  };
}

     
