#ifndef _LIBAGF__SVM2CLASS__H_
#define _LIBAGF__SVM2CLASS__H_

#include "agf_defs.h"
#include "svm_multi.h"

namespace libagf {

  template <class real, class cls_t>
  class linbin:public binaryclassifier<real, cls_t> {
    protected:
      real *coef;		//coefficients
    public:
      linbin();
      linbin(char *modfile); 		//file containing LIBSVM model
      virtual ~svm2class();

      virtual real R(real *x, real *praw=NULL);
      virtual cls_t class_list(cls_t *cls);

      virtual int ltran_model(real **mat, real *b, dim_ta d1, dim_ta d2);

      //return difference in conditional prob. (R=P(2|x)-P(1|x)) plus derivative:
      real R_deriv(real *x, 			//test point
		      real *drdx);		//gradient of R

  };

  //for training class borders:
  template <class real, class cls_t>
  real svmrfunc(real *x, void *param, real *deriv);

}

#endif

