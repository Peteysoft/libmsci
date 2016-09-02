#ifndef _LIBAGF__SVM2CLASS__H_
#define _LIBAGF__SVM2CLASS__H_

#include "agf_defs.h"
#include "svm_multi.h"

namespace libagf {

  template <class real, class cls_t>
  class svm2class:public binaryclassifier<real, cls_t> {
    protected:
      svm_multi<real, cls_t> *classifier;
      int dflag;		//whether or not to delete the classifier when done
      cls_t ind1;
      cls_t ind2;
      cls_t label1;		//stupid bookkeeping...
      cls_t label2;
    public:
      svm2class();
      svm2class(char *modfile); 		//file containing LIBSVM model
      //convert 1 vs. 1 multi-class classifier into binary classifier:
      svm2class(svm_multi<real, cls_t> *other, 		//1 vs. 1 multi-class classifier
		      cls_t i, 				//index of first class
		      cls_t j,				//index of second class
		      int cflag=0);			//copy the multi-class classifier?
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

