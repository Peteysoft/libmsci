#ifndef _LIBAGF__SVM2CLASS__H_
#define _LIBAGF__SVM2CLASS__H_

#include "agf_defs.h"
#include "binaryclassifier.h"
#include "svm_multi.h"

namespace libagf {

  template <class real, class cls_t>
  class svm2class:public binaryclassifier<real, cls_t> {
    protected:
      svm_multi<real, cls_t> *classifier;
      int ttype;		//is it worth adding back in?
      int dflag;		//whether or not to delete the classifier when done
      cls_t ind1;
      cls_t ind2;
      cls_t label1;		//stupid bookkeeping...
      cls_t label2;
    public:
      svm2class();
      svm2class(char *modfile, 		//file containing LIBSVM model
		      int tc=0);	//how to transform decision value:
      					//-1: f (return raw decision value)
      					//0:  1 - 2/(1+exp(A*f+B)
					//1:  -tanh(f)
      //convert 1 vs. 1 multi-class classifier into binary classifier:
      svm2class(svm_multi<real, cls_t> *other, 		//1 vs. 1 multi-class classifier
		      cls_t i, 				//index of first class
		      cls_t j,				//index of second class
		      int cflag=0);			//copy the multi-class classifier?
      virtual ~svm2class();

      virtual real R(real *x, real *praw=NULL);
      virtual cls_t classify(real *x, real *p, real *praw=NULL);
      virtual cls_t classify(real *x, real &p, real *praw=NULL);
      virtual cls_t class_list(cls_t *cls);

      //return difference in conditional prob. (R=P(2|x)-P(1|x)) plus derivative:
      real R_deriv(real *x, 			//test point
		      real *drdx);		//gradient of R
      //return numerical derivative to validate above:
      void R_deriv_num(real *x, 		//test point
		      real dx, 			//absolute different in x
		      				//(dr/dx=(R(x+dx)-R(x))/dx)
		      real *drdx);		//approximate gradient of R

  };

  template <class real, class cls_t>
  inline real svm2class<real, cls_t>::R(real *x, real *praw) {
    return classifier->R(x, ind1, ind2, praw);
  }

  template <class real, class cls_t>
  inline real svm2class<real, cls_t>::R_deriv(real *x, real *drdx) {
    return classifier->R_deriv(x, ind1, ind2, drdx);
  }

  //for training class borders:
  template <class real, class cls_t>
  real svmrfunc(real *x, void *param, real *deriv);

}

#endif

