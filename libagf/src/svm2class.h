#ifndef _LIBAGF__SVM2CLASS__H_
#define _LIBAGF__SVM2CLASS__H_

#include "agf_defs.h"
#include "svm_multi.h"

namespace libagf {

  //new paradigm:
  template <typename real>
  class svm_helper {
    protected:
      nel_ta nsv;		//total support
      real **sv;		//support vectors
      dim_ta D;			//dimension of features data
      real **kval;		//kernel values
      real *test;		//current test point

      //in order for this optimization to work, all the binary classifiers
      //must have been trained with the same parameters:
      real *param;			//parameters for kernel function
      //kernel function:
      real (* kernel) (real *, real *, dim_ta, void *);
      real (* kernel_deriv) (real *, real *, dim_ta, void *, real *);
    public:
      svm_helper(FILE *fs);
      ~svm_helper();

      void register_point(real *x);
      real get_kernel(nel_ta index);
      int ltran_model(real **mat, real *b, dim_ta d1, dim_ta d2);
  };

  template <typename real, typename cls_t>
  class svm2class2:public binaryclassifier<real, cls_t> {
    protected:
      svm_helper<real> *helper;		//contains support vectors
      nel_ta *ind;			//indexes into the support vectors
      real *coef;			//coefficients
      real probA;			//for calculating probabilities
      real probB;

      //kernel function:
      real (* kernel) (real *, real *, dim_ta, void *);
      real (* kernel_deriv) (real *, real *, dim_ta, void *, real *);
    public:
      svm2class2(char *name, void *param);
      ~svm2class2();

      virtual real R(real *x, real *praw=NULL);

      virtual int ltran_model(real **mat, real *b, dim_ta d1, dim_ta d2);

      real R_deriv(real *x, 			//test point
		      real *drdx);		//gradient of R
  };

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
      svm2class(svm_multi<real, cls_t> *other, 	//1 vs. 1 multi-class classifier
		      cls_t i, 			//index of first class
		      cls_t j,			//index of second class
		      int cflag=0);		//copy the multi-class classifier?
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

