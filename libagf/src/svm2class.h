#ifndef _LIBAGF__SVM2CLASS__H_
#define _LIBAGF__SVM2CLASS__H_

#include "agf_defs.h"
#include "binaryclassifier.h"

namespace libagf {

  template <class real, class cls_t>
  class svm_multi;

  template <class real, class cls_t>
  class svm2class:public binaryclassifier<real, cls_t> {
    //*sniff*, *sniff*...
    friend svm_multi<real, cls_t>;

    protected:
      real **sv;		//support vectors
      real *coef;		//coefficients
      nel_ta nsv;		//number of support vectors
      real rho;			//constant term

      //kernel function:
      real (* kernel) (real *, real *, dim_ta, void *);
      //parameters for kernel function:
      real *param;
      //kernel function with derivatives:
      real (* kernel_deriv) (real *, real *, dim_ta, void *, real *);

      //for calculating probabilities:
      real probA;
      real probB;

      //book-keeping:
      cls_t label1;
      cls_t label2;

      //type of transformation for decision value:
      int ttype;
    public:
      svm2class();
      svm2class(char *modfile, 		//file containing LIBSVM model
		      int tc=0);	//how to transform decision value:
      					//-1: f (return raw decision value)
      					//0:  1 - 2/(1+exp(A*f+B)
					//1:  -tanh(f)
      virtual ~svm2class();
      int init(char *modfile, int tc=0);

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

  //for training class borders:
  template <class real, class cls_t>
  real svmrfunc(real *x, void *param, real *deriv);

}

#endif

