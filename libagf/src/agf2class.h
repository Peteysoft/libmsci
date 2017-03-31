#ifndef _LIBAGF__AGF2CLASS__H
#define _LIBAGF__AGF2CLASS__H 1

#include <stdio.h>
#include <gsl/gsl_matrix.h>

//we are working in single precision for no other reason than to prove
//that my algorithms are stable:
#define SIGFUN_TYPE double

#include "binaryclassifier.h"
#include "direct_classifier.h"

namespace libagf {
  template <class real, class cls_t>
  class svm2class;

  template <typename real, typename cls_t>
  class agf2class:public binaryclassifier<real, cls_t> {
    private:
      //diagnostics:
      real min_f;
      real max_f;
      real total_f;
      real min_W;
      real max_W;
      real total_W;
      iter_ta min_nd;
      iter_ta max_nd;
      iter_ta total_nd;
      nel_ta ntrial;
      nel_ta ntrial_k;
    protected:
      real *unsort;		//unsorted training data
      real **trn0;		//base of pointer
      real **train;		//features data sorted by class
      nel_ta ntrain;		//number of training samples
      nel_ta cind;		//index of 0 to 1 transition in class labels
      real W;			//W parameter
      nel_ta k;			//k-nearest-neighbours to use in calculation
      real var0[2];		//variance brackets (unmodified)
      real var[2];		//variance brackets
    public:
      agf2class(agf_classifier<real, cls_t> *other,	//multi-class agf classifier
		      cls_t *part, cls_t npart);	//class partition
      ~agf2class();
      //difference in conditional probabilities:
      virtual real R(real *x, real *praw=NULL);
      //difference in conditional prob. plus derivatives:
      real R_deriv(real *x, 			//test point
		      real *drdx);		//gradient vector
      //reset filter variance brackets:
      void reset_var();
      //initialize parameters for training class borders:
      int param_init(nel_ta ns,				//desired number of border samples
	      int (*sfunc)(void *, real *, real *),	//function to sample opposite classes
	      bordparam<real> *param);			//parameters to pass to sfunc

  };

  //the "engine" upon which all the others are built:
  template <class real, class cls_t>
  class borders_classifier:public binaryclassifier<real, cls_t> {
    protected:
      real **brd;		//border samples
      real **grd;		//gradients at the samples
      real *gd;			//lengths of all the gradient vectors
      nel_ta n;			//number of samples
      //function to transform decision value to approximate probabilities:
      SIGFUN_TYPE (*sigmoid_func) (SIGFUN_TYPE);
      void calc_grad_len();
    public:
      borders_classifier();
      //initialized from a pair of files
      //(<fbase>.brd contains border samples, <fbase>.bgd contains 
      //corresponding gradients):
      //
      //this version allows you to pass code directly from command line:
      borders_classifier(const char *fbase, 		//base name of model files
		      int sigtype=0);		//0=tanh
      						//1=erf
						//2=2/(1-exp(..))
						//3=atan
      borders_classifier (const char *fbase, 		//base name of model files
		SIGFUN_TYPE (*sigfun)(SIGFUN_TYPE));	//function to transform decision values

      //convert from AGF classifier:
      void train(agf2class<real, cls_t> *other, 
		      nel_ta nsamp,		//number of border samples
		      real tol);		//tolerance of border samples

      //convert from svm binary classifier:
      void train(svm2class<real, cls_t> *other, 
		      real **x,			//training samples
		      cls_t *cls,		//class data
		      dim_ta nvar,		//number of variables
		      nel_ta ntrain,		//number of training samples
		      nel_ta nsamp,		//number of border samples
		      real tol,			//tolerance of border samples
		      int tflag=0);		//copy transformation

      virtual ~borders_classifier();
      int init(const char *fbase, SIGFUN_TYPE (*sigfun)(SIGFUN_TYPE));
      //transformation matrix is not copied, only the pointer is stored
      //--do not delete original before classifier class instance:
      virtual int ltran_model(real **mat1, real *b1, dim_ta d1, dim_ta d2);
      virtual real R(real *x, real *praw=NULL);
      virtual void print(FILE *fs, char *fbase=NULL, int depth=0);

      //load and save all in one as ASCII:
      virtual int load(FILE *fs);
      int load(FILE *fs, int vflag);
      virtual int save(FILE *fs);
  };

  template <class real, class cls_t>
  class borders_calibrated:public borders_classifier<real, cls_t> {
    protected:
      //calibration coefficients:
      real *coef;
      int order;
    public:
      borders_calibrated();
      borders_calibrated(const char *fbase);
      ~borders_calibrated();

      //compare calculated probabilities to actual and derive 
      //calibration coefficients:
      void calibrate(real **train, 		//training data
		      cls_t *cls, 
		      nel_ta ntrain, 		//number of samples
		      int O=3,			//order
		      int nhist=10);		//number of divisions

      virtual real R(real *x, real *praw);
  };
}

#endif

