#ifndef _LIBAGF__BINARYCLASSIFIER__H
#define _LIBAGF__BINARYCLASSIFIER__H 1

#include <stdio.h>
#include <gsl/gsl_matrix.h>

//we are working in single precision for no other reason than to prove
//that my algorithms are stable:
#define SIGFUN_TYPE double

#include "classifier_obj.h"

namespace libagf {
  //much simpler and more sensible way of doing things:
  template <typename real, typename cls_t>
  binaryclassifier<real, cls_t> * binclass_init(char *name, int typecode);

  //forward declaration:
  template <typename real, typename cls_t>
  class svm2class;

  //general 2-class classifier (non-abstract) class:
  template <typename real, typename cls_t>
  class binaryclassifier:public classifier_obj<real, cls_t> {
    protected:
      int id;			//unique index simply them  
				//in the order in which they divide the classes

      int order;		//order of calibration [1]
      real *calcoef;		//calibration coefficients

      //function to transform decision value to approximate probabilities:
      SIGFUN_TYPE (*sigmoid_func) (SIGFUN_TYPE);

    public:
      binaryclassifier();
      binaryclassifier(char *nm);	//have to give it a name!
      virtual ~binaryclassifier();

      virtual real decision(real *x);	//decision function

      //returns the difference in cond. prob.
      virtual real R(real *x, 		//test point
		real *praw=NULL); 	//sticks R in correct location
      //linearly transform test point before performing classification:
      real R_t(real *x, 		//test point
		real *praw=NULL); 	//sticks R in correct location
      //the next two are defined from R:
      virtual cls_t classify(real *x, real &p, real *praw=NULL);
      virtual cls_t classify(real *x, real *p, real *praw=NULL);
      virtual void batchR(real **x, real *R, nel_ta n, dim_ta nvar);
      //linearly transform test point before performing classification:
      void batchR_t(real **x, real *R, nel_ta n, dim_ta nvar);
      //these two are defined from batchR:
      virtual void batch_classify(real **x, cls_t *cls, real *p, nel_ta n, dim_ta nvar);
      virtual void batch_classify(real **x, cls_t *cls, real **p, nel_ta n, dim_ta nvar);
      //return numerical derivative to validate above:
      void R_deriv_num(real *x, 			//test point
		      real dx, 			//absolute different in x
		      				//(dr/dx=(R(x+dx)-R(x))/dx)
		      real *drdx);		//approximate gradient of R
      virtual void print(FILE *fs, char *fbase=NULL, int depth=0);
      virtual int commands(multi_train_param &param, cls_t **clist, char *fbase);
      virtual int get_code(cls_t **clist, int **code, char **model, int &nmodel, char *fbase);
      virtual int get_code(int **code, char **model);
      virtual void set_id(cls_t *id1);
      virtual cls_t collect_binary_classifiers(binaryclassifier<real, cls_t> **list);
  };

  //if we want to use an external binary classifier:
  template <typename real, typename cls_t>
  class general2class:public binaryclassifier<real, cls_t> {
    protected:
      char *model;		//file containing model data
      char *command;		//command for doing the classification

      int Mflag;		//use LIBSVM format
      int Kflag;		//keep temporary files

      //pre-allocated "workspace" variables:
      char *syscall;		//system call
      int insert_pt;		//where in syscall to start writing feature data
      char *infile;		//pass to external classification routine
      char *outfile;		//output from external classification routine
    public:
      //initialized from a command name and a file name:
      general2class (const char *mod,		//model file
			const char *com);	//command to return predictions
      general2class (const char *mod,		//model file
			const char *com,	//command to return predictions
			int mf=0,			//use libsvm format
			int kf=0);			//keep temporary files
      virtual ~general2class();
      virtual real R(real *x, real *praw=NULL);
      virtual void batchR(real **x, real *R, nel_ta n, dim_ta nvar);
      virtual void batch_classify(real **x, cls_t *cls, real *p, nel_ta n, dim_ta nvar);
      virtual void batch_classify(real **x, cls_t *cls, real **p, nel_ta n, dim_ta nvar);
      virtual int commands(multi_train_param &param, cls_t **clist, char *fbase);
      virtual int get_code(cls_t **clist, int **code, char **model, int &nmodel,
		      char *fbase=NULL);
      virtual dim_ta n_feat();		//always -1 since classifications are
      					//done by calling an external executable
  };

  /*
  //a recursive construct like this = weird kind of neural network
  template <typename real, typename cls_t>
  class binary_calibrated<real, cls_t>:public binaryclassifier<real, cls_t> {
    protected:
      binaryclassifier<real, cls_t> *uncal;	//un-calibrated classifier

      int order;		//order 
      real *coef;		//calibration coefficients

      real (*sigfun) (real);
    public:
      virtual real R(real *x, real *praw=NULL);
      real decision(real *x);

  };
  */

}

#endif

