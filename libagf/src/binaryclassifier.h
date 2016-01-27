#ifndef _LIBAGF__BINARYCLASSIFIER__H
#define _LIBAGF__BINARYCLASSIFIER__H 1

#include <stdio.h>
#include <gsl/gsl_matrix.h>

//we are working in single precision for no other reason than to prove
//that my algorithms are stable:
#define SIGFUN_TYPE double

#include "agf_defs.h"

namespace libagf {
  template <class real, class cls_t>
  class svm2class;

  //general 2-class classifier (non-abstract) class:
  template <class real, class cls_t>
  class binaryclassifier:public classifier_obj<real, cls_t> {
    protected:
      real **mat;		//transformation matrix
      int id;			//unique index simply them  
				//in the order in which they divide the classes
    public:
      binaryclassifier();
      binaryclassifier(char *nm);	//have to give it a name!
      virtual ~binaryclassifier();
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
      virtual void set_id(cls_t *id1);

  };

  //if we want to use an external binary classifier:
  template <class real, class cls_t>
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
      virtual void print(FILE *fs, char *fbase=NULL, int depth=0);
      virtual int commands(multi_train_param &param, cls_t **clist, char *fbase);
      virtual dim_ta n_feat();		//always -1 since classifications are
      					//done by calling an external executable
  };

  template <class real>
  real logistic_function(real x);

  //the "engine" upon which all the others are built:
  template <class real, class cls_t>
  class agf2class:public binaryclassifier<real, cls_t> {
    protected:
      real **brd;		//border samples
      real **grd;		//gradients at the samples
      real *gd;			//lengths of all the gradient vectors
      nel_ta n;			//number of samples
      //function to transform decision value to approximate probabilities:
      SIGFUN_TYPE (*sigmoid_func) (SIGFUN_TYPE);
      //real *ave;		//need these to condtion the test data
      void calc_grad_len();
    public:
      agf2class();
      //initialized from a pair of files
      //(<fbase>.brd contains border samples, <fbase>.bgd contains 
      //corresponding gradients):
      //
      //this version allows you to pass code directly from command line:
      agf2class(const char *fbase, 		//base name of model files
		      int sigtype=0);		//0=tanh
      						//1=erf
						//2=2/(1-exp(..))
						//3=atan
      agf2class (const char *fbase, 		//base name of model files
		SIGFUN_TYPE (*sigfun)(SIGFUN_TYPE));	//function to transform decision values
      //convert from svm binary classifier:
      agf2class(svm2class<real, cls_t> *other, 
		      real **x,			//training samples
		      cls_t *cls,		//class data
		      dim_ta nvar,		//number of variables
		      nel_ta ntrain,		//number of training samples
		      nel_ta nsamp,		//number of border samples
		      real tol,			//tolerance of border samples
		      int tflag=0);		//copy transformation

      virtual ~agf2class();
      int init(const char *fbase, SIGFUN_TYPE (*sigfun)(SIGFUN_TYPE));
      //transformation matrix is not copied, only the pointer is stored
      //--do not delete original before classifier class instance:
      virtual int ltran_model(real **mat1, real *b1, dim_ta d1, dim_ta d2);
      virtual real R(real *x, real *praw=NULL);
      virtual void print(FILE *fs, char *fbase=NULL, int depth=0);

      //load and save all in one as ASCII:
      int load(FILE *fs, int vflag=0);
      int save(FILE *fs);
  };

}

#endif

