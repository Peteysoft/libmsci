#ifndef MULTICLASS_H
#define MULTICLASS_H 1

#include <stdio.h>
#include <gsl/gsl_matrix.h>

#include "agf_defs.h"
#include "classifier_obj.h"
#include "binaryclassifier.h"

namespace libagf {

  //non-hierarchical multi-class classification:
  template <class real, class cls_t>
  class multiclass:public classifier_obj<real, cls_t> {
    private:
      //use this function to solve for the conditional probabilities:
      void (* solve_class) (real **, int, int, real *, real *);

      //maps the actual results to the conditional probabilities:
      //(singular value decomposition of matrix mapping of final cond. prob. to raw cond. prob.)
      gsl_matrix *u;
      gsl_matrix *vt;
      gsl_vector *s;

      //pre-compute constraint coefficients:
      gsl_matrix *cnorm;
      gsl_vector *cthresh;

      //internal routines:
      //run each of the binary classifiers and collect result into vector:
      void raw_classify(real *x, gsl_vector *b);

      //old methods:
      cls_t solve_class_old(gsl_vector *b, real *p1);

      //basic least squares solution:
      cls_t classify_basic(gsl_vector *b, real *p);
      cls_t vote_label(gsl_vector *b, real *tly);
      cls_t vote_pdf(gsl_vector *b, real *tly);

      //special 1 versus rest solution:
      cls_t classify_1vR(gsl_vector *b, real *p);
      //1 versus 1 solution:
      cls_t classify_1v1(gsl_vector *b, real *p);

      //more experimental versions:
      //least squares including "non-strict" coding matrices:
      cls_t classify_scratch(gsl_vector *b, real *p);
      //constrained least squares (not brute force but can return sub-optimal
      //sol'n):
      cls_t classify_special(gsl_vector *b, real *p);
      //assumes orthogonal coding matrix (multiple re-normalization steps):
      cls_t vote_pdf2(gsl_vector *b, real *tly);

      //prepare SVD of coding matrix:
      int code_svd();

      //initialize constraints:
      void init_constraint();

      //type of classification:
      int type;

      //weight for normalization constraint:
      real constraint_weight;

      //"polarity" of each binary classifier:
      int *pol;
    protected:
      //the partitions:
      binaryclassifier<real, cls_t> **twoclass;
      int nmodel;		//number of partitions

      //raw mapping:
      gsl_matrix *map;
      //in two formats:
      real ** code;		//(want to keep only this format...)

      //is the mapping "strict"--i.e. each row includes all classes?
      int strictflag;
    public:
      multiclass(int ct=0);
      //initialized from a control file:
      multiclass(const char *file, 	//control file
		int clstyp=0,		//type of classification result
					//  0=constrained inverse;
					//  1=raw inverse;
					//  2=voting from probabilities; 
					//  3=voting from classes
					//  4=inverse/voting from probabilities
					//  5=inverse/voting from classes
					//  6=inverse weighted by raw prob.
					//  7=constrained inverse (may be inefficient)
					//  8=voting from pdf, corrected and normalized
					//    (designed for orthogonal coding matrix)
					// 10=special for 1 v rest
					// 11=special for 1 v. 1
					// 12=special for adjacent
		const char *com=NULL,	//command for binary classifier
		int Mflag=0,		//LIBSVM format for external classifiers
		int Kflag=0,		//keep temporary files
		int sigcode=0,		//code for sigmoid func. to transform
     					//decision values
		int Zflag=0);		//use in house SVM codes
		
      virtual ~multiclass();

      //int init_standard(char **fname, cls_t **part, int npart, const char *com=NULL, int mf=0, int kf=0);

      //initialize using giant structure used by all the other parsing routines:
      int init(multi_parse_param &param);

      //if we want to initialize with a list of files and partitions:
      int init(char **fname, 		//name of each of the binary models
		cls_t **part, 		//partitions
		int npart, 		//number of partitions
      		char *prefix=NULL,	//path to data files
 		int tflag=0,		//for training purposes
		char *com=NULL,		//external binary classification command
		int Mflag=0,		//external command uses LIBSVM format
 		int Kflag=0,		//keep temporary files
		int sigcode=0,		//code for sigmoid trans. func.
		int Zflag=0);		//use in house SVM codes

      //transformation matrix is not copied, only the pointer is stored
      //--do not delete original before classifier class instance:
      virtual int ltran_model(real **mat, real *b, dim_ta d1, dim_ta d2);
      virtual cls_t classify(real *x, real *p, real *praw=NULL);
      virtual dim_ta n_feat();

      //disastrous "regularized" version want to keep around for a while:
      //virtual cls_t classify_reg(real *x, real *pdf);

      virtual void batch_classify(real **x, cls_t *cls, real **p, nel_ta n, dim_ta nvar);

      //for parsing:
      virtual void print(FILE *fs, char *fbase=NULL, int depth=0);
      virtual int commands(multi_train_param &param, cls_t **clist, char *fbase);

      virtual void set_id(cls_t *id); 		//set id's of each binary classifier

      int detect_type();	//detects if it is one of the special cases:
      				//0 = 1v1
				//1 = 1 vs. rest
				//2 = adjacent

      virtual int load(FILE *fs);
      virtual int save(FILE *fs);

      //fracking hack:
      virtual void train(real **train, cls_t *cls, nel_ta ntrain,
		      int type, real *param);
  };

}

#endif


