#ifndef MULTICLASS_H
#define MULTICLASS_H 1

#include <stdio.h>
#include <gsl/gsl_matrix.h>

#include "agf_defs.h"
#include "classifier_obj.h"
#include "binaryclassifier.h"

namespace libagf {

  //non-hierarchical multi-class classification:
  template <typename real, typename cls_t>
  class multiclass:public classifier_obj<real, cls_t> {
    private:
      //use this function to solve for the conditional probabilities:
      void (* solve_class) (real **, int, int, real *, real *);

      //set the function based on a code:
      void set_solve_type(int ct);

      //maps the actual results to the conditional probabilities:
      //(singular value decomposition of matrix mapping of final cond. prob. to raw cond. prob.)
      gsl_matrix *u;
      gsl_matrix *vt;
      gsl_vector *s;

      //internal routines:
      //run each of the binary classifiers and collect result into vector:
      void raw_classify(real *x, gsl_vector *b);

      //prepare SVD of coding matrix:
      int code_svd();

      //type of classification:
      int type;

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
					//  see init method for codes
					//  and multiclass_methods.h
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

      int init(char **fname,
		      int **code,
		      int npart,
		      cls_t ncls,
		      binaryclassifier<real, cls_t> * (* binit) (char *, void *),
		      void *param);

      //if we want to initialize with a list of files and partitions:
      int init(char **fname, 		//name of each of the binary models
		int **coding_matrix,	//coding matrix
		int npart, 		//number of partitions
		cls_t ncls,		//number of classes
      		char *prefix=NULL,	//path to data files
		int method=0,		//solution method
 		int tflag=0,		//for training purposes
		char *com=NULL,		//external binary classification command
		int Mflag=0,		//external command uses LIBSVM format
 		int Kflag=0,		//keep temporary files
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

      //these two do almost the same thing:
      //(design tree data structure, implement as generalized operator...)
      virtual void set_id(cls_t *id); 		//set id's of each binary classifier
      virtual cls_t collect_binary_classifiers(binaryclassifier<real, cls_t> **list);

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


