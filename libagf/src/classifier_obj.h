#ifndef CLASSIFIER_OBJ_H
#define CLASSIFIER_OBJ_H 1

#include <stdio.h>
#include <gsl/gsl_matrix.h>

#include "agf_defs.h"
#include "multi_parse.h"

namespace libagf {
  //"abstract" base class (now supplying default behaviour, so i guess it's
  //not so abstract after all):
  template <class real, class cls_t>
  class classifier_obj {
    protected:
      char *name;		//the file in which it is contained etc.

      cls_t ncls;		//number of classes
      dim_ta D;			//number of dimensions in test points (for original, non-transformed system only)

      real **mat;		//transformation matrix
      real *b;			//translation vector
      nel_ta D1;		//dimension of transformed system

      //real *xtran;		//transformed feature vector
      real* do_xtran(real *x);	//transforms the feature vector

    public:
      classifier_obj();
      virtual ~classifier_obj();

      //classification returning conditional prob:
      virtual cls_t classify(real *x, 	//test point
		real *p,		//returned conditional probabilities
		real *praw=NULL);	//returned "raw" conditional probabilities

      //classification returning conditional prob:
      virtual cls_t classify(real *x, 	//test point
		real &p,		//returned winning conditional probability
		real *praw=NULL);	//returned "raw" conditional probabilities

      //linear transformation on features (plus translation...):
      //(transformations are not copied but stored as references)
      int ltran(real **mat, 		//tranformation matrix
		real *b,		//constant term
		dim_ta d1,		//first dimension of trans. mat.
		dim_ta d2, 		//second 	"
		int flag);		//transform model?

      //transform the model:
      virtual int ltran_model(real **mat, 	//tranformation matrix
		real *b,		//constant term
		dim_ta d1,		//first dimension of trans. mat.
		dim_ta d2); 		//second 	"

      //copy transformation from another classifier:
      int copy_ltran(classifier_obj<real, cls_t> *other);

      //classification including linear transformation of feature space:
      cls_t classify_t(real *x, 		//test point
		real *p,		//returned conditional probabilities
		real *praw=NULL);	//returned "raw" conditional probabilities

      //classification returning conditional prob:
      //(including linear transformation of feature space)
      cls_t classify_t(real *x, 	//test point
		real &p,		//returned winning conditional probability
		real *praw=NULL);	//returned "raw" conditional probabilities

      //return list of classes (default is 0..n_cls-1):
      virtual cls_t class_list(cls_t *cls);

      //depth of hierarchy (default is 1):
      virtual int max_depth(int cur=0);

      virtual cls_t n_class();		//return number of classes (default is ncls)
      virtual dim_ta n_feat();		//return number of features for transformed problem (default is D1)
      dim_ta n_feat_t();		//return number of features for non-transformed problem

      //for use with an external command:
      //only winning probability is returned:
      virtual void batch_classify(real **x, 	//list of test points
		cls_t *cls, 			//returned classes
		real *p, 			//returned probabilities
		nel_ta n, 			//number of test points
		dim_ta nvar);			//number of variables 
      //one probability is returned for each class:
      virtual void batch_classify(real **x, 
		cls_t *cls, 
		real **p, 		//pre-allocated matrix of probabilities
		nel_ta n, 
		dim_ta nvar);

      //only winning probability is returned:
      //test point is linearly transformed first
      void batch_classify_t(real **x, 		//list of test points
		cls_t *cls, 			//returned classes
		real *p, 			//returned probabilities
		nel_ta n, 			//number of test points
		dim_ta nvar);			//number of variables 
      //one probability is returned for each class:
      //test point is linearly transformed first
      void batch_classify_t(real **x, 
		cls_t *cls, 
		real **p, 		//pre-allocated matrix of probabilities
		nel_ta n, 
		dim_ta nvar);

      //we're roping these classes in to do double duty:
      virtual void print(FILE *fs, 		//output file stream
		char *fbase=NULL,		//base name of model files
		int depth=0);			//for pretty, indented output
      //generate commands for training:
      virtual int commands(multi_train_param &param, //see structure in multi_parse.h
		cls_t **clist, 			//list of classes, partioned if necessary
		char *fbase);			//base name of model files

      //for extension to continuum retrievals (setting raw prob.):
      virtual void set_id(cls_t *id);

      //load and save files all in one as ASCII:
      virtual int load(FILE *fs);
      virtual int save(FILE *fs);

      //sure, lets do it, a general "training" method:
      virtual void train(real **train,		//features data
		      cls_t *cls,		//class labels
		      nel_ta ntrain,		//number of samples
		      int type,			//what type of training?
		      real *param);		//parameters
  };

  //lame, doesn't really do anything, but makes a lot of the methods in 
  //multiclass_hier a lot more elegant:
  //(need to do the same for the dendrogram...)
  template <class real, class cls_t>
  class oneclass:public classifier_obj<real, cls_t> {
    protected:
      cls_t cls;
    public:
      oneclass(cls_t cl);
      virtual ~oneclass();
      virtual cls_t classify(real *x, real &p, real *praw=NULL);
      virtual cls_t classify(real *x, real *p, real *praw=NULL);
      virtual int ltran_model(real **mat, real *b, dim_ta d1, dim_ta d2);
      virtual cls_t class_list(cls_t *cls);
      virtual int max_depth(int cur=0);
      virtual void batch_classify(real **x, cls_t *cls, real *p, nel_ta n, dim_ta nvar);
      virtual void batch_classify(real **x, cls_t *cls, real **p, nel_ta n, dim_ta nvar);
      virtual void print(FILE *fs, char *fbase=NULL, int depth=0);
      virtual int commands(multi_train_param &param, 
		cls_t **clist, char *fbase);
      virtual int load(FILE *fs);
      virtual int save(FILE *fs);

  };

}

#endif

