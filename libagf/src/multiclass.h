#ifndef MULTICLASS_H
#define MULTICLASS_H 1

#include <stdio.h>
#include <gsl/gsl_matrix.h>

#include "agf_defs.h"
#include "classifier_obj.h"
#include "binaryclassifier.h"

namespace libagf {
  //write the most common kinds of control files:
  void one_against_all(FILE *fs, int ncls, const char *options=NULL);
  void partition_adjacent(FILE *fs, int ncls, const char *options=NULL);
  void random_coding_matrix(FILE *fs, int ncls, int ntrial, int strictflag=0);
  void exhaustive_coding_matrix(FILE *fs, int ncls);

  //non-hierarchical multi-class classification:
  template <class real, class cls_t>
  class multiclass:public classifier_obj<real, cls_t> {
    private:
      //maps the actual results to the conditional probabilities:
      //(singular value decomposition of matrix mapping of final cond. prob. to raw cond. prob.)
      gsl_matrix *u;
      gsl_matrix *vt;
      gsl_vector *s;

      //empirically derived mapping (from raw to final):
      real **imap;

      //internal routines:
      void raw_classify(real *x, gsl_vector *b);
      cls_t classify_basic(gsl_vector *b, real *p);
      cls_t vote_label(gsl_vector *b, real *tly);
      cls_t vote_pdf(gsl_vector *b, real *tly);
      cls_t solve_class(gsl_vector *b, real *p1);

      //more experimental versions:
      cls_t classify_map(gsl_vector *b, real *p);
      cls_t classify_scratch(gsl_vector *b, real *p);
      cls_t classify_special(gsl_vector *b, real *p);

      //prepare SVD of coding matrix:
      int code_svd();

      //type of classification:
      int type;

      //weight for normalization constraint:
      real constraint_weight;

    protected:
      //the partitions:
      binaryclassifier<real, cls_t> **twoclass;
      int nmodel;		//number of partitions

      //raw mapping:
      gsl_matrix *map;

    public:
      multiclass();
      //initialized from a control file:
      multiclass(const char *file, 	//control file
		int clstyp=0,		//type of classification result
					// -1=derived ("empirical") mapping
					//  0=constrained inverse;
					//  1=raw inverse;
					//  2=voting from probabilities; 
					//  3=voting from classes
					//  4=inverse/voting from probabilities
					//  5=inverse/voting from classes
					//  6=inverse weighted by raw prob.
		real cw=1.,		//constraint weight (sum of cond. prob.)
		const char *com=NULL,	//command for binary classifier
		int Mflag=0,		//LIBSVM format for external classifiers
		int Kflag=0);		//keep temporary files
		
      virtual ~multiclass();

      //int init_standard(char **fname, cls_t **part, int npart, const char *com=NULL, int mf=0, int kf=0);

      //initialize using giant structure used by all the other parsing routines:
      int init(multi_parse_param &param);

      //if we want to initialize with a list of files and partitions:
      int init(char **fname, 		//name of each of the binary models
		      cls_t **part, 	//partitions
		      int npart, 	//number of partitions
		      int tflag=0,	//for training purposes
		      char *com=NULL,	//external binary classification command
		      int Mflag=0,	//external command uses LIBSVM format
		      int Kflag=0);	//keep temporary files

      //use actual examples to train the mapping:
      int train_map(real **train, cls_t *cls, nel_ta n);
      int write_map(FILE *fs);
      int load_map(FILE *fs);

      //transformation matrix is not copied, only the pointer is stored
      //--do not delete original before classifier class instance:
      virtual int ltran(real **mat, real *b, dim_ta d1, dim_ta d2, int flag);
      virtual cls_t classify(real *x, real *p, real *praw=NULL);
      virtual dim_ta n_feat();

      //disastrous "regularized" version want to keep around for a while:
      //virtual cls_t classify_reg(real *x, real *pdf);

      virtual void batch_classify(real **x, cls_t *cls, real **p, nel_ta n, dim_ta nvar);

      //for parsing:
      virtual void print(FILE *fs, char *fbase=NULL, int depth=0);
      virtual int commands(multi_train_param &param, cls_t **clist, char *fbase);

      virtual void set_id(cls_t *id); 		//set id's of each binary classifier
  };

}

#endif


