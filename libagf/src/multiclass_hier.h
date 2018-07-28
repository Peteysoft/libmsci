#ifndef MULTICLASS_HIER_H
#define MULTICLASS_HIER_H 1

#include <stdio.h>
#include <gsl/gsl_matrix.h>

#include "agf_defs.h"
#include "classifier_obj.h"
#include "binaryclassifier.h"
#include "multiclass.h"

namespace libagf {

  //hierarchical multi-class classification:
  template <typename real, typename cls_t> 
  class multiclass_hier:public classifier_obj<real, cls_t> {
    protected:
      //decides which branch to ascend:
      classifier_obj<real, cls_t> *classifier;
      //"branches":
      classifier_obj<real, cls_t> **children;

      int nonh_flag;		//only one level: not hierarchical

    public:
      multiclass_hier();
      //initialize from a control file for classification:
      multiclass_hier(const char *file, 	//control file
		int type=0,		//type of result in non-hier classification
					//(just passed to above object-class)
		char *prefix=NULL,	//path to data files
		const char *com=NULL,	//command for external binary classifier
		int Mflag=0,		//use LIBSVM format for external classifier
		int Kflag=0,		//keep temporary files
		int Zflag=0);		//use in house SVM codes

      //initialize from a control file for training only:
      multiclass_hier(const char *file, 	//control file
		int argc,		//number of option arguments
		char **argv,		//default option arguments
		int maxstacksize=MAXNOPTSTACK);	//size of option stack

      //control files passed as streams:
      multiclass_hier(FILE *fs, int type=0, char *prefix=NULL, 
		const char *com=NULL, int Mflag=0, int Kflag=0, int Zflag=0);
      multiclass_hier(FILE *fs, int argc, char **argv, 
		int maxstacksize=MAXNOPTSTACK);

      virtual ~multiclass_hier();

      //move initialization to special methods:
      //high level initialization:
      int init(FILE *fs, int type, char *prefix, const char *com,
		int Mflag, int Kflag=0, int Zflag=0);
      int init(FILE *fs, int argc, char **argv,	int maxstacksize=MAXNOPTSTACK);

      //low level initialization:
      int init(multi_parse_param &param);

      //transformation matrix is not copied, only the pointer is stored
      //--do not delete original before classifier class instance:
      virtual int ltran_model(real **mat, real *b, dim_ta d1, dim_ta d2);

      //classification:
      virtual cls_t classify(real *x, real &p, real *praw=NULL);
      virtual cls_t classify(real *x, real *p, real *praw=NULL);

      //informational:
      virtual cls_t n_class();
      virtual dim_ta n_feat();
      virtual cls_t class_list(cls_t *cls);
      virtual int max_depth(int cur=0);

      //"batch" routines:
      virtual void batch_classify(real **x, cls_t *cls, real *p, nel_ta n, dim_ta nvar);
      virtual void batch_classify(real **x, cls_t *cls, real **p, nel_ta n, dim_ta nvar);

      //parsing:
      virtual void print(FILE *fs, char *fbase=NULL, int depth=0);
      //more user friendly interface to commands method:
      int generate_commands(FILE *fs, 		//write commands to file stream
		char *train,			//file containing training data
		char *fbase,			//base file name for model data
		char *command=NULL,		//trains binary classifier
		char *partcom=NULL,	//command for partitioning classes
		char *concom=NULL,	//command for file conversion
		int Kflag=0,		//keep temporary files
		unsigned long session_id=0);
      void generate_commands2(FILE *fs, 	//write commands to file stream
		char *train,			//file containing training data
		char *fbase,			//base file name for model data
		char *command=NULL,		//trains binary classifier
		char *partition=NULL,	//command for partitioning classes
		char *fconv=NULL,	//command for file conversion
		int Kflag=0,		//keep temporary files
		unsigned long session_id=0);
      virtual int commands(multi_train_param &param, cls_t **clist,
		char *fbase);
      void get_code(int **&code, char **&model, int &nmodel, cls_t &ncls);
      virtual int get_code(int **code, char **model);

      virtual void set_id(cls_t *id);
      //collect the binary classifiers into a single list:
      virtual cls_t collect_binary_classifiers(binaryclassifier<real, cls_t> **list);

      //load and save to single ASCII file:
      int load(FILE *fs, int ct);
      virtual int load(FILE *fs);
      virtual int save(FILE *fs);

      //one more ugly hack:
      virtual void train(real **train,
		      cls_t *cls,
		      nel_ta ntrain,
		      int type,
		      real *param);
  };
}

#endif

