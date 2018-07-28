#ifndef MULTI_PARSE_H
#define MULTI_PARSE_H 1

#include <stdio.h>
#include <gsl/gsl_matrix.h>

#include "randomize.h"

#include "agf_defs.h"

//default stack limit:
#define MAXNOPTSTACK 10001

//because we love arbitrary limitations:
#define MAXNPART 10000

namespace libagf {

  //for generating training commands:
  struct multi_train_param {
    char *train;		//base name for training data
    char *commandname;		//command name
    FILE *commandfs;		//write commands to this file stream (stdout, normally)

    //optional parameters:
    char *partcom;		//command for partitioning class labels (optional)
    int32_t session_id;		//for creating unique names for temporary files
    char *concom;		//optional command for file conversion
    char *precom;		//optional command: use current model to generate new one
    int Kflag;				//keep temporary files
  };

  //for parsing control files:
  struct multi_parse_param {
    FILE *infs;				//control file input stream
    int lineno;				//line number
    int trainflag;			//1=training, 0=classification

    char *prefix;		//add this prefix to the file names/options

    //these parameters used for training only:
    //char *train;			//base name for training data
    char *commandname;			//command name
    //FILE *outfs;			//output control file stream
    //FILE *commandfs;			//write commands to this file stream (stdout, normally)
    //vector<char *> optstack;		//list of previous options
    char **optstack;			//list of previous options
    int stackptr;
    int maxnstack;

    int depth;		//stack depth--not used

    //for use in initializing class structures (others will be phased out):
    int type;			//type of result in non-hierarchical multi-class
    int Mflag;			//LIBSVM format
    int Kflag;			//keep temporary files
    int Zflag;			//use "in house" SVM codes
    int sigcode;	//code for sigmoid func. to transform decision values

    //mother of all hacks:
    //binaryclassifier<real, cls_t> * (* bininit) (char *name, void *options);
  };

  //initializes internal values and concatinates bottom-level options;
  //names and file streams must be initialized directly:
  //(e.g. we may want the output control file to go to stdout instead of
  // the training commands, etc...)
  void multi_train_begin(multi_parse_param *param,
		int argc,			//# command line arguments to pass
		char **argv,			//command line options
		int maxstacksize=MAXNOPTSTACK);	//I know, the stack is just an array

  void multi_train_end(multi_parse_param *param);

  template <typename code_t>
  void parse_multi_partitions(multi_parse_param *param, 	//parameters
		  char** &model, 		//name of each binary classifier
		  code_t ** &code, 		//coding matrix
		  int &nmodel, 			//number of models found
		  int &ncls);			//number of classes

  template <typename code_t>
  int test_parse_multi_partitions(int nmodel, int ncls);

  //parse the partitions for non-hierarchical multi-class:
  template <class cls_t>
  int parse_multi_partitions(multi_parse_param *param, 	//parameters--see above
		char **model, 			//list of binary partition models
		cls_t **partition, 		//classes in each partition
		int maxn);			//parse at most maxn partitions (array size)

  //test for "strict" partitioning: same number of classes in each partition, 
  //labels cannot exceed size of partitions
  template <class cls_t>
  int multi_partition_strict(
		char **model, 			//list of binary partition models
		cls_t **partition, 		//classes in each partition
		int npart);			//parse at most maxn partitions (array size)

  //these are needed internally and to initialize the classifier objects:
  int scan_tochar(FILE *fs, const char *charlist, int &lineno);
  int scan_nowhitespace(FILE *fs, int &lineno);
  char *parse_multi_start(multi_parse_param *param, int &flag, char *&options);
  char *scan_class_model1(FILE *fs, int &lineno);
  char *scan_class_model2(FILE *fs, int &lineno);
  char *scan_class_label(FILE *fs, int &lineno);

  //for debugging:
  template <class cls_t>
  void print_partition_list(FILE *fs, cls_t **part, cls_t npart);

}

#endif

