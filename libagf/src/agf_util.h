//
// This software is released under the following terms:
//
// 1. No commercial use.
// 2. Copies and derivative works are free to use and modify.
// 3. Attribution must be given to all contributors of both original and derivative works.
//
// Authors:
//
// 2017-07-16 Peter Mills: added license information 
//

//
// Math and API (such as it is...) utilities for libAGF
//

#ifndef AGF_UTIL_H_INCLUDED
#define AGF_UTIL_H_INCLUDED 1

#include <stdio.h>

#include "gsl/gsl_rng.h"
#include "randomize.h"

#include "agf_defs.h"
#include "agf_metric2.h"

namespace libagf {

  //random number generator:
  extern const gsl_rng_type *agf_gsl_rng_type;

  //to prevent the gsl library from aborting if it encounters an error:
  void agf_gsl_handler(const char *reason, const char * file, int line, int gsl_errno);

  //this rubbish should be redone...
  struct agf_command_opts {
    char *normfile;	//-a file containing normalization values
    char *ofile;	//-o output file
    char *multicommand;	//-O command to use for multi-class class.
    char *path;		//-y path
    real_a missing;	//-E value for missing data
    real_a ftest;	//-f fraction of data to use for testing
    real_a hrel;	//-h relative difference for numerical derivatives
    real_a pmin;	//-p threshold prob for clustering algorithm
    real_a rthresh;	//-r location of class border
    real_a tol;		//-t tolerance of border samples
    real_a var[2];	//-v 1st variance bracket
  			//-V 2nd variance bracket/initial filter variance
    real_a W1;		//-w minimum total weight
    real_a W2;		//-W target/maximum total weight
    real_a pratio;	//-X ratio between sample classes
    cls_ta cl_thresh;	//-T class threshold: multiclass->2 class
    dim_ta svd;		//-S perform singular-value-decomposition (number to keep)
    nel_ta k;		//-k number of nearest-neighbours
    nel_ta n;		//-s number of border samples
    nel_ta nt;		//-q number of trials (AGF opt)
    nel_ta div;		//-d number of divisions ~= -f
    enum_a algtype;	//-c algorithm to use in n-fold cross-validation
    enum_a metrictype;	//-m type of metric
    enum_a Qtype;	//-Q algorithm selection 2
    flag_a stdinflag;	//-0 read from stdin
    flag_a stdoutflag;	//-1 write to stdout
    flag_a asciiflag;	//-A operate on ascii files from stdin & stdout
    flag_a Bflag;	//-B sort by ordinate
    flag_a Cflag;	//-C no class data
    flag_a errflag;	//-e return error estimates for interpolation
    flag_a fflag;	//-f
    flag_a selectflag;	//-F select features
    flag_a Gflag;	//-G "non-strict" partitioning in non-hier multi-class
    flag_a Hflag;	//-H no header/strip header
    flag_a jointflag;	//-j return joint probabilities (instead of cond.) to sdout
    flag_a Kflag;	//-K keep temporary files
    flag_a Lflag;	//-L floating point ordinates
    flag_a Mflag;	//-M SVM format
    flag_a normflag;	//-n flag to normalize data
    flag_a Nflag;	//-N take data from stdin
    flag_a Pflag;	//-P calculate cross-correlation matrix
    flag_a Rflag;	//-R randomly sample
    flag_a uflag;	//-u un-normalize borders and gradients
    flag_a Uflag;	//-U re-label classes to go from [0-nc).
    flag_a xflag;	//-x run in background
    flag_a Yflag;	//-Y top row of orthogonal control matrix implicitly al  1's
    flag_a zflag;	//-z randomize
    flag_a Zflag;	//-Z use "in-house" SVM estimator
  };

  //selects out the nearest distances:
  template <class real>
  void kinearest(real **xvec,		//sample vectors 
		nel_ta n, 		//number of samples
		dim_ta D, 		//dimension of each sample
		real *xtest, 		//test point
		nel_ta k,		//number to select
		real *knearest,		//k-nearest
		long *ind,		//indices
		//metric to use:
  		real (*metric2) (real *, real *, dim_ta)=&metric2);
  template <class real>
  void knearest(real **xvec, nel_ta n, dim_ta D, real *xtest, nel_ta k, real *knearest,
  		real (*metric2) (real *, real *, dim_ta)=&metric2);


  //normalizes a set of vectors from averages and standard deviations:
  template <class real>
  void norm_vec(real **mat, 		//list of vectors
		dim_ta D, 		//dimension
		nel_ta n, 		//number of vectors
		real *ave, 		//returned averages (size D)
		real *std);		//returned std. dev. (size D)

  //un-normalizes a set of vectors from averages and standard deviations:
  template <class real>
  void unnorm_vec(real **mat, 		//vectors (test or training)
		dim_ta D, 		//dimension of vectors
		nel_ta n, 		//number of vectors
		real *ave, 		//averages
		real *std);		//standard deviations

  //parse command line options:
  int agf_parse_command_opts (int &argc, 	//number of arguments
		char **&argv, 			//arguments from command line
		const char * optlist, 		//list of options
		agf_command_opts *opt_args);	//returned options

  //sorts a set of vectors by their classes:
  //
  //returns a longword array of length ncl+1 indexing into the rearranged mat
  //	giving the starting point of each class
  template <class real, class cls_t>
  nel_ta * sort_classes(real **mat, 	//coordinate data
		nel_ta n,		//number of samples
		cls_t *cl, 		//classes
		cls_t ncl);		//number of classes

  //calculates averages and standard deviations of coordinate data:
  template <class real>
  void calc_norm(real **mat, 		//coordinate samples
		dim_ta D, 		//coordinate dimension
		nel_ta n, 		//number of samples
		real *ave, 		//D returned averages
		real *std);		//D returned standard deviations

  //normalizes both test and training coordinates with averages and std. dev.:
  template <class real>
  void norm_vec_std(real **mat1,	//training vectors 
		dim_ta D, 		//number of dimensions
		nel_ta n1, 		//number of training vectors
		real **mat2,		//test vectors 
		nel_ta n2);		//number of test vectors

  //randomly permutes training data:
  template <class real, class cls_t>
  void randomize_vec(real **mat,	//coordinates
		dim_ta D, 		//number of dimensions
		nel_ta n, 		//number of samples
		cls_t *cls);		//classes

  //removes duplicate vectors and returns a set of weights with the
  //number of duplicates:
  template <class real>
  real ** remove_duplicates(real **mat, //set of vectors
		dim_ta D, 		//dimension
		nel_ta n, 		//number of vectors
		real *wt, 		//returned weights
		nel_ta &nnew);		//new number of samples

  //make class labels go from 0-nc-1
  //including negative labels
  template <typename cls_t>
  cls_t compress_labels(cls_t *cls, 	//classes -in -out
		  nel_ta n);		//number of classes

}

#endif
