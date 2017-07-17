//
// Not much of this shit is used any more. Kept for backwards compatability.
//

#ifndef CLASS_BORDERS_H_INCLUDED
#define CLASS_BORDERS_H_INCLUDED 1

#include "agf_defs.h"

//**note: find_class_borders now just wraps new, more general routines in
//sample_class_borders

namespace libagf {

  //constants and global variables:
  extern long agf_global_borders_maxiter;

  //parameter structure specifically for AGF borders training
  //to be used by agfrfunc
  template <class real>
  struct agfparam {
    //parameters for AGF classifiation:
    real W;			//total weight
    real var0[2];		//variance brackets (fixed)
    real var[2];		//variance brackets (floating)
    nel_ta k;			//number of nearest neighbours

    //diagnostic parameters:
    iter_ta min_nd;		//minimum number of iterations in weights calculation
    iter_ta max_nd;		//maximum		"
    iter_ta total_nd;		//total (average) 	"
    real min_f;			//minimum minimum weight
    real max_f;			//maximum minimum weight
    real total_f;		//total of minimum weights
    real min_W;
    real max_W;
    real total_W;
    iter_ta ncall;		//number of calls to agf_calc_w

  };

  //mainly to zero all those bloody diagnostic parameters:
  template <class real>
  int agfparam_init(agfparam<real> *param, real var[2], nel_ta k, real W,
		  int stickyflag=0);

  //stick the training data plus AGF parameters into the structure:
  template <class real>
  int agfbordparam_init(bordparam<real> *param,
		real **train, 		//training samples
		dim_ta D, 		//number of dimensions
		nel_ta n, 		//number of samples
		nel_ta clind, 		//index of first point in second class
		real var[2], 		//variance bracket
		nel_ta k, 		//number of nearest neighbours
		real W,			//total weights
		int smallflag=0);	//small dataset?

  //reset variance brackets if necessary:
  template <typename real>
  void agfparam_reset_var(agfparam<real> *param);

  template <class real>
  void agfbordparam_clean(bordparam<real> *param);

  //print out AGF diagnostic parameters:
  template <class real>
  void print_agfborddiag(bordparam<real> *param, FILE *fs);

  //returns the estimated difference in conditional probabilities:
  template <class real>
  real agfrfunc(real *x, 		//test point
		void *param, 		//parameters: training samples etc.
		real *drdx);		//gradient of difference in cond. prob.

  //for borders classification:
  //(returns number of samples actually found)
  template <class real>
  nel_ta find_class_borders(real **x,	//location of samples, sorted by class
		dim_ta nvar,		//number of dimensions
		nel_ta ntrain,		//number of training vectors
		nel_ta clind,		//index of start of second class
		nel_ta nsamp,		//number of samples to take
		real var[2],		//filter variance brackets
		nel_ta k,		//number of nearest neighbours to use
		real wc,		//"critical" weight
		real tol,		//desired tolerance
		real ** border,		//returned border vectors       
		real ** gradient,	//returned gradient vectors
		real rthresh=0);	//threshold for discrimination b

  //for small datasets:
  template <class real>
  nel_ta find_class_borders_small(real **x, dim_ta D, nel_ta ntrain, nel_ta clind,
			nel_ta &nsamp, real var[2], nel_ta k, real wc, real tol,
			real ** border, real ** gradient, real rthresh=0);

  //for "classifying" classes:
  //parses a partition specified on the command line:
  //returns number of new classes
  template <class cls_t>
  int parse_partition(int argc,		//number of arguments
		char **argv, 		//arguments
		cls_t ncls, 		//number of classes in original data
		cls_t *map);		//mapping from original class labels to new ones

  //applies a partition to a set of class labels:
  template <class cls_t>
  void apply_partition(cls_t *cls, 	//class labels
		nel_ta n, 		//number of samples
		cls_t *map);		//mapping from original to new class labels

}

#endif

