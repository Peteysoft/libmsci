#ifndef AGF_TYPES_H_INCLUDED
#define AGF_TYPES_H_INCLUDED 1

#include <stdint.h>

#include "randomize.h"

//command filenames:
#define AGF_OPT_VER ""
#define AGF_COMMAND_PREFIX ""
#define AGF_BINARY_CLASSIFIER "class_borders"
#define AGF_LTRAN_COM "agf_precondition"
#define AGF_PARTCOM "agf_preprocess"

//names of important functions
#define KLEAST_FUNC kleast_quick	//function for finding nearest neighbours
#define AGF_CALC_W_FUNC agf_calc_w2	//function for calculating AGF weights

namespace libagf {

  typedef int32_t cls_ta;		//class data
  typedef int32_t dim_ta;		//dimension data
  typedef int32_t nel_ta;		//sample sizes
  typedef int32_t iter_ta;		//number of iterations
  typedef float real_a;			//floating point data
  typedef int flag_a;			//boolean data
  typedef int enum_a;			//enumerated data

  //defaults and constants:
  //global default parameters:
  const iter_ta BORDERS_MAXITER=100;	//maximum iterations in borders calculation
  const iter_ta WEIGHTS_MAXITER=100;	//       "	in AGF weights calculation
  const real_a WEIGHTS_TOL=0.001;	//tolerance of AGF weights calculation

  //global defaults:
  //AGF and KNN:
  const real_a W_DEFAULT=11;		//W parameter in standard AGF
  const nel_ta K_DEFAULT_KNN=11;	//KNN nearest neighbours
  const nel_ta K_DEFAULT_AGF=-1;	//AGF nearest neighbours (-1=all data)
  const real_a F_DEFAULT=0.1;		//default test proportion

  //for agf borders:
  const nel_ta NBORD_DEFAULT=100;	//number of border samples
  const real_a TOL_DEFAULT=0.00001;	//tolerance of border samples
  const real_a W_DEFAULT_BORDERS=20;	//W parameter for AGF borders
  const nel_ta SMALL_DATASET=200;	//threshold find_class_borders vs. find_class_borders_small
  const real_a RTHRESH_DEFAULT=0;	//location of Bayesian border
  const real_a HREL_DEFAULT=0.0001;	//relative difference for numerical derivatives.

  //for "optimal" AGF:
  const iter_ta NT_DEFAULT=10;		//number of trials
  const real_a WMIN_DEFAULT=2.;		//minimum weights total
  const real_a WMAX_DEFAULT=100.;	//maximum weights total

  //number of histogram bins when comparing accuracy with confidence:
  const int32_t NCONHIST=10;

  //partition symbol for multi-class module:
  const char PARTITION_SYMBOL='/';

  //returned diagnostics for AGF:
  struct agf_diag_param {
    iter_ta nd;           //number of iterations in weights calculation
    real_a f;             //ratio of min. weight to max.
    real_a W;             //total weight
  };

  //simple routines for class selection
  //all of them return random results in the case of a tie to reduce bias
  template <class cls_t, class number>
  inline cls_t choose_class(number *p, cls_t n) {
    cls_t cls[n];
    cls_t ntie=1;
    cls_t result;
    cls[0]=0;
    for (cls_t i=1; i<n; i++) {
      if (p[i]>p[cls[0]]) {
        cls[0]=i;
        ntie=1;
      } else if (p[i]==p[cls[0]]) {
        cls[ntie]=i;
        ntie++;
      }
    }
    if (ntie==1) result=cls[0]; else result=cls[(cls_t) libpetey::ranu()*ntie];
    return result;
  }

  template <class real>
  inline int convertR(real r) {
    if (r<0) {
      return 0;
    } else if (r>0) {
      return 1;
    } else {
      return libpetey::ranu()*2;
    }
  }

  template <class real>
  inline int convertR(real r, real *p) {
    int result=convertR(r);
    p[0]=(1-r)/2;
    p[1]=1-p[0];
    return result;
  }

  template <class real>
  inline int convertR(real r, real &p) {
    int result;
    if (r<0) {
      result=0;
      p=(1-r)/2;
    } else if (r>0) {
      result=1;
      p=(1+r)/2;
    } else {
      result=libpetey::ranu()*2;
      p=0;
    }
    return result;
  }
}

#endif
