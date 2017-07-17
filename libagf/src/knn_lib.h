//
// Library for k-nearest neighbours: classification, interpolation and density
// estimation.
//

#ifndef KNN_H_INCLUDED
#define KNN_H_INCLUDED 1

#include "agf_defs.h"

namespace libagf {

  //K-nearest neighbours (KNN) classifier:
  template <class real, class cls_t>
  cls_t knn(real (* metric) (real *, real *, dim_ta), 	//distance function
		real **mat, 		//matrix of coordinate data
		dim_ta D, 		//dimension of coordinate data
		nel_ta n, 		//number of training samples
		cls_t *cl, 		//class numbers [0, ncls-1]
		cls_t ncl, 		//number of classes
		real *vec,		//test vector
		nel_ta k, 		//number of nearest neighbours
		real *pdf, 		//returns conditional/joint prob.
		int joint=0);		//return joint probabilities?

  //estimate probability densities using KNN:
  template <class real>
  real knn_pdf(real **mat, 		//coordinate samples
		dim_ta D, 		//dimension of samples
		nel_ta n, 		//number of samples
		real *vec, 		//test point
		nel_ta k);		//number of nearest neighbours

  //perform kernel interpolation with KNN:
  template <class real>
  real int_knn(real (* metric) (real *, real *, dim_ta), //distance metric
		real **mat, 		//coordinate locations
		dim_ta D, 		//coordinate dimension
		nel_ta n, 		//number of samples
		real *ord, 		//ordinates
		real *vec, 		//test point
		nel_ta k);		//number of nearest neighbours
 
  //perform kernel interpolation with KNN--return RMS error estimates:
  template <class real>
  real int_knn(real (* metric) (real *, real *, dim_ta), //distance metric 
		real **mat, 
		dim_ta D, 
		nel_ta n, 
		real *ord,
		real *vec, 
		nel_ta k, 
		real &rms);		//error estimate from root-mean-square
}

#endif

