#ifndef MULTI_CONTROL_EXPERIMENTAL__H
#define MULTI_CONTROL_EXPERIMENTAL__H 1

#include <stdio.h>

namespace libagf {

  template <typename scalar>
  scalar ** ortho_coding_matrix_nqbf(int ncls, 
		  int strictflag=1, 		//zeroes not allowed
		  int toprow1=0);		//top row is all ones

  //find partitions with maximum distance between each side:
  template <typename scalar, typename real>
  scalar ** optimal_coding_matrix(int n,	//number of classes
		  int nrow,			//number of rows
		  real *d);			//distance triangle btw. classes

  //find partitions with maximum distance between each side:
  template <typename scalar, typename real, typename cls_t>
  scalar ** optimal_coding_matrix(real *d,	//distance triangle (all pts.)
		  cls_t *cls,			//classes
		  nel_ta n,			//number of points
		  int nrow);			//number of rows

}

#endif

