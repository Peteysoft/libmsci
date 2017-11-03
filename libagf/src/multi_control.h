#ifndef MULTI_CONTROL__H
#define MULTI_CONTROL__H 1

#include <stdio.h>

namespace libagf {
  //calculate hausdorff metric between each pair of classes:
  template <class real, class cls_t>
  real * class_triangle(real **x, cls_t *cls, int n, int D);

  //generate the most common kinds of coding matrices:
  //(non-hierarchical component only)
  template <typename scalar>
  scalar ** one_against_all(int ncls);
  template <typename scalar>
  scalar ** one_against_one(int ncls);
  template <typename scalar>
  scalar ** partition_adjacent(int ncls);
  template <typename scalar>
  scalar ** random_coding_matrix(int ncls, int &ntrial, int strictflag=0);
  template <typename scalar>
  scalar ** exhaustive_coding_matrix(int ncls);
  template <typename scalar>
  scalar ** ortho_coding_matrix_nqbf(int ncls, 
		  int strictflag=1, 		//zeroes not allowed
		  int toprow1=0);		//top row is all ones
  template <typename scalar>
  scalar ** ortho_coding_matrix_brute_force(int ncls, int toprow1=0);
  template <typename scalar>
  scalar ** hierarchical_nonhierarchical(int ncls);

  //generate common control files:
  //(complete control file)
  void print_control_hier(FILE *fs, int ncls, int c0=0, int depth=0);
  template <typename scalar>
  void print_control_nonhier(FILE *fs, 		//output file stream
		  scalar **coding_matrix, 	//coding matrix
		  int n, 			//number of rows (models)
		  int ncls, 			//number of columns (classes)
		  const char *opt=NULL);

}

#endif

