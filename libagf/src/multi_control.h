#ifndef MULTI_CONTROL__H
#define MULTI_CONTROL__H 1

#include <stdio.h>

namespace libagf {
  //calculate hausdorff metric between each pair of classes:
  template <class real, class cls_t>
  real * class_triangle(real **x, cls_t *cls, int n, int D);

  //given the class distance triangle, calculate the Hausdorff distance
  //between each side of a partition in a single row of the coding matrix:
  template <typename real, typename scalar>
  real partition_distance(real *d,		//distance triangle
		  int n,			//number of classes
		  scalar *row);			//coding row

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
  scalar ** orthogonal_coding_matrix(int ncls, int ncode);
  template <typename scalar>
  scalar ** exhaustive_coding_matrix(int ncls);
  template <typename scalar>
  scalar ** ortho_coding_matrix_nqbf(int ncls, 
		  int strictflag=1, 		//zeroes not allowed
		  int toprow1=0);		//top row is all ones
  template <typename scalar>
  scalar ** ortho_coding_matrix_greedy(int m, int n, int toprow1=0);
  template <typename scalar>
  scalar ** ortho_coding_matrix_brute_force(int ncls);
  template <typename scalar>
  scalar ** hierarchical_nonhierarchical(int ncls);

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

  //generate common control files:
  //(complete control file)
  void print_control_hier(FILE *fs, int ncls, int c0=0, int depth=0);

  template <typename scalar>
  void print_control_nonhier(FILE *fs, 		//output file stream
		  scalar **coding_matrix, 	//coding matrix
		  int n, 			//number of rows (models)
		  int ncls, 			//number of columns (classes)
		  char **opt=NULL,		//list of options/models
		  cls_ta *label=NULL);		//list of class labels

}

#endif

