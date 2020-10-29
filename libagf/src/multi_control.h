#ifndef MULTI_CONTROL__H
#define MULTI_CONTROL__H 1

#include <stdio.h>

namespace libagf {
  template <class vector_t>
  void random_coding_row(vector_t &coding_row, int n, int strictflag);

  template <class vector_t>
  int check_coding_row(vector_t &coding_row, int n);

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
  scalar ** ortho_coding_matrix_greedy(int m, int n, int toprow1=0);
  template <typename scalar>
  scalar ** hierarchical_nonhierarchical(int ncls);


  template <typename scalar>
  scalar ** orthogonal_coding_matrix_with_degenerates(int ncls, int ncode, int strictflag=0);
  template <typename scalar>
  scalar ** orthogonal_coding_matrix(int ncls, int &ncode, int strictflag=0);

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

