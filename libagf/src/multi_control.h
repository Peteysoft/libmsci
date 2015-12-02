#ifndef MULTI_CONTROL__H
#define MULTI_CONTROL__H 1

#include <stdio.h>

namespace libagf {
  //calculate hausdorff metric between each pair of classes:
  template <class real, class cls_t>
  void class_triangle(real **x, cls_t *cls, int n, int D);

  //generate the most common kinds of coding matrices:
  //(non-hierarchical component only)
  int ** one_against_all(int ncls);
  int ** one_against_one(int ncls);
  int ** partition_adjacent(int ncls);
  int ** random_coding_matrix(int ncls, int &ntrial, int strictflag=0);
  int ** exhaustive_coding_matrix(int ncls);
  int ** ortho_coding_matrix_nqbf(int ncls, int strictflag=1);
  int ** ortho_coding_matrix_brute_force(int ncls);

  //generate common control files:
  //(complete control file)
  void print_control_hier(FILE *fs, int ncls, int c0=0, int depth=0);
  void print_control_nonhier(FILE *fs, int **coding_matrix, int n, int ncls, const char *opt=NULL);

}

#endif

