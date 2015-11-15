#ifndef MULTI_CONTROL__H
#define MULTI_CONTROL__H 1

#include <stdio.h>

namespace libagf {
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

