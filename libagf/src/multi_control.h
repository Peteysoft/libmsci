#ifndef MULTI_CONTROL__H
#define MULTI_CONTROL__H 1

#include <stdio.h>

namespace libagf {
  //write the most common kinds of control files:
  //(non-hierarchical component only)
  void one_against_all(FILE *fs, int ncls, const char *options=NULL);
  void one_against_one(FILE *fs, int ncls, const char *options=NULL);
  void partition_adjacent(FILE *fs, int ncls, const char *options=NULL);
  void random_coding_matrix(FILE *fs, int ncls, int ntrial, int strictflag=0);
  void exhaustive_coding_matrix(FILE *fs, int ncls);
  void ortho_coding_matrix_brute_force(FILE *fs, int ncls);

  //generate common control files:
  //(complete control file)
  void print_control_hier(FILE *fs, int ncls, int c0=0, int depth=0);
  void print_control_1vsall(FILE *fs, int ncls, const char *opt=NULL);
  void print_control_1vs1(FILE *fs, int ncls, const char *opt=NULL);
  void print_control_adj(FILE *fs, int ncls, const char *opt=NULL);
  void print_control_random(FILE *fs, int ncls, int nrow, int strictflag=0);
  void print_control_exhaustive(FILE *fs, int ncls);
  void print_control_ortho(FILE *fs, int ncls);

}

#endif

