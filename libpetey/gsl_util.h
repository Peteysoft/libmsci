#ifndef GSL_UTIL__H
#define GSL_UTIL__H

#include <stdio.h>
#include <stdio.h>

#include <gsl/gsl_matrix.h>

namespace libpetey {

  //random debugging code...
  void print_gsl_matrix(FILE *fs, gsl_matrix * mat);
  void print_gsl_matrix_complex(FILE *fs, gsl_matrix_complex * mat);
  gsl_matrix *random_gsl_matrix(int m, int n);
  void random_gsl_matrix(gsl_matrix *a);

  //throws error using C++ throw operator:
  void gsl_throw_handler(const char *reason, const char *file, int line, int gsl_err);
  //prints out message but doesn't otherwise interrupt the program:
  void gsl_continue_handler(const char *reason, const char *file, int line, int gsl_err);

}

#endif

