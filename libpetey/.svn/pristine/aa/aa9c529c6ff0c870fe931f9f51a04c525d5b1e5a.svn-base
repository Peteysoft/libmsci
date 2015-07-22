#ifndef FULL_UTIL_INCLUDED
#define FULL_UTIL_INCLUDED 1

#include <stdio.h>

namespace libpetey {

template <class real, class integer>
real ** allocate_matrix(integer m, integer n);

template <class real, class integer>
void zero_matrix(real ** mat, integer m, integer n);

template <class real, class integer>
real ** zero_matrix(integer m, integer n);

template <class real, class integer>
void identity_matrix(real ** mat, integer m, integer n);

template <class real, class integer>
real ** identity_matrix(integer m, integer n);

template <class real>
void delete_matrix(real ** mat);

template <class real, class integer>
void copy_matrix(real **m1, real **m2, integer m, integer n);

template <class real, class integer>
real ** copy_matrix(real **mat, integer m, integer n);

//should be more efficient:
template <class real, class integer>
void matrix_mult_t(real **plier, real **cand, real **result, integer m, integer p, integer n);

template <class real, class integer>
void matrix_mult(real **plier, real **cand, real **result, integer m, integer p, integer n);

template <class real, class integer>
real ** matrix_mult(real **plier, real **cand, integer m, integer p, integer n);

template <class real, class integer>
void vector_mult(real **plier, real *cand, real *result, integer m, integer n);

template <class real, class integer>
real * vector_mult(real **plier, real *cand, integer m, integer n);

template <class real, class integer>
void left_vec_mult(real *plier, real **cand, real *result, integer m, integer n);

template <class real, class integer>
real * left_vec_mult(real *plier, real **cand, integer m, integer n);

template <class real, class integer>
void matrix_add(real **mat1, real ** mat2, integer m, integer n);

//for square matrices (inplace):
template <class real, class integer>
void matrix_transpose(real **mat, integer m);

template <class real, class integer>
real ** matrix_transpose(real **mat, integer m, integer n);

template <class real, class integer>
real ** scan_matrix(FILE *fptr, integer &m, integer &n, int flag=0);

template <class integer>
void print_matrix(FILE *fptr, float **mat, integer m, integer n);
template <class integer>
void print_matrix(FILE *fptr, double **mat, integer m, integer n);

template <class real, class integer>
real ** read_matrix(FILE *fptr, integer &m, integer &n);

template <class real, class integer>
size_t write_matrix(FILE *fptr, real **mat, integer m, integer n);

template <class real, class integer>
real matrix_norm(real **mat, integer m, integer n);

} //end namespace libpetey

#endif
