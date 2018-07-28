#include <gsl/gsl_eigen.h>
#include <gsl/gsl_blas.h>
#include <math.h>
#include <gsl/gsl_linalg.h>

#include "gsl_util.h"
#include "randomize.h"

using namespace libpetey;

int main () {
  ran_init();
  int n=5;
  gsl_eigen_nonsymm_workspace *work=gsl_eigen_nonsymm_alloc(n);
  gsl_matrix *A=random_gsl_matrix(n, n);
  gsl_matrix *B=random_gsl_matrix(n, n);
  gsl_matrix *prod=gsl_matrix_alloc(n, n);
  gsl_vector_complex *eval1=gsl_vector_complex_alloc(n);
  gsl_vector_complex *eval2=gsl_vector_complex_alloc(n);

  //multiply two matrices:
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, A, B, 0, prod);

  //find the eigenvalues:
  gsl_eigen_nonsymm_params(0, 1, work);
  gsl_eigen_nonsymm(prod, eval1, work);

  //multiply two matrices:
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, B, A, 0, prod);

  //find the eigenvalues:
  gsl_eigen_nonsymm_params(0, 1, work);
  gsl_eigen_nonsymm(prod, eval2, work);

  printf("eigen value comparison:\n");
  for (int i=0; i<n; i++) {
    printf("%g + %g i   %g + %g i\n", 
		      GSL_REAL(gsl_vector_complex_get(eval1, i)),
		      GSL_IMAG(gsl_vector_complex_get(eval1, i)),
		      GSL_REAL(gsl_vector_complex_get(eval2, i)),
		      GSL_IMAG(gsl_vector_complex_get(eval2, i)));
  }
  printf("\n");

  //check same thing for singular values:
  gsl_matrix *v=gsl_matrix_alloc(n, n);
  gsl_vector *s1=gsl_vector_alloc(n);
  gsl_vector *s2=gsl_vector_alloc(n);
  gsl_vector *w2=gsl_vector_alloc(n);
  gsl_eigen_symm_workspace *w3=gsl_eigen_symm_alloc(n);
  gsl_matrix *p2=gsl_matrix_alloc(n, n);

  //multiply two matrices:
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, A, B, 0, prod);
  //gsl_linalg_SV_decomp_jacobi(prod, v, s1);
  gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1, prod, prod, 0, p2);
  //gsl_eigen_symm(p2, s1, w3);
  gsl_eigen_nonsymm(p2, eval1, work);

  //multiply two matrices:
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, B, A, 0, prod);
  //gsl_linalg_SV_decomp_jacobi(prod, v, s2);
  gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1, prod, prod, 0, p2);
  //gsl_eigen_symm(p2, s2, w3);
  gsl_eigen_nonsymm(p2, eval2, work);
/*
  printf("eigen value comparison:\n");
  for (int i=0; i<n; i++) {
    printf("%g  %g\n", gsl_vector_get(s1, i), gsl_vector_get(s2, i));
  }
*/
  printf("eigen value comparison:\n");
  for (int i=0; i<n; i++) {
    printf("%g + %g i   %g + %g i\n", 
		      GSL_REAL(gsl_vector_complex_get(eval1, i)),
		      GSL_IMAG(gsl_vector_complex_get(eval1, i)),
		      GSL_REAL(gsl_vector_complex_get(eval2, i)),
		      GSL_IMAG(gsl_vector_complex_get(eval2, i)));
  }
  printf("\n");

  random_gsl_matrix(A);

  gsl_eigen_symmv_workspace *w4=gsl_eigen_symmv_alloc(n);

  printf("Singular values over time:\n");
  for (int i=0; i<100; i++) {
    gsl_matrix *swp;
    random_gsl_matrix(B);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, B, A, 0, prod);
    gsl_matrix_memcpy(A, prod);
    gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1, prod, prod, 0, p2);
    gsl_linalg_SV_decomp(prod, B, s1, w2);
    //gsl_eigen_symmv(p2, s1, B, w4);
    //for (int j=0; j<n; j++) printf("%12.6g ", gsl_vector_get(s1, j));
    for (int j=0; j<n; j++) printf("%12.6g ", log(gsl_vector_get(s1, j))/i);
    printf("\n");
    //swp=A;
    //A=prod;
    //prod=A;
  }
  printf("\n");
  printf("final product:\n");
  print_gsl_matrix(stdout, A);

  printf("\n");
  printf("final eigenvectors:\n");
  print_gsl_matrix(stdout, B);

}

