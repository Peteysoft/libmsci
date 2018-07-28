#include <assert.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_complex_math.h>

#include "peteys_no_templates.h"
#include "gsl_util.h"
#include "randomize.h"
#include "solve_lode.h"

namespace libpetey {

  int compare_gsl_complex(gsl_complex *c1, gsl_complex *c2) {
    if (GSL_REAL(*c1) > GSL_REAL(*c2)) return 1;
    if (GSL_REAL(*c1) < GSL_REAL(*c2)) return -1;
    if (GSL_IMAG(*c1) > GSL_IMAG(*c2)) return 1;
    if (GSL_IMAG(*c1) < GSL_IMAG(*c2)) return -1;
    return 0;
  }

  void check_eig_decomp(gsl_matrix_complex *v, gsl_vector_complex *eval,
		  gsl_matrix_complex *vinv, gsl_matrix_complex *result) {
    int n=eval->size;
    int err;
    gsl_matrix_complex *emat=gsl_matrix_complex_alloc(n, n);

    gsl_matrix_complex_set_zero(emat); 
    for (int i=0; i<n; i++) gsl_matrix_complex_set(emat, i, i, 
		gsl_vector_complex_get(eval, i));
    err=gsl_blas_zgemm (CblasNoTrans, CblasNoTrans, 
		gsl_complex_rect(1, 0), emat, vinv,
                gsl_complex_rect(0, 0), result);
    err=gsl_blas_zgemm (CblasNoTrans, CblasNoTrans, 
		gsl_complex_rect(1, 0), v, result, 
                gsl_complex_rect(0, 0), emat);
    gsl_matrix_complex_memcpy(result, emat);

    gsl_matrix_complex_free(emat);
  }

  int test_lode(int n,				//size of test case
			double tol) {		//expected tolerance
    gsl_matrix_complex *v;
    gsl_matrix_complex *vinv;
    gsl_vector_complex *eval;
    gsl_matrix_complex *result;
    int err=0;

    //random matrix:
    gsl_matrix *A=random_gsl_matrix(n, n);

    //allocate matrices:
    v=gsl_matrix_complex_alloc(n, n);
    vinv=gsl_matrix_complex_alloc(n, n);
    result=gsl_matrix_complex_alloc(n, n);
    eval=gsl_vector_complex_alloc(n);

    diagonalize_matrix(A, v, eval, vinv);

    err=gsl_blas_zgemm (CblasNoTrans, CblasNoTrans, 
		gsl_complex_rect(1, 0), vinv, v,
                gsl_complex_rect(0, 0), result);
    //check each of the values:
    for (int i=0; i<n; i++) {
      for (int j=0; j<n; j++) {
        if (i==j) {
          if (gsl_complex_abs(gsl_complex_sub_real(gsl_matrix_complex_get(result, i, j), 1))>tol) {
            fprintf(stderr, "test_lode: first test, V^-1*V, failed, element (%d, %d)\n", i, j);
	    err=1;
	  } 
	} else {
          if (gsl_complex_abs(gsl_matrix_complex_get(result, i, j))>tol) {
            fprintf(stderr, "test_lode: first test, V^-1*V, failed, element (%d, %d)\n", i, j);
	    err=1;
	  }
	}
      }
    }
    if (err!=0) {
      printf("Should be identity:\n");
      print_gsl_matrix_complex(stdout, result);
    }

    err=gsl_blas_zgemm (CblasNoTrans, CblasNoTrans, 
		gsl_complex_rect(1, 0), v, vinv,
                gsl_complex_rect(0, 0), result);
    for (int i=0; i<n; i++) {
      for (int j=0; j<n; j++) {
        if (i==j) {
          if (gsl_complex_abs(gsl_complex_sub_real(gsl_matrix_complex_get(result, i, j), 1))>tol) {
            fprintf(stderr, "test_lode: second test, V*V^-1, failed, element (%d, %d)\n", i, j);
	    err=1;
	  } 
	} else {
          if (gsl_complex_abs(gsl_matrix_complex_get(result, i, j))>tol) {
            fprintf(stderr, "test_lode: second test, V*V^-1, failed, element (%d, %d)\n", i, j);
	    err=1;
	  }
	}
      }
    }

    if (err!=0) {
      printf("Should be identity:\n");
      print_gsl_matrix_complex(stdout, result);
    }

    check_eig_decomp(v, eval, vinv, result);

    for (int i=0; i<n; i++) {
      for (int j=0; j<n; j++) {
        gsl_complex val1=gsl_matrix_complex_get(result, i, j);
	double val2=gsl_matrix_get(A, i, j);
        if (gsl_complex_abs(gsl_complex_div(gsl_complex_sub_real(val1, val2), gsl_complex_add_real(val1, val2)))>tol/2) {
          fprintf(stderr, "test_lode: third test, V*lambda*V^-1==A, failed, element (%d, %d)\n", i, j);
          err=1;
        } 
      }
    }

    if (err!=0) {
      printf("The following matrices should be the same to round-off error:\n");
      print_gsl_matrix(stdout, A);
      printf("\n");
      print_gsl_matrix_complex(stdout, result);
    }

    gsl_matrix_complex_free(result);
    gsl_matrix_complex_free(vinv);
    gsl_matrix_complex_free(v);
    gsl_matrix_free(A);
    gsl_vector_complex_free(eval);
  }

  //performs the eigenvalue decomposition: A=v*eval*vinv
  //where vinv=v^(-1)
  int diagonalize_matrix(gsl_matrix *A0, gsl_matrix_complex *v, 
		gsl_vector_complex *eval, gsl_matrix_complex *vinv) {
    int err;
    int n;
    gsl_eigen_nonsymmv_workspace *work;
    gsl_matrix *A;
    gsl_matrix *AT;			//transposed matrix
    gsl_vector_complex *eval2;		//left eigen-values (should be same as right)
    gsl_vector_complex *eunsrt;		//unsorted eigen-values
    gsl_matrix_complex *vunsrt;		//unsorted eigen-vectors
    gsl_complex *esrt[A0->size1];	//for sorting eigen values
    long sind[A0->size1];		//sorting indices
    gsl_matrix_complex *id;		//for scaling eigenvectors

    assert(A0->size1==A0->size2);
    n=A0->size1;

    work=gsl_eigen_nonsymmv_alloc(n);
    gsl_eigen_nonsymmv_params(1, work);

    eunsrt=gsl_vector_complex_alloc(n);
    vunsrt=gsl_matrix_complex_alloc(n, n);
    eval2=gsl_vector_complex_alloc(n);
    id=gsl_matrix_complex_alloc(n, n);

    A=gsl_matrix_alloc(n, n);
    AT=gsl_matrix_alloc(n, n);
    err=gsl_matrix_memcpy(A, A0);
    if (err!=0) goto finish;
    err=gsl_matrix_transpose_memcpy(AT, A0);
    if (err!=0) goto finish;

    //get right eigenvectors:
    err=gsl_eigen_nonsymmv(A, eunsrt, vunsrt, work);
    if (err!=0) goto finish;
    //sort the eigenvalues alongside the eigenvectors:
    for (int i=0; i<n; i++) esrt[i]=gsl_vector_complex_ptr(eunsrt, i);
    heapsort((void **) esrt, sind, n, (void *) &compare_gsl_complex);
    for (int i=0; i<n; i++) {
      gsl_vector_complex_set(eval, i, gsl_vector_complex_get(eunsrt, sind[i]));
    }
    for (int i=0; i<n; i++) {
      for (int j=0; j<n; j++) {
        gsl_matrix_complex_set(v, i, j, 
			gsl_matrix_complex_get(vunsrt, i, sind[j]));
      }
    }
    //printf("Right eigenvectors, unsorted:\n");
    //print_gsl_matrix_complex(stdout, vunsrt);
    //printf("Right eigenvectors, sorted:\n");
    //print_gsl_matrix_complex(stdout, v);

    //get left eigenvectors:
    err=gsl_eigen_nonsymmv(AT, eunsrt, vunsrt, work);
    if (err!=0) goto finish;
    //sort these fuckers too:
    for (int i=0; i<n; i++) esrt[i]=gsl_vector_complex_ptr(eunsrt, i);
    heapsort((void **) esrt, sind, n, (void *) &compare_gsl_complex);
    for (int i=0; i<n; i++) {
      gsl_vector_complex_set(eval2, i, gsl_vector_complex_get(eunsrt, sind[i]));
    }
    for (int i=0; i<n; i++) {
      for (int j=0; j<n; j++) {
        gsl_matrix_complex_set(vinv, j, i, 
			gsl_matrix_complex_get(vunsrt, i, sind[j]));
      }
    }

    //printf("Left eigenvectors, unsorted, un-transposed:\n");
    //print_gsl_matrix_complex(stdout, vunsrt);
    //printf("Left eigenvectors, sorted, transposed:\n");
    //print_gsl_matrix_complex(stdout, vinv);

    //gsl_matrix_complex_transpose(vinv);
    //*** should do something with second set of eigenvalues...
    
    printf("eigen values:\n");
    printf("left         right  \n");
    for (int i=0; i<n; i++) {
      printf("%g + %g i   %g + %g i\n", 
		      GSL_REAL(gsl_vector_complex_get(eval, i)),
		      GSL_IMAG(gsl_vector_complex_get(eval, i)),
		      GSL_REAL(gsl_vector_complex_get(eval2, i)),
		      GSL_IMAG(gsl_vector_complex_get(eval2, i)));
    }
    printf("\n");

    //scale the eigenvectors so that left and right are the inverse of one another:
    err=gsl_blas_zgemm (CblasNoTrans, CblasNoTrans, 
		gsl_complex_rect(1, 0), vinv, v,
                gsl_complex_rect(0, 0), id);

    for (int i=0; i<n; i++) {
      gsl_complex norm=gsl_complex_sqrt(gsl_matrix_complex_get(id, i, i));
      for (int j=0; j<n; j++) {
        gsl_complex val=gsl_matrix_complex_get(v, j, i);
	gsl_matrix_complex_set(v, j, i, gsl_complex_div(val, norm));
        val=gsl_matrix_complex_get(vinv, i, j);
	gsl_matrix_complex_set(vinv, i, j, gsl_complex_div(val, norm));
      }
    }

    finish:
      gsl_matrix_free(A);
      gsl_matrix_free(AT);
      gsl_eigen_nonsymmv_free(work);
      gsl_vector_complex_free(eval2);
      gsl_vector_complex_free(eunsrt);
      gsl_matrix_complex_free(vunsrt);
      gsl_matrix_complex_free(id);

    return err;
  }

  //computes: x=vinv*exp(t*eval)*v*x0
  int solve_lode(gsl_vector *x0, gsl_matrix_complex *v, 
		  gsl_vector_complex *eval,  gsl_matrix_complex *vinv,
		  double t, gsl_vector *x) {
    gsl_vector_complex *x0c;		//initial conditions as complex number
    gsl_vector_complex *x0t;		//transformed initial conditions
    gsl_vector_complex *xf;		//final result as complex vector
    gsl_vector_complex *b;		//intermediate result
    int err;

    x0c=gsl_vector_complex_alloc(eval->size);
    x0t=gsl_vector_complex_alloc(eval->size);
    xf=gsl_vector_complex_alloc(eval->size);
    b=gsl_vector_complex_alloc(eval->size);

    for (int i=0; i<eval->size; i++) gsl_vector_complex_set(x0c, i, 
		    gsl_complex_rect(gsl_vector_get(x0, i), 0));

    //use the blas routines:
    err=gsl_blas_zgemv (CblasNoTrans, gsl_complex_rect(1, 0), vinv, x0c, 
                     gsl_complex_rect(0, 0), x0t);
    if (err!=0) goto finish;
    gsl_vector_complex_memcpy(b, eval);
    gsl_blas_zscal(gsl_complex_rect(t, 0), b);
    for (int i=0; i<eval->size; i++) {
      gsl_vector_complex_set(b, i, 
		      gsl_complex_exp(gsl_vector_complex_get(b, i)));
    }
    err=gsl_vector_complex_mul(b, x0t);
    if (err!=0) goto finish;
    err=gsl_blas_zgemv (CblasNoTrans, gsl_complex_rect(1, 0), v, b, 
                     gsl_complex_rect(0, 0), xf);
    if (err!=0) goto finish;

    //convert to real (should check imaginary part which should be small):
    for (int i=0; i<eval->size; i++) {
      gsl_vector_set(x, i, GSL_REAL(gsl_vector_complex_get(xf, i)));
    }
    
    printf("final result:\n");
    for (int i=0; i<eval->size; i++) {
      printf("%g + %g i\n", 
		      GSL_REAL(gsl_vector_complex_get(xf, i)),
		      GSL_IMAG(gsl_vector_complex_get(xf, i)));
    }

    finish:
      gsl_vector_complex_free(x0t);
      gsl_vector_complex_free(xf);
      gsl_vector_complex_free(b);

    return err;
  }

} //end namespace libpetey

