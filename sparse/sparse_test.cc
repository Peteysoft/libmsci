#include <math.h>
#include <stdio.h>

#include "full_util.h"
#include "sparse.h"
#include "randomize.h"


//test cases:
//  y=A*b		vector multiplication
//  y=b*A
//  C=A*B		matrix multiplication
//  C=A*B		multiply sparse with full
//  C=B*A
//  C=A+B		matrix addition
//  C=A+B		add sparse to full
namespace libpetey {
  namespace libsparse {

    using namespace libpetey;
    using namespace libsparse;

    template <typename real>
    real test_sparse_arithmetic(
		int m,		//size of matrices
		int n,
		int p,
		real sparsity,	//sparsity of matrix
		real eps,
	 	int rantype) {	//uniform (0) or Gaussian (1)?

      double (*rangen) ();
      FILE *logfs=stderr;
      real res;		//residual
      real max_res=0;	//maximum size of residual
      int err=0;		//error code (number of bad values)
      real **A;
      real **B;
      real **C;
      real *v;
      real *y;

      if (rantype) {
        rangen=&rang;
      } else {
        rangen=&ranu;
      }

      //full matrices:
      fprintf(logfs, "Allocating matrices\n");
      A=allocate_matrix<real>(m, n);
      zero_matrix(A, m, n);
      B=allocate_matrix<real>(n, p);
      zero_matrix(B, n, p);
      C=allocate_matrix<real>(m, p);

      v=new real[n];
      y=new real[m];

      fprintf(logfs, "Initializing full matrices\n");

      for (int i=0; i<n; i++) v[i]=(*rangen)();

      for (int i=0; i<m; i++) {
        for (int j=0; j<n; j++) {
          if (ranu() > sparsity) A[i][j]=(*rangen)();
        }
      }
      for (int i=0; i<n; i++) {
        for (int j=0; j<p; j++) {
          if (ranu() > sparsity) B[i][j]=(*rangen)();
        }
      }

      //sparse matrices:
      fprintf(logfs, "Initializing sparse matrices\n");
      sparse<int32_t, real> A1(A, m, n);
      sparse<int32_t, real> B1(B, n, p);
      sparse<int32_t, real> C1(m, p);
      real **C2=allocate_matrix<real>(m, p);
      real y1[m];

      //beginning of test cases:

      //y=A*v
      fprintf(logfs, "y=A*v\n");
      vector_mult(A, v, y, m, n);		//full version
      A1.vect_mult(v, y1);		//sparse version

      for (int i=0; i<m; i++) {
        res=fabs(2*(y[i]-y1[i])/(y[i]+y1[i]));
        if (res > max_res) max_res=res;
	if (res > eps) {
          fprintf(stderr, "Residual of %g found; %g max\n", res, eps);
	  fprintf(stderr, "%g %g\n", y[i], y1[i]);
	  err++;
	}
      }
 
      fprintf(logfs, "y=v*A\n");
      left_vec_mult(v, A, y, m, n);	//full version
      A1.left_mult(v, y1);		//sparse version

      for (int i=0; i<m; i++) {
        res=fabs(2*(y[i]-y1[i])/(y[i]+y1[i]));
        if (res > max_res) max_res=res;
	if (res > eps) {
          fprintf(stderr, "Residual of %g found; %g max\n", res, eps);
	  fprintf(stderr, "%g %g\n", y[i], y1[i]);
	  err++;
	}
      }

      fprintf(logfs, "C=A*B\n");
      matrix_mult(A, B, C, m, n, p);	//full version
      A1.mat_mult(B1, C1);			//sparse version

      for (int i=0; i<m; i++) {
        for (int j=0; j<p; j++) {
          res=fabs(2*(C[i][j]-C1(i,j))/(C[i][j]+C1(i,j)));
          if (res > max_res) max_res=res;
	  if (res > eps) {
            fprintf(stderr, "Residual of %g found; %g max\n", res, eps);
	    fprintf(stderr, "%g %g\n", C[i][j], C1(i,j));
            err++;
          }
        }
      }

      fprintf(logfs, "C=A*B (sparse multiplicor)\n");
      A1.mat_mult(B, C2, p);		//sparse version

      for (int i=0; i<m; i++) {
        for (int j=0; j<p; j++) {
          res=fabs(2*(C[i][j]-C2[i][j])/(C[i][j]+C2[i][j]));
          if (res > max_res) max_res=res;
        }
      }

      fprintf(logfs, "C=A*B (sparse multiplicand)\n");
      B1.left_m_mult(A, C2, m);		//sparse version

      for (int i=0; i<m; i++) {
        for (int j=0; j<p; j++) {
          res=fabs(2*(C[i][j]-C2[i][j])/(C[i][j]+C2[i][j]));
          if (res > max_res) max_res=res;
	  if (res > eps) {
            fprintf(stderr, "Residual of %g found; %g max\n", res, eps);
	    fprintf(stderr, "%g %g\n", C[i][j], C2[i][j]);
            err++;
          }
        }
      }

      //matrix addition:
      fprintf(logfs, "Resizing matrics in preparation for addition\n");
      delete_matrix(B);
      delete_matrix(C);
      delete_matrix(C2);
      B=allocate_matrix<real>(m, n);
      B1.reset(m, n);
      for (int i=0; i<m; i++) {
        for (int j=0; j<n; j++) {
          if (ranu() > sparsity) {
            B[i][j]=(*rangen)();
            B1.add_el(B[i][j], i, j);
          }
        }
      }

      fprintf(logfs, "C=A+B (one sparse term\n");
      C=libpetey::copy_matrix(A, m, n);	//full version
      matrix_add(C, B, m, n);
      C2=libpetey::copy_matrix(B, m, n);	//sparse version
      A1.full_add(C2);

      for (int i=0; i<m; i++) {
        for (int j=0; j<n; j++) {
          res=fabs(2*(C[i][j]-C2[i][j])/(C[i][j]+C2[i][j]));
          if (res > max_res) max_res=res;
	  if (res > eps) {
            fprintf(stderr, "Residual of %g found; %g max\n", res, eps);
	    fprintf(stderr, "%g %g\n", C[i][j], C2[i][j]);
            err++;
          }
        }
      }

      fprintf(logfs, "C=A+B (two sparse terms)\n");
      C1=B1;
      C1.sparse_add(A1);			//sparse version

      for (int i=0; i<m; i++) {
        for (int j=0; j<n; j++) {
          res=fabs(2*(C[i][j]-C1(i, j))/(C[i][j]+C1(i, j)));
          if (res > max_res) max_res=res;
	  if (res > eps) {
            fprintf(stderr, "Residual of %g found; %g max\n", res, eps);
	    fprintf(stderr, "%g %g\n", C[i][j], C1(i,j));
            err++;
          }
        }
      }

      //clean up:
      fprintf(logfs, "Cleaning up\n");
      delete_matrix(A);
      delete_matrix(B);
      delete_matrix(C);
      delete_matrix(C2);
      delete [] v;
      delete [] y;

      return max_res;
    }

    template float test_sparse_arithmetic<float>(int, int, int, float, float, int);
    template double test_sparse_arithmetic<double>(int, int, int, double, double, int);
  }
}
