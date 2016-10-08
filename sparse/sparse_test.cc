#include <math.h>
#include <stdio.h>

#include "full_util.h"
#include "sparse.h"
#include "randomize.h"

using namespace libpetey;
using namespace libsparse;

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
      real res;			//residual
      real max_res=0;		//largest residual
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
/*
      print_matrix(stdout, A, m, n);
      printf("\n");
      print_matrix(stdout, B, n, p);
      printf("\n");
      A1.print(stdout);
      printf("\n");
      B1.print(stdout);
      printf("\n");
*/

      //beginning of test cases:

      //y=A*v
      fprintf(logfs, "y=A*v\n");
      vector_mult(A, v, y, m, n);		//full version
      A1.vect_mult(v, y1);			//sparse version
/*
      for (int i=0; i<m; i++) printf("%g ", y[i]);
      printf("\n\n");
      for (int i=0; i<m; i++) printf("%g ", y1[i]);
      printf("\n\n");
*/
      for (int i=0; i<m; i++) {
        if (y[i]==0 && y1[i]==0) res=0; 
		else res=fabs(2*(y[i]-y1[i])/(y[i]+y1[i]));
        if (res > max_res) max_res=res;
	if (res > eps) {
          fprintf(stderr, "Residual of %g found; %g max\n", res, eps);
	  fprintf(stderr, "%g %g\n", y[i], y1[i]);
	  err++;
	}
      }

      for (int i=0; i<m; i++) y[i]=(*rangen) ();
      fprintf(logfs, "v=y*A\n");
      left_vec_mult(y, A, v, m, n);	//full version
      real v1[m];
      A1.left_mult(y, v1);		//sparse version
/*
      for (int i=0; i<n; i++) printf("%g ", v[i]);
      printf("\n\n");
      for (int i=0; i<n; i++) printf("%g ", v1[i]);
      printf("\n\n");
*/
      for (int i=0; i<n; i++) {
        if (v[i]==0 && v1[i]==0) res=0; 
        	else res=fabs(2*(v[i]-v1[i])/(v[i]+v1[i]));
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
/*
      print_matrix(stdout, C, m, p);
      printf("\n");
      C1.print(stdout);
      printf("\n");
*/
      for (int i=0; i<m; i++) {
        for (int j=0; j<p; j++) {
          if (C[i][j]==0 && C1(i,j)==0) res=0; 
          	else res=fabs(2*(C[i][j]-C1(i,j))/(C[i][j]+C1(i,j)));
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
          if (C[i][j]==0 && C2[i][j]==0) res=0; 
          	else res=fabs(2*(C[i][j]-C2[i][j])/(C[i][j]+C2[i][j]));
          if (res > max_res) max_res=res;
        }
      }

      fprintf(logfs, "C=A*B (sparse multiplicand)\n");
      B1.left_m_mult(A, C2, m);		//sparse version

      for (int i=0; i<m; i++) {
        for (int j=0; j<p; j++) {
          if (C[i][j]==0 && C2[i][j]==0) res=0; 
          	else res=fabs(2*(C[i][j]-C2[i][j])/(C[i][j]+C2[i][j]));
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
      zero_matrix(B, m, n);
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
      //print_matrix(stdout, C, m, n);
      //print_matrix(stdout, B, m, n);
      matrix_add(C, B, m, n);
      C2=libpetey::copy_matrix(B, m, n);	//sparse version
      A1.full_add(C2);

      //print_matrix(stdout, C, m, n);
      //print_matrix(stdout, C2, m, n);

      for (int i=0; i<m; i++) {
        for (int j=0; j<n; j++) {
          if (C[i][j]==0 && C2[i][j]==0) res=0; 
          	else res=fabs(2*(C[i][j]-C2[i][j])/(C[i][j]+C2[i][j]));
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

      //print_matrix(stdout, C, m, n);
      //C1.print(stdout);

      for (int i=0; i<m; i++) {
        for (int j=0; j<n; j++) {
          if (C[i][j]==0 && C1(i,j)==0) res=0; 
          	else res=fabs(2*(C[i][j]-C1(i, j))/(C[i][j]+C1(i, j)));
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

    template <typename real>
    real test_sparse_mem(
		int m,		//size of matrices
		int n,
		real sparsity,
		real tol,
		int gdev) {
      double (*rangen) ();
      FILE *logfs=stderr;
      FILE *fs;
      sparse<int32_t, real> S1(m, n);
      sparse<int32_t, real> S2;
      int i, j;
      int errcount=0;
      real res;

      if (gdev) {
        rangen=&rang;
      } else {
        rangen=&ranu;
      }

      for (int k=0; k<n*m*(1-sparsity); k++) {
        i=ranu()*m;
	j=ranu()*n;
	S1.add_el((*rangen)(), i, j);
      }

      fprintf(logfs, "Assignment: S2=S1\n");
      S2=S1;

      /* there is no comparison operator!:
      if (S1==S2) {
        fprintf(logfs, "S1==S2: passed\n");
      } else {
        fprintf(stderr, "S1==S2: failed\n");
	errcount++;
      }
      */

      for (int i=0; i<n; i++) {
        for (int j=0; j<m; j++) {
          if (S1(i,j)!=S2(i,j)) {
            fprintf(stderr, "S1(%d,%d)==S2(%d,%d): failed\n", i, j, i, j);
	    errcount++;
	  }
	}
      }

      fprintf(logfs, "Binary write: S1.write(fs); S2.read(fs)\n");

      S1.reset(m, n);
      for (int k=0; k<n*m*(1-sparsity); k++) {
        i=ranu()*m;
	j=ranu()*n;
	S1.add_el((*rangen)(), i, j);
      }

      fs=tmpfile();
      S1.write(fs);
      rewind(fs);
      S2.read(fs);
      fclose(fs);

      for (int i=0; i<n; i++) {
        for (int j=0; j<m; j++) {
          if (S1(i,j)!=S2(i,j)) {
            fprintf(stderr, "S1(%d,%d)==S2(%d,%d): failed\n", i, j, i, j);
	    errcount++;
	  }
	}
      }

      S1.reset();
      for (int k=0; k<n*m*(1-sparsity); k++) {
        i=ranu()*m;
	j=ranu()*n;
	S1.add_el((*rangen)(), i, j);
      }

      fprintf(logfs, "ASCII write: S1.print(fs); S2.scan(fs)\n");
      fs=tmpfile();
      S1.print(fs);
      rewind(fs);
      S2.scan(fs);

      fclose(fs);

      for (int i=0; i<n; i++) {
        for (int j=0; j<m; j++) {
          res=fabs(2*(S1(i, j)-S2(i, j))/(S1(i, j)+S2(i, j)));
          if (res>tol) {
            fprintf(stderr, "S1(%d,%d)~=S2(%d,%d): failed\n", i, j, i, j);
	    errcount++;
	  }
	}
      }

      return errcount;
    }


    template float test_sparse_arithmetic<float>(int, int, int, float, float, int);
    template double test_sparse_arithmetic<double>(int, int, int, double, double, int);

    template float test_sparse_mem<float>(int, int, float, float, int);
    template double test_sparse_mem<double>(int, int, double, double, int);
  }
}
