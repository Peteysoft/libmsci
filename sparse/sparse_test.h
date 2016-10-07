#ifndef __SPARSE_TEST_H
#define __SPARSE_TEST_H

namespace libpetey {
  namespace libsparse {

    //test cases:
    //  y=A*b		vector multiplication
    //  y=b*A
    //  C=A*B		matrix multiplication
    //  C=A*B		multiply sparse with full
    //  C=B*A
    //  C=A+B		matrix addition
    //  C=A+B		add sparse to full

    template <typename real>
    real test_sparse_arithmetic(
	int m,		//size of matrices
	int n,
	int p,
	real sparsity,	//sparsity of matrix
	real eps,	//acceptable relative error (close to machine size)
	int rantype);	//uniform (0) or Gaussian (1) deviates?
  }
}

#endif

