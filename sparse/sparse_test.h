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

    //test cases:
    //  fill S1
    //  S2=S1:			assign S1 to S2 and compare
    //  S1.write(fs); S2.read(S2); write to binary file and read back again
    //  S1.print(fs); S2.scan(S2); write to ASCII file and read back again
    //  
    template <typename real>
    real test_sparse_mem(
	int m,		//size of matrices
	int n,
	real sparsity,	//sparsity of matrix
	real eps,	//acceptable relative error (less than precision of ASCII output)
	int rantype);	//uniform (0) or Gaussian (1) deviates?
  }
}

#endif

