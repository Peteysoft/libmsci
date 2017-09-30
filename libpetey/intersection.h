#ifndef LIBPETEY_INTERSECTION_H
#define LIBPETEY_INTERSECTION_H

namespace libpetey {

  template <typename scalar>
  int compare_vectors(scalar *v1, scalar *v2, int n);

  template <typename scalar>
  int unify_vectors(scalar ***v, 		//list of lists of vectors
		  int *n, 			//number of elements in each list
		  int n0, 			//number of lists
		  int D,			//dimension of each vector
		  scalar **result,		//intersection
		  int **ind);			//indexes into result

  template <typename tp>
  int32_t intersection(tp *l1, int32_t n1, tp *l2, int32_t n2, tp *lr);

  //returns 0 for success:
  template <typename scalar>
  int test_unify_vectors(int nv, 		//total number of vectors(union)
		  int D, 			//dimension of each vector
		  int nss, 			//number of sub sets
		  float frac);			//fraction for each ss

}

#endif

