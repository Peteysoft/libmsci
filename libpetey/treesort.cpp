//really a tree-sort, but fills the role for the time being...

#include "tree_tmp.h"

namespace libpetey {

  //long * heapsort(float * data, long n);

  template <class dt>
  long * treesort(dt *data, long n) {
	tree<dt> data_tree;
	long *indices;
	long i;
	indices=new long[n];
	
	for (i=0; i<n; i++) data_tree.add(data[i], i);

	i=0;
	data_tree.decompose(data, indices, n, i);

	//for (i=0; i<n; i++) cout << indices[i] << " ";
	//cout << "\n";

	return indices;

  } 
}

