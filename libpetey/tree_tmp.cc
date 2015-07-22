#include <stdint.h>

#include "tree_tmp.h"
#include "string_petey.h"
#include "time_class.h"
#define TEST 0

#if TEST
  #include <iostream>
#endif

namespace libpetey {

template <class primitive>
tree<primitive>::tree() {
  value=NULL;
  left=NULL;
  right=NULL;
}

//adds an element to the tree
//If the element is a duplicate, it is not added
//and the function returns the index value of the duplicate tree element
//Otherwise the method returns a copy of the index
template <class primitive>
long tree<primitive>::add_member(primitive new_val, long new_ind) {
  if (value==NULL) {
    value=new primitive;
    *value=new_val;
    index=new_ind;
    //cout<<*value<<"\n";
  } else {
    if (new_val == *value) {
      return index;
    } else if (new_val < *value) {
      if (left==NULL) left=new tree;
      //cout<<" ";
      left->add(new_val, new_ind);
    } else {
      if (right==NULL) right=new tree;
      //cout <<" ";
      right->add(new_val, new_ind);
    }
  }
  return index;
}

//adds an element to the tree
//If the element is a duplicate, it is not added
//and the function returns the index value of the duplicate tree element
//Otherwise the method returns a copy of the index
template <class primitive>
long tree<primitive>::add(primitive new_val, long new_ind) {
  if (value==NULL) {
    value=new primitive;
    *value=new_val;
    index=new_ind;
    //cout<<*value<<"\n";
  } else {
    if (new_val < *value) {
      if (left==NULL) left=new tree;
      //cout<<" ";
      left->add(new_val, new_ind);
    } else {
      if (right==NULL) right=new tree;
      //cout <<" ";
      right->add(new_val, new_ind);
    }
  }
  return index;
}

template <class primitive>
//decomposes the tree into an array
//the size of the array must be known in advance
//so the number of elements in the tree must be kept track of
//externally
void tree<primitive>::decompose(primitive *data_arr, long *index_arr, long n, long &i) {

  if (left!=NULL) {
    left->decompose(data_arr, index_arr, n, i);
  }
  if (i >=n) return;
  data_arr[i]=*value;
  index_arr[i]=index;
  i++;
  if (right!=NULL) {
    right->decompose(data_arr, index_arr, n, i);
  }
}

template <class primitive>
tree<primitive>::~tree() {
  if (left != NULL) delete left;
  if (right != NULL) delete right;
  if (value != NULL) delete value;
}

//now to test the tree:
template <class primitive>
primitive *treesort(primitive *data, long n, long **index) {
  primitive *sorted_data;
  long *sorted_index;
  tree<primitive> sorter;
  long i;

  sorted_data=new primitive[n];
  sorted_index=new long[n];

  for (i=0; i<n; i++) {
    //cout<<data[i];
    sorter.add(data[i], i);
  }
  i=0;
  sorter.decompose(sorted_data, sorted_index, n, i);

  //for (i=0;i<n;i++) cout<<sorted_index[i]<<" ";
  //cout<<"\n";

  *index=sorted_index;

  return sorted_data;
}

template class tree<float>;
template class tree<int32_t>;
template class tree<int64_t>;
template class tree<double>;
template class tree<time_class>;
template class tree<string_petey>;

} //end namespace libpetey

#if TEST
#define N 5L

int main() {
  float data[N]={4.0, 2.0, 3.0, 1.0, 6.0};
  float *sorted;
  long *index;
  long i;

  for (i=0;i<N;i++) cout<<data[i]<<" ";
  cout<<"\n";

  sorted=sort(data, N, &index);

  for (i=0;i<N;i++) cout<<sorted[i]<<" ";
  cout<<"\n";
  for (i=0;i<N;i++) cout<<index[i]<<" ";
  cout<<"\n";
//  cout<<index;

  delete sorted;
  delete index;

  return 1;
}
#endif

