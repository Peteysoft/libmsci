#include <math.h>
#include <iostream.h>
#include "dependent_dataset.h"
#include "simple_temp.h"

#pragma implementation

template class simple<float>;
template class simple<double>;
template class simple<long>;

template <class primitive>
simple<primitive>::simple() {
  data=NULL;
  n_data=0;
}

template <class primitive>
simple<primitive>::simple(const simple<primitive> &other) {
  long i;
  n_data=other.n_data;
  data=new primitive[n_data];
  for (i=0;i<n_data;i++) data[i]=other.data[i];
}

template <class primitive>
simple<primitive> simple<primitive>::operator =(const simple &other) {
  long i;
  if (data != NULL) delete data;

  n_data=other.n_data;
  data=new primitive[n_data];
  for (i=0;i<n_data;i++) data[i]=other.data[i];
}

template <class primitive>
simple<primitive>::~simple() {
  delete data;
}

template <class primitive>
long simple<primitive>::get (primitive &value, long ind) {

  if (ind < 0 || ind >= n_data) return -1;
  value=data[ind];

  return 0;
}

//adds a new value to the dataset
template <class primitive>
long simple<primitive>::add_el(primitive value) {
  long ind, i;
  primitive *new_data;
  dependent_dataset * this_dependent;

  if (n_data==0) {
    data=new primitive[0];
    data[0]=value;
    ind = -1;
  } else {
    ind = search(value);
    if (ind == -1) {
      //add the data value at the beginning of the array:
      new_data=new primitive[n_data+1];
      for (i=0;i<n_data;i++) new_data[i+1]=data[i];
      new_data[0]=value;
      delete data;
      data=new_data;
    } else if (data[ind]!=value) {
      //add the data value in the middle of the array:
      new_data=new primitive[n_data+1];
      for (i=0;i<=ind;i++) new_data[i]=data[i];
      new_data[ind+1]=value;
      for (i=ind+1;i<n_data;i++) new_data[i+1]=data[i];
      delete data;
      data=new_data;
    } else {
      return ind;
    }
  }

  n_data++;
  ind++;

  for (i=0;i<ndep;i++) {
    dependents[i]->insert(rank[i], ind);
  }

  return ind;
}

template <class primitive>
long simple<primitive>::del(primitive value) {
  long ind, i;
  primitive *new_data;

  ind=search(value);
  if (ind == -1) return -1;
  if (data[ind] != value) return -1;
  new_data=new primitive[n_data-1];
  for (i=0;i<ind;i++) new_data[i]=data[i];
  for (i=ind+1;i<n_data;i++) new_data[i-1]=data[i];
  delete data;
  data=new_data;
  n_data--;

  //delete a slab in each of the dependents:
  for (i=0;i<ndep;i++) {
    dependents[i]->del(rank[i], ind);
  }

  return ind;
}

//gets the zero-based index of the array element containing value
template <class primitive>
long simple<primitive>::search(primitive value) {
  //do a simple binary search:
  long first, last, mid;

  if (n_data == 0) return -1L;

  last=n_data-1;
  first=0;
  if (value>=data[last]) return last;
  if (value<data[first]) return -1L;
  if (value==data[first]) return first;
  do {
    mid=(last+first)/2;
    //cout << first << mid << last << "\n";
    if (value==data[mid]) return mid;
    if (value > data[mid]) {
      first=mid;
    } else {
      last=mid;
    }
  } while (last-first > 1L);
  //if (fabs(data[last]-value) > fabs(data[first]-value)) {
  //  last=first;
  //}
  return first;
}

simple<float>::simple() {
  n_data=0;
  data=NULL;
  type=SIMPLE_FLOAT;
  ndep=0;
  dependents=NULL;
  rank=NULL;
}

simple<double>::simple() {
  n_data=0;
  data=NULL;
  type=SIMPLE_DOUBLE;
  ndep=0;
  dependents=NULL;
  rank=NULL;
}

simple<long>::simple() {
  n_data=0;
  data=NULL;
  type=SIMPLE_LONG;
  ndep=0;
  dependents=NULL;
  rank=NULL;
}

