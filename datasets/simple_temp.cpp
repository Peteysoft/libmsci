#include <math.h>
#include <stdint.h>
//#include <stdio.h>
//#include <iostream.h>

#include <peteys_tmpl_lib.h>
#include <time_class.h>
#include <string_petey.h>

#include "dependent_dataset.h"
#include "simple_temp.h"

#pragma implementation

namespace libpetey {
namespace datasets {

template <class primitive>
ind_type simple<primitive>::search(primitive value) {
  last_search=(ind_type) bin_search(data, n_data, value, last_search);
  return last_search;
}

template <class primitive>
interpol_index simple<primitive>::interp(primitive value) {
  interpol_index ind;
//  printf("Last index searched: %d\n", last_search);
  ind=interpolate(data, n_data, value, last_search);
  last_search=(ind_type) ind;
  return ind;
}

template <class primitive>
ind_type simple<primitive>::get (primitive &value, ind_type ind) {

  if (ind < 0 || ind >= n_data) return SUBSCRIPT_OUT_OF_RANGE;
  value=data[ind];

  return 0;
}

template <class primitive>
simple<primitive>::simple() {
  set_type();
  data=NULL;
  n_data=0;
  this->ndep=0;
  dependents=NULL;
  this->rank=NULL;

  last_search=-1;
}

template<>
void simple<float>::set_type() {
  type=SIMPLE_FLOAT;
}

template<>
void simple<double>::set_type() {
  type=SIMPLE_DOUBLE;
}

/*
template<>
void simple<long>::set_type() {
  type=SIMPLE_LONG;
}
*/

template<>
void simple<int32_t>::set_type() {
  type=SIMPLE_INT32;
}

template<>
void simple<int64_t>::set_type() {
  type=SIMPLE_INT64;
}

template<>
void simple<time_class>::set_type() {
  type=SIMPLE_TIME;
}

template<>
void simple<string>::set_type() {
  type=SIMPLE_STRING;
}

template <class primitive>
simple<primitive>::simple(const simple<primitive> &other) {
  ind_type i;

  this->ndep=0;
  dependents=NULL;
  this->rank=NULL;

  n_data=other.n_data;
  data=new primitive[n_data];
  type=other.type;
  for (i=0;i<n_data;i++) data[i]=other.data[i];

  last_search=-1;
}

template <class primitive>
template <class p2>
simple<primitive> simple<primitive>::operator =(const simple<p2> &other) {
  ind_type i;

  if (data != NULL) delete [] data;

  if (this->ndep != 0) {
    delete [] dependents;
    delete [] this->rank;
  }
  this->ndep=0;
  dependents=NULL;
  this->rank=NULL;

  n_data=other.n_data;
  data=new primitive[n_data];
  //type=other.type;
  for (i=0;i<n_data;i++) data[i]=(primitive) other.data[i];
}

template <class primitive>
simple<primitive>::~simple() {
//  printf("Simple template destructor called\n");
/*  if (this->ndep != 0) {
    delete dependents;
    delete this->rank;
  }
*/
  if (data != NULL) delete [] data;
}

template <class primitive>
long simple<primitive>::read(FILE *fileptr) {
  long nread;

  nread=fread(&n_data, sizeof(n_data), 1, fileptr);
  if (data != NULL) delete [] data;
  data=new primitive[n_data];
  nread+=fread(data, sizeof(primitive)*n_data, 1, fileptr);

  return nread;
}

template <class primitive>
long simple<primitive>::write(FILE *fileptr) {
  long nwritten;

  nwritten=fwrite(&n_data, sizeof(n_data), 1, fileptr);
  nwritten+=fwrite(data, sizeof(primitive)*n_data, 1, fileptr);

  return nwritten;
}

template <>
void simple<int32_t>::print() {
  printf("Type: %d\n", (int32_t) type);
  printf("%ld dependents:\n", (int64_t) ndep);
  for (long i=0; i<ndep; i++) {
    printf("%d %d\n", (int32_t) rank[i], (int32_t) dependents[i]->type_of());
  }

  for (long i=0; i<n_data; i++) printf("%d ", data[i]);
  printf("\n");
}

template <>
void simple<int64_t>::print() {
  printf("Type: %d\n", (int32_t) type);
  printf("%ld dependents:\n", (int64_t) ndep);
  for (long i=0; i<ndep; i++) {
    printf("%d %d\n", (int32_t) rank[i], (int32_t) dependents[i]->type_of());
  }

  for (long i=0; i<n_data; i++) printf("%ld ", data[i]);
  printf("\n");
}

template <>
void simple<float>::print() {
  printf("Type: %d\n", (int32_t) type);
  printf("%ld dependents:\n", (int64_t) ndep);
  for (long i=0; i<ndep; i++) {
    printf("%d %d\n", (int32_t) rank[i], (int32_t) dependents[i]->type_of());
  }

  for (long i=0; i<n_data; i++) printf("%g ", data[i]);
  printf("\n");
}

template <>
void simple<double>::print() {
  printf("Type: %d\n", (int32_t) type);
  printf("%ld dependents:\n", (int64_t) ndep);
  for (long i=0; i<ndep; i++) {
    printf("%d %d\n", (int32_t) rank[i], (int32_t) dependents[i]->type_of());
  }

  for (long i=0; i<n_data; i++) printf("%lg ", data[i]);
  printf("\n");
}

template <>
void simple<time_class>::print() {
  char tstr[30];

  printf("Type: %d\n", (int32_t) type);
  printf("%ld dependents:\n", (int64_t) ndep);
  for (long i=0; i<ndep; i++) {
    printf("%d %d\n", (int32_t) rank[i], (int32_t) dependents[i]->type_of());
  }

  for (long i=0; i<n_data; i++) {
    data[i].write_string(tstr);
    printf("%s ", tstr);
  }
  printf("\n");
}

template <>
void simple<string_petey>::print() {
  char *str;

  printf("Type: %d\n", (int32_t) type);
  printf("%ld dependents:\n", (int64_t) ndep);
  for (long i=0; i<ndep; i++) {
    printf("%d %d\n", (int32_t) rank[i], (int32_t) dependents[i]->type_of());
  }

  for (long i=0; i<n_data; i++) {
    str=(char *) data[i];
    printf("%s ", str);
    delete [] str;
  }
  printf("\n");
}

//adds a new value to the dataset
template <class primitive>
ind_type simple<primitive>::add_el(primitive value) {
  ind_type ind, i;
  primitive *new_data;
  dependent_dataset * this_dependent;

  if (n_data==0) {
    data=new primitive[1];
    data[0]=value;
    ind = -1;
  } else {
    ind = search(value);
    if (ind == -1) {
      //add the data value at the beginning of the array:
      new_data=new primitive[n_data+1];
      for (i=0;i<n_data;i++) new_data[i+1]=data[i];
      new_data[0]=value;
      delete [] data;
      data=new_data;
    } else if (data[ind]!=value) {
      //add the data value in the middle of the array:
      new_data=new primitive[n_data+1];
      for (i=0;i<=ind;i++) new_data[i]=data[i];
      new_data[ind+1]=value;
      for (i=ind+1;i<n_data;i++) new_data[i+1]=data[i];
      delete [] data;
      data=new_data;
    } else {
      return ind;
    }
  }

  n_data++;
  ind++;
  last_search++;

  for (i=0;i<this->ndep;i++) {
    dependents[i]->insert(this->rank[i], ind, 1);
  }

  return ind;
}

template <class primitive>
ind_type simple<primitive>::del(primitive value) {
  ind_type ind, i;
  primitive *new_data;

  ind=search(value);
  if (ind == -1) return DATA_ELEMENT_NOT_FOUND;
  if (data[ind] != value) return DATA_ELEMENT_NOT_FOUND;
  new_data=new primitive[n_data-1];
  for (i=0;i<ind;i++) new_data[i]=data[i];
  for (i=ind+1;i<n_data;i++) new_data[i-1]=data[i];
  delete [] data;
  data=new_data;
  n_data--;

  //delete a slab in each of the dependents:
  for (i=0;i<this->ndep;i++) {
    dependents[i]->del(this->rank[i], ind, 1);
  }

  return ind;
}

template <class primitive>
ind_type simple<primitive>::get(primitive &val, interpol_index ind) {
  ind_type low, high;
  primitive val1;
  primitive val2;
  primitive val3;

  if (data == NULL) return NO_DATA;

  low=(ind_type) ind;
  if (low < 0) {
    low=0;
  } else if (low >= n_data-1) {
    low=n_data-2;
  }
  high=low+1;

  val1=data[low];
  val2=data[high];
  val3=val2-val1;
  val2=(primitive) (val3*(ind-(interpol_index) low));
  val=val1+val2;
  
  return low;
}

template<>
ind_type simple<string>::get(string &val, interpol_index ind) {
  ind_type median;

  if (data == NULL) return NO_DATA;

  median=(ind_type) (ind+0.5);
  if (median < 0) median=0; else if (median >= n_data) median=n_data-1;

  val=data[median];

  return median;
}

template <class primitive>
simple<primitive>::simple(primitive *new_data, ind_type n, char sorted) {
  ind_type i;
//  printf("Templated simple preload constructor called\n");

  set_type();

  n_data=n;
  data=new primitive[n];
  for (i=0; i<n; i++) data[i]=new_data[i];
  if (sorted != 1) heapsort_inplace(data, n_data);

  fprintf(stderr, "%d\n", type);

  last_search=-1;
}

template <class primitive>
simple<primitive>::simple(const primitive &x0, const primitive &xf, ind_type n) {
  ind_type i;
//  printf("Templated simple preload constructor called\n");

  set_type();

  n_data=n;
  data=new primitive[n];
  for (i=0; i<n; i++) data[i]=(primitive) (((double) xf- (double) x0)/(n-1)*i + (double) x0);

  fprintf(stderr, "%d\n", type);

  last_search=-1;
}

template <>
simple<time_class>::simple(const time_class &x0, const time_class &xf, ind_type n) {
  ind_type i;
  double diff;
//  printf("Templated simple preload constructor called\n");

  set_type();

  n_data=n;
  diff=xf.diff(x0);
  data=new time_class[n];
  for (i=0; i<n; i++) {
    data[i]=x0;
    data[i].add(i*diff/(n-1));
  }

  fprintf(stderr, "%d\n", type);

  last_search=-1;
}

template <>
simple<string_petey>::simple(const string_petey &x0, const string_petey &xf, ind_type n) {
}


template <class primitive>
template <class p2>
int simple<primitive>::operator == (const simple<p2> &other) {
  if (other.n_data != n_data) return 0;

  for (long i=0; i<n_data; i++) {
    if (data[i] != (primitive) other.data[i]) return 0;
  }

  return 1;
}

template class simple<float>;
template class simple<double>;
template class simple<int32_t>;
template class simple<int64_t>;
template class simple<time_class>;
template class simple<string_petey>;

template int simple<float>::operator == <float> (const simple<float> &other);
//template int simple<float>::operator == <double> (const simple<double> &other);
//template int simple<float>::operator == <long> (const simple<long> &other);

template int simple<double>::operator == <double> (const simple<double> &other);
//template int simple<double>::operator == <float> (const simple<float> &other);
//template int simple<double>::operator == <long> (const simple<long> &other);
//template int simple<double>::operator == <time> (const simple<time> &other);

//template int simple<long>::operator == <double> (const simple<double> &other);
//template int simple<long>::operator == <float> (const simple<float> &other);
template int simple<long>::operator == <long> (const simple<long> &other);

template int simple<time_class>::operator == <time_class> (const simple<time_class> &other);
//template int simple<time>::operator == <double> (const simple<double> &other);

template int simple<string>::operator == <string> (const simple<string> &other);

} //end namespace datasets
} //end namespace libpetey


