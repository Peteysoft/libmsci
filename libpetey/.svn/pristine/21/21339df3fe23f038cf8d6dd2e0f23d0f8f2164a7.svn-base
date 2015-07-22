#include <assert.h>
#include <stdlib.h>

#include "time_class.h"
#include "string_petey.h"
#include "linked.h"

#pragma implementation

namespace libpetey {

template <class ds>
linked_list<ds>::linked_list() {
  head=NULL;
  tail=head;
  current=head;
  n_elements=0;
}

template <class ds>
linked_list<ds>::~linked_list() {
  reset();
}

//reset the linked list:
template <class ds>
void linked_list<ds>::reset() {
  ll_element<ds> *prev;
  if (n_elements==0) return;

  //recursive algorithm eats up stack space
  //(can cause segmentation fault for large lists...):
  //head->del();
  //delete head;

  prev=head;
  current=head;
  for (long i=1; i<n_elements; i++) {
    current=current->next;
    delete prev;
    prev=current;
  }
  assert(current==tail);
  delete current;

  n_elements=0;
  head=NULL;
  tail=head;
  current=head;
}

template <class ds>
//returns the number of elements in the list if successful
//otherwise, returns -1
long linked_list<ds>::push(ds data) {
  ll_element<ds> *tween;
  if (head==NULL) {
    head=new ll_element<ds>(data);
    if (head==NULL) return -1;
    tail=head;
  } else {
    tween=head;
    head=new ll_element<ds>(data);
    if (head==NULL) {
      head=tween;
      return -1;
    }
    head->next=tween;
  }
  current=head;
  last=0;
  return ++n_elements;
}

template <class ds>
long linked_list<ds>::pop(ds &data) {
  ll_element<ds> *next;
  if (n_elements==0L) {
    return -1L;
  } else {
    data=head->data;
    next=head->next;
    delete head;
    head=next;
    current=head;
    last=0;
    return --n_elements;
  }
}

template <class ds>
long linked_list<ds>::add(ds data) {
  if (head==NULL) {
    head=new ll_element<ds>(data);
    if (head==NULL) return -1;
    tail=head;
  } else {
    tail->next=new ll_element<ds>(data);
    if (tail->next==NULL) return -1;
    tail=tail->next;
  }
  current=tail;
  last=++n_elements;
  return n_elements;
}

template <class ds>
ds *linked_list<ds>::make_array(long &n) {
  ds *array;
  ll_element<ds> *current;
  long i;
  if (n_elements==0) {
    n=0;
    return NULL;
  }
  array=new ds[n_elements];
  current=head;
  for (i=0;i<n_elements;i++) {
    array[i]=current->data;
    current=current->next;
  }
  n=n_elements;
  return array;
}

template <class ds>
ds linked_list<ds>::del (long n) {
  ds data;
  long i, n1;  
  ll_element<ds> * next;

  if (head==NULL) return data;

  if (n==0) {
    data=head->data;
    next=head->next;
    delete head;
    head=next;
    n_elements--;
    return data;
  }
 
  if (n>=last+1) {
    n1=n-last-1;
  } else {
    n1=n-1;
    current=head;
    last=0;
  }
  i=0;
  while (i < n1 && current->next != NULL) {
    current=current->next;
    i++;
  }
  last+=i;

  if (current->next != NULL) {
    if (current->next == tail) tail=current;
    next=current->next->next;
    data=current->next->data;
    delete current->next;
    current->next=next;
    n_elements--;
  } else {
    printf("Warning: attempting to delete out-of-bounds element\n");
  }

  return data;

}

template <class ds>
ds linked_list<ds>::operator [] (long n) {
  long i, n1;  
 
  if (n>=last) {
    n1=n-last;
  } else {
    n1=n;
    current=head;
    last=0;
  }
  i=0;
  while (i < n1 && current->next != NULL) {
    current=current->next;
    i++;
  }

  last+=i;

  return current->data;
}

template class linked_list<float>;
template class linked_list<long>;
template class linked_list<int>;
template class linked_list<double>;
template class linked_list<time_class>;
template class linked_list<string_petey>;
template class linked_list<char *>;
template class linked_list<char>;
/*
template class ll_element<float>;
template class ll_element<long>;
template class ll_element<double>;
template class ll_element<time_class>;
template class ll_element<string_petey>;
*/

} //end namespace libpetey
