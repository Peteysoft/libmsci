#ifndef LINKED_H_INCLUDED
#define LINKED_H_INCLUDED 1

#pragma interface

namespace libpetey {

  template <class ds>
  class ll_element {
    public:
      ds data;
      ll_element<ds> *next;
      ll_element(ds new_data);	//constructor--adds data
      ~ll_element();
      void del();			//deletes all elements ahead of it
  };

  template <class ds>
  ll_element<ds>::ll_element(ds new_data) {
    data=new_data;
    next=NULL;
  }

  template <class ds>
  ll_element<ds>::~ll_element() {
  }

  template <class ds>
  void ll_element<ds>::del() {
    if (next != NULL) next->del();
    delete next;
  }

  template <class ds>
  class linked_list {
    protected:
      ll_element<ds> *head;
      ll_element<ds> *tail;
      ll_element<ds> *current;	//last element searched
      long n_elements;
      long last;			//last index searched
    public:
      linked_list();
      ~linked_list();
      long size_of();
      long push(ds data);		//pushes the first element
      long pop(ds &data);		//pops the first element
      long add(ds data);		//adds an element at the end
      ds *make_array(long &n);	//converts to an array

      ds del(long n);	//deletes specified element

      ds operator [] (long n);	//gets the specified element
      void reset();
  };

  template <class ds>
  inline long linked_list<ds>::size_of() {
    return n_elements;
  }

} //end namespace libpetey

#endif

