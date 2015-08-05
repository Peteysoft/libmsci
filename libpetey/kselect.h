#ifndef KSELECT_H__DEFINED

#include <vector>
#include "vector_s.h"
#include "tree_lg.h"
#include "tree_lgi.h"

using namespace std;

namespace libpetey {

  template <class type>
  class kselect_base {
    protected:
      long k;
    public:
      kselect_base();
      kselect_base(long k1);
      virtual ~kselect_base();
      //add an element, returns current size of data structure:
      virtual long add(type val)=0;
      //marshal out the k-least values:
      virtual void get(type * kleast)=0;
  };
    
  template <class type>
  class kselect_naive:public kselect_base<type> {
    protected:
      type *select;
      long ncur;
    public:
      kselect_naive(long k1);
      virtual ~kselect_naive();
      //add an element, returns current size of data structure:
      virtual long add(type val);
      //marshal out the k-least values:
      virtual void get(type * kleast);
  };

  template <class type>
  class kselect_tree:public kselect_base<type> {
    protected:
      tree_lg<type> data;
    public:
      kselect_tree(long k1);
      virtual ~kselect_tree();
      virtual long add(type val);
      virtual void get(type * kleast);
  };

  template <class type>
  class kselect_heap:public kselect_base<type> {
    protected:
      vector<type> *data;
    public:
      kselect_heap(long k1);
      virtual ~kselect_heap();
      virtual long add(type val);
      virtual void get(type * kleast);
  };

  template <class type>
  class kselect_quick:public kselect_base<type> {
    protected:
      vector_s<type> *data;
      long ncur;
    public:
      kselect_quick(long k1);
      //since quick-sort version must store all the data, pre-sets data structure to size n:
      kselect_quick(long k1, long n);
      virtual ~kselect_quick();
      virtual long add(type val);
      virtual void get(type * kleast);
  };

  template <class type>
  class kiselect_base {
    protected:
      long ncur;
      long k;
    public:
      kiselect_base();
      kiselect_base(long k1);
      virtual ~kiselect_base();
      //add an element, set associated index to element count; 
      //return current size of data structure:
      virtual long add(type val)=0;
      //add an element and set the index, returns current size of data structure:
      virtual long add(type val, long ind)=0;
      //marshall out the k-least and associated indices:
      virtual void get(type * kleast, long *ind)=0;

      //test the selection algo:
      int test(long n);		//number of elements

  };
    
  template <class type>
  class kiselect_naive:public kiselect_base<type> {
    protected:
      type *select;
      long *ind;
    public:
      kiselect_naive(long k1);
      virtual ~kiselect_naive();
      //add an element, returns current size of data structure:
      virtual long add(type val);
      virtual long add(type val, long i);
      //marshal out the k-least values:
      virtual void get(type * kleast, long *j);
  };

  template <class type>
  class kiselect_tree:public kiselect_base<type> {
    protected:
      tree_lgi<type> data;
    public:
      kiselect_tree(long k1);
      virtual ~kiselect_tree();
      virtual long add(type val);
      virtual long add(type val, long ind);
      virtual void get(type * kleast, long *ind);
  };

  template <class type>
  struct heap_ind_el {
    type data;
    long ind;
  };

  template <class type>
  int heap_ind_el_comp(const heap_ind_el<type> &v1, heap_ind_el<type> &v2);

  template <class type>
  class kiselect_heap:public kiselect_base<type> {
    protected:
      vector<heap_ind_el<type> > *data;
    public:
      kiselect_heap(long k1);
      virtual ~kiselect_heap();
      virtual long add(type val);
      virtual long add(type val, long ind);
      virtual void get(type * kleast, long *ind);
  };

  template <class type>
  class kiselect_quick:public kiselect_base<type> {
    protected:
      vector_s<type> *data;
      vector_s<long> *ind;
    public:
      kiselect_quick(long k1);
      //since quick-sort version must store all the data, pre-sets data structure to size n:
      kiselect_quick(long k1, long n);
      virtual ~kiselect_quick();
      virtual long add(type val);
      virtual long add(type val, long ind2);
      virtual void get(type * kleast, long *ind2);
  };

} //end namespace libpetey

#endif

