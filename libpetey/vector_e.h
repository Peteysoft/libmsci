#ifndef VECTOR_S_H
#define VECTOR_S_H 1

#include <assert.h>

#include "petey_pointer.h"

namespace libpetey {

  //an extensible vector class:
  //(why isn't this a parent to vector_s?  this-> would have to be added to all
  //the fields plus all the constructors will need to be duplicated anyway)
  template <class type>
  class vector_e {
    protected:
      int nel;
      int array_size;
      type * data;

    public:
      type missing;

      vector_e();
      vector_e(int n, type m);
      vector_e(type *dt, int n, int ncp_flag=0);
      vector_e(type **dt, int n);
      ~vector_e();

      void resize (int n);
      void resize (int n, type m);

      vector_e<type> & operator = (vector_e<type> &other);
      vector_e(vector_e<type> &other);

      type & operator [] (int ind);

      int size();

      operator type * ();

  };

  template <class type>
  vector_e<type>::vector_e() {
    nel=0;
    array_size=1;
    data=new type[array_size];
    //missing=(type) 0.;
  }

  template <class type>
  vector_e<type>::vector_e(int n, type m) {
    nel=n;
    array_size=n;
    data=new type[nel];
    missing=m;
  }

  template <class type>
  vector_e<type>::vector_e(type *dt, int n, int ncp_flag) {
    nel=n;
    array_size=n;
    if (ncp_flag) {
      data=dt;
    } else {
      data=new type[nel];
      for (int i=0; i<n; i++) {
        data[i]=dt[i];
      }
    }
  }

  template <class type>
  vector_e<type>::~vector_e() {
    delete [] data;
  }

  template <class type>
  vector_e<type>::vector_e(vector_e<type> &other) {
    data=new type[other.nel];
    nel=other.nel;
    array_size=nel;
    missing=other.missing;
    for (int i=0; i<nel; i++) data[i]=other.data[i];
  }
    
  template <class type>
  vector_e<type> & vector_e<type>::operator = (vector_e<type> &other) {
    delete [] data;
    data=new type[other.nel];
    nel=other.nel;
    array_size=nel;
    missing=other.missing;
    for (int i=0; i<nel; i++) data[i]=other.data[i];
    return *this;		//can this work??
  }

  template <class type>
  void vector_e<type>::resize(int n) {
    if (n>array_size) {
      type *newdata=new type[n];
      array_size=n;
      for (int i=0; i<nel; i++) newdata[i]=data[i];
      for (int i=nel; i<array_size; i++) newdata[i]=missing;
      delete [] data;
      data=newdata;
    }
    nel=n;
  }
    
  template <class type>
  void vector_e<type>::resize(int n, type m) {
    missing=m;
    resize(n);
  }

  template <class type>
  type & vector_e<type>::operator [] (int ind) {
    assert(ind>-1);
    if (ind>=nel) {
      if (ind>=array_size) {
        type *newdata;
        array_size=(ind/array_size+1)*array_size;
        newdata=new type[array_size];
        for (int i=0; i<nel; i++) newdata[i]=data[i];
        for (int i=nel; i<=ind; i++) newdata[i]=missing;
        delete [] data;
        data=newdata;
      }
      nel=ind+1;
    }
    return data[ind];
  }

  template <class type>
  int vector_e<type>::size () {
    return nel;
  }

  template <class type>
  vector_e<type>::operator type * () {
    return data;
  }
 
} //end fucking namespace libwhathisface...

#endif

