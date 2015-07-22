#include <assert.h>
#include <stdio.h>
#include <math.h>

namespace libpetey {

  template <class type>
  vector_e<type>::vector_e() {
    nel=0;
    array_size=1;
    data=new type[array_size];
    set_default_missing();
    //missing=(type) 0.;
  }

  template <class type>
  vector_e<type>::vector_e(int n) {
    nel=n;
    array_size=n;
    data=new type[nel];
    set_default_missing();
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
  vector_e<type>::vector_e(vector_s<type> &other) {
    data=new type[other.nel];
    nel=other.nel;
    array_size=nel;
    missing=other.missing;
    for (int i=0; i<nel; i++) data[i]=other.data[i];
  }
    
  template <class type>
  vector_e<type> & vector_e<type>::operator = (vector_s<type> &other) {
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

