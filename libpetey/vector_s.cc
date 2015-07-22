#include <assert.h>
#include <stdio.h>
#include <math.h>

#include "time_class.h"
#include "error_codes.h"

#include "vector_s.h"

//#define TEST_VECTOR_S 

//great, a "sortable" vector

namespace libpetey {

  template <class type>
  vector_s<type>::vector_s() {
    nel=0;
    array_size=1;
    data=new type[array_size];
    set_default_missing();
    //missing=(type) 0.;
  }

  template <class type>
  vector_s<type>::vector_s(int n) {
    nel=n;
    array_size=n;
    data=new type[nel];
    set_default_missing();
  }

  template <class type>
  vector_s<type>::vector_s(int n, type m) {
    nel=n;
    array_size=n;
    data=new type[nel];
    missing=m;
  }

  template <class type>
  vector_s<type>::vector_s(type *dt, int n, int ncp_flag) {
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
  vector_s<type>::~vector_s() {
    delete [] data;
  }

  template <class type>
  vector_s<type>::vector_s(vector_s<type> &other) {
    data=new type[other.nel];
    nel=other.nel;
    array_size=nel;
    missing=other.missing;
    for (int i=0; i<nel; i++) data[i]=other.data[i];
  }
    
  template <class type>
  vector_s<type> & vector_s<type>::operator = (vector_s<type> &other) {
    delete [] data;
    data=new type[other.nel];
    nel=other.nel;
    array_size=nel;
    missing=other.missing;
    for (int i=0; i<nel; i++) data[i]=other.data[i];
    return *this;		//can this work??
  }

  template <class type>
  void vector_s<type>::resize(int n) {
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
  void vector_s<type>::resize(int n, type m) {
    missing=m;
    resize(n);
  }

  template <class type>
  type & vector_s<type>::operator [] (int ind) {
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
  vector_s<type> & vector_s<type>::operator + (vector_s<type> &other) {
    vector_s<type> *result;
    if (nel != other.nel) exit(SAMPLE_COUNT_MISMATCH);
    result=new vector_s<type>(nel, missing);
    for (int i=0; i<nel; i++) {
      if (data[i]!=missing && other.data[i]!=missing) {
        result->data[i]=data[i]+other.data[i];
      } else {
        result->data[i]=missing;
      }
    }
    return *result;		//can this work??
  }
    
  template <class type>
  vector_s<type> & vector_s<type>::operator - (vector_s<type> &other) {
    vector_s<type> *result;
    if (nel != other.nel) exit(SAMPLE_COUNT_MISMATCH);
    result=new vector_s<type>(nel, missing);
    for (int i=0; i<nel; i++) {
      if (data[i]!=missing && other.data[i]!=missing) {
        result->data[i]=data[i]-other.data[i];
      } else {
        result->data[i]=missing;
      }
    }
    return *result;		//can this work??
  }

  template <class type>
  vector_s<type> & vector_s<type>::operator + (type &c) {
    vector_s<type> *result;
    result=new vector_s<type>(nel, missing);
    for (int i=0; i<nel; i++) {
      if (data[i]!=missing) {
        result->data[i]=data[i]+c;
      } else {
        result->data[i]=missing;
      }
    }
    return *result;		//can this work??
  }
    
  template <class type>
  vector_s<type> & vector_s<type>::operator - (const type &c) {
    vector_s<type> *result;
    result=new vector_s<type>(nel, missing);
    for (int i=0; i<nel; i++) {
      if (data[i]!=missing) {
        result->data[i]=data[i]-c;
      } else {
        result->data[i]=missing;
      }
    }
    return *result;		//can this work??
  }

  template <class type>
  vector_s<type> & vector_s<type>::operator * (const double &c) {
    vector_s<type> *result;
    result=new vector_s<type>(nel, missing);
    for (int i=0; i<nel; i++) {
      if (data[i]!=missing) {
        result->data[i]=data[i]*c;
      } else {
        result->data[i]=missing;
      }
    }
    return *result;		//can this work??
  }
    
  template <class type>
  vector_s<type> & vector_s<type>::operator / (const double &c) {
    vector_s<type> *result;
    result=new vector_s<type>(nel, missing);
    for (int i=0; i<nel; i++) {
      if (data[i]!=missing) {
        result->data[i]=data[i]/c;
      } else {
        result->data[i]=missing;
      }
    }
    return *result;		//can this work??
  }

  template <class type>
  int vector_s<type>::operator == (vector_s<type> &other) {
    assert(nel == other.nel);
    for (int i=0; i<nel; i++) {
      if (other.data[i]==other.missing || data[i]==missing) break;
      if (other.data[i] != data[i]) return 0;
    }
    return 1;
  }

  template <class type>
  int vector_s<type>::operator > (vector_s<type> &other) {
    assert(nel == other.nel);
    for (int i=0; i<nel; i++) {
      if (other.data[i]==other.missing || data[i]==missing) break;
      if (data[i] > other.data[i]) return 1;
		else if (data[i] < other.data[i]) return 0;
    }
    return 0;
  }

  template <class type>
  int vector_s<type>::operator != (vector_s<type> &other) {
    assert(nel == other.nel);
    for (int i=0; i<nel; i++) {
      if (other.data[i]==other.missing || data[i]==missing) break;
      if (other.data[i] != data[i]) return 1;
    }
    return 0;
  }

  template <class type>
  int vector_s<type>::operator < (vector_s<type> &other) {
    assert(nel == other.nel);
    for (int i=0; i<nel; i++) {
      if (other.data[i]==other.missing || data[i]==missing) break;
      if (data[i] < other.data[i]) return 1;
		else if (data[i] > other.data[i]) return 0;
    }
    return 0;
  }

  template <class type>
  int vector_s<type>::operator >= (vector_s<type> &other) {
    assert(nel == other.nel);
    for (int i=0; i<nel; i++) {
      if (other.data[i]==other.missing || data[i]==missing) break;
      if (data[i] > other.data[i]) return 1;
		else if (data[i] < other.data[i]) return 0;
    }
    return 1;
  }

  template <class type>
  int vector_s<type>::operator <= (vector_s<type> &other) {
    assert(nel == other.nel);
    for (int i=0; i<nel; i++) {
      if (other.data[i]==other.missing || data[i]==missing) break;
      if (data[i] < other.data[i]) return 1;
		else if (data[i] > other.data[i]) return 0;
    }
    return 1;
  }

  template <class type>
  int vector_s<type>::size () {
    return nel;
  }

  template <class type>
  vector_s<type>::operator type * () {
    return data;
  }

  template <>
  void vector_s<float>::print (FILE *fs) {
    for (int i=0; i<nel; i++) {
      printf("%g ", data[i]);
    }
    printf("\n");
  }

  template <>
  void vector_s<double>::print (FILE *fs) {
    for (int i=0; i<nel; i++) {
      printf("%lg ", data[i]);
    }
    printf("\n");
  }

  template <>
  void vector_s<float>::set_default_missing () {
    missing=NAN;
  }

  template <>
  void vector_s<double>::set_default_missing () {
    missing=NAN;
  }

  template <>
  void vector_s<long>::set_default_missing () {
    missing=1;
    for (int i=1; i<sizeof(long); i++) missing=missing*256;
    missing*=128;
  }

  template <>
  void vector_s<int64_t>::set_default_missing () {
    missing=1;
    for (int i=1; i<sizeof(int64_t); i++) missing=missing*256;
    missing*=128;
  }

  template <>
  void vector_s<int32_t>::set_default_missing () {
    missing=2147483648;
    missing=1;
    for (int i=1; i<sizeof(int32_t); i++) missing=missing*256;
    missing*=128;
  }

  template <>
  void vector_s<time_class>::set_default_missing () {
    missing.init(-1000, 1, 1, 0, 0, 0);
  }


  template class vector_s<int32_t>;
  template class vector_s<long>;
  template class vector_s<int64_t>;
  template class vector_s<float>;
  template class vector_s<double>;
  template class vector_s<time_class>;

#ifdef TEST_VECTOR_S

  int main () {
    vector_s<float> a;
    vector_s<float> b;
    vector_s<float> c;
    vector_s<float> d;

    a[0]=1;
    a[1]=2;
    a[2]=3;
    b[0]=2;
    b[1]=1;
    b[2]=3;
    c[0]=1;
    c[1]=1;
    c[2]=1;

    d=a*b+c;

    printf("%f %f %f\n", d[0], d[1], d[2]);

  }
#endif
 
} //end fucking namespace libwhathisface...
