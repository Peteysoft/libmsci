#ifndef VECTOR_S_H
#define VECTOR_S_H 1

#include <stdio.h>

#include "petey_pointer.h"

namespace libpetey {

  //sortable vector class

  //a vector class, designed for (at the moment) removing duplicates
  //within a list through sorting...
  template <class type>
  class vector_s {
    protected:
      int nel;
      int array_size;
      type * data;

    public:
      type missing;

      vector_s();
      vector_s(int n);
      vector_s(int n, type m);
      vector_s(type *dt, int n, int ncp_flag=0);
      vector_s(type **dt, int n);
      ~vector_s();

      void set_default_missing();

      void resize (int n);
      void resize (int n, type m);

      vector_s<type> & operator = (vector_s<type> &other);
      vector_s(vector_s<type> &other);

      type & operator [] (int ind);

      int operator == (vector_s<type> &other);
      int operator != (vector_s<type> &other);
      int operator > (vector_s<type> &other);
      int operator < (vector_s<type> &other);
      int operator >= (vector_s<type> &other);
      int operator <= (vector_s<type> &other);

      vector_s<type> & operator + (vector_s<type> &other);
      vector_s<type> & operator - (vector_s<type> &other);
      //vector_s<type> & operator * (vector_s<type> &other);
      //vector_s<type> & operator / (vector_s<type> &other);

      vector_s<type> & operator + (type &b);
      vector_s<type> & operator - (const type &b);
      vector_s<type> & operator * (const double &k);
      vector_s<type> & operator / (const double &k);

      int size();

      operator type * ();

      void print(FILE *fs);
  };

}

#endif

