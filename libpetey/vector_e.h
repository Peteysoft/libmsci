#ifndef VECTOR_E_H
#define VECTOR_E_H 1

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

      void set_default_missing();
    public:
      type missing;

      vector_e();
      vector_e(int n);
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

      void print(FILE *fs);

  };

  template <class integer>
  inline integer missing_val_signed_int() {
    integer missing=1;
    missing=1;
    for (int i=1; i<sizeof(long); i++) missing=missing*256;
    missing*=128;
  }
 
} //end namespace libpetey

#endif

