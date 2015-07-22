#ifndef PETEY_POINTER_H
#define PETEY_POINTER_H

#include "stdlib.h"

namespace libpetey {

  //class is not thread safe.
  //"smart" (ha, ha) pointer
  template <class type>
  class petey_pointer {

    protected:
      type *obj;
      int *nptr;

    public:
      petey_pointer() {
        obj=NULL;
        nptr=NULL;
      }

      //seems a bit stupid...  shouldn't we just type, *<this>=n instead?
      petey_pointer(const type &n) {
        nptr=new int;
        *nptr=1;
        obj=new type(n);
      }

      petey_pointer(type *n) {
        nptr=new int;
        *nptr=1;
        obj=n;
      }

      ~petey_pointer() {
        clean();
      }

      void clean() {
        if (obj!=NULL) {
          (*nptr)--;
          if ((*nptr) == 0) {
            delete obj;
	    delete nptr;
          }
        }
        obj=NULL;
        nptr=NULL;
      }

      type & operator * () {
        return *obj;
      }

      //if any changes are made to the original, make a copy:
      type * operator ->() {
        if ((*nptr)>1) {
          obj=new type(*obj);
          (*nptr)--;
          nptr=new int;
          *nptr=1;
        }
        return obj;
      }
    
      petey_pointer (const petey_pointer<type> &ot) {
        nptr=ot.nptr;
        obj=ot.obj;
        (*nptr)++;
      }

      petey_pointer<type> & operator = (const petey_pointer<type> &ot) {
        clean();
        nptr=ot.nptr;
        obj=ot.obj;
        (*nptr)++;
        return *this;
      }

      petey_pointer<type> & operator = (type *ot) {
        clean();
        obj=ot;
        nptr=new int;
        *nptr=1;
        return *this;
      }

      operator type * () {
        (*nptr)++;
        return obj;
      }

      petey_pointer<type> & operator = (const type &ot) {
        clean();
        obj=new type(ot);
        nptr=new int;
        *nptr=1;
        return *this;
      }

      void print() {
        printf("%d %p %p\n", *nptr, nptr, obj);
      }

  };

} //end namespace libstupid

#endif

