#ifndef KEXTREME_INCLUDED
#define KEXTREME_INCLUDED

namespace libpetey {

  template<class type>
  void kleast(type *data, long n, long k, type *result);

  template<class type>
  void kgreatest(type *data, long n, long k, type *result);

  template<class type>
  void kleast(type *data, long n, long k, type *result, long *ind);

  template<class type>
  void kgreatest(type *data, long n, long k, type *result, long *ind);

}

#endif

