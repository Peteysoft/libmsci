#include <stdlib.h>
#include "tree_lg.h"
#include "tree_lgi.h"

#include "peteys_tmpl_lib.h"

namespace libpetey {

  template<class type>
  void kleast(type *data, long n, long k, type *result) {
    tree_lg<type> sorter;

    for (long i=0L; i<k; i++) {
      sorter.add(data[i]);
    }

    for (long i=k; i<n; i++) {
      sorter.add(data[i]);
      sorter.delete_greatest();
    }

    sorter.decompose(result, k);
  }

  template<class type>
  void kgreatest(type *data, long n, long k, type *result) {
    tree_lg<type> sorter;

    for (long i=0L; i<k; i++) {
      sorter.add(data[i]);
    }

    for (long i=k; i<n; i++) {
      sorter.add(data[i]);
      sorter.delete_least();
    }

    sorter.decompose(result, k);
  }

  template void kleast<float>(float *data, long n, long k, float *result);
  template void kgreatest<float>(float *data, long n, long k, float *result);

  template void kleast<double>(double *data, long n, long k, double *result);
  template void kgreatest<double>(double *data, long n, long k, double *result);

  template<class type>
  void kleast(type *data, long n, long k, type *result, long *ind) {
    tree_lgi<type> sorter;

    for (long i=0L; i<k; i++) {
      sorter.add(data[i], i);
    }

    for (long i=k; i<n; i++) {
      sorter.add(data[i], i);
      sorter.delete_greatest();
    }

    sorter.decompose(result, ind, k);

  }

  template<class type>
  void kgreatest(type *data, long n, long k, type *result, long *ind) {
    tree_lgi<type> sorter;

    for (long i=0L; i<k; i++) {
      sorter.add(data[i], i);
    }

    for (long i=k; i<n; i++) {
      sorter.add(data[i], i);
      sorter.delete_least();
    }

    sorter.decompose(result, ind, k);

  }

  template void kleast<float>(float *data, long n, long k, float *result, long *ind);
  template void kgreatest<float>(float *data, long n, long k, float *result, long *ind);

  template void kleast<double>(double *data, long n, long k, double *result, long *ind);
  template void kgreatest<double>(double *data, long n, long k, double *result, long *ind);
  
  template <typename type>
  int verify_kleast(type *list, long n, type *kl, long k) {
    type gl;		//greatest of the least
    type count=0;		//number in list greater than gl
    int err=0;

    heapsort_inplace(list, n);
    heapsort_inplace(kl, k);

    for (long i=0; i<k; i++) {
      if (list[i]!=kl[i]) {
        err=1;
	break;
      }
    }

    return err;
  }

  template int verify_kleast<float>(float *, long, float *, long);
  template int verify_kleast<double>(double *, long, double *, long);

} //end namespace libpetey

