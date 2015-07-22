#include "string_petey.h"

namespace libpetey {

  template <class dt>
  long * heapsort_ptr(dt **data, long n) {

  //***********************************************************************
  //
  // purpose:	Performs a heap-sort on an array of variables.  Can be used
  //		to sort user defined data types (obviously) by defining '>'
  //		(greater than) and '<' (less than) operators for the desired
  //		class.
  //
  //		This version an array of longword integers giving the indices
  //		of the sorted array, while leaving the original array untouched.
  //
  // usage:	ind=heapsort(data, n)
  //
  // input/output:	data: an array of variables.  Can be used with any
  // 		user defined type, so long as the comparator operators are
  // 		defined for it.
  //
  //		n: the number of elements in the array.
  //
  // written by:	Peter Mills (peteymills@hotmail.com)
  //
  // history:	Created 2003-2-25 based on a heapsort function in IDL which
  //		in turn was loosely based on that in Numerical Recipes.
  //
  //************************************************************************

    long i, j, jold, k;
    long temp;
    long * ind;

    ind=new long[n];
    for (i=0;i<n;i++) ind[i]=i;

    k=n/2;

    //build the heap:
    for (i=k; i>=0; i--) {
      jold=i;
      j=i*2+1;

      while (j < n) {
        if (j < n-1) {
          if (*data[ind[j]] < *data[ind[j+1]]) j=j+1;
        }
        if (*data[ind[jold]] < *data[ind[j]]) {
          temp=ind[jold];
          ind[jold]=ind[j];
          ind[j]=temp;
        } else {
          break;
        }
        jold=j;
        j=j*2+1;
      }
    }

    //pull each element off the heap in turn:
    for (i=n-1; i>=1; i--) {
      temp=ind[i];
      ind[i]=ind[0];
      ind[0]=temp;

      jold=0;
      j=1;
      while (j < i) {
        if (j < i-1) {
          if (*data[ind[j]] < *data[ind[j+1]]) j=j+1;
        }
        if (*data[ind[jold]] < *data[ind[j]]) {
          temp=ind[jold];
          ind[jold]=ind[j];
          ind[j]=temp;
        } else {
          break;
        }
        jold=j;
        j=j*2+1;
      }
    }

    return ind;

  }

  template long *heapsort_ptr<string>(string_petey **data, long n);

} //end namespace libpetey

