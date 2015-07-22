#include <stdio.h>

template<class type>
void kleast(type *data, long n, long k, type *result) {
  long ind;
  type swap;

  result[0]=data[0];
  for (long i=1; i<k; i++) {
    if (data[i] > result[0]) {
      result[i]=result[0];
      result[0]=data[i];
    } else {
      result[i]=data[i];
    }
    //printf("%f\n", result[i]);
  }

  for (long i=k; i<n; i++) {
    if (data[i]<result[0]) {
      result[0]=data[i];
      ind=0;
      for (long j=1; j<k; j++) {
        if (result[j]>result[ind]) {
	  ind=j;
	}
        //printf("%f\n", result[j]);
      }
      swap=result[0];
      result[0]=result[ind];
      result[ind]=swap;
    }
  }
}

template<class type>
void kgreatest(type *data, long n, long k, type *result) {

  result[0]=data[0];
  for (long i=1; i<k; i++) {
    if (data[i] < result[0]) {
      result[i]=result[0];
      result[0]=data[i];
    } else {
      result[i]=data[i];
    }
  }

  for (long i=k; i<n; i++) {
    if (data[i]>result[0]) {
      result[0]=data[i];
    }
  }
}

template<class type>
void kleast(type *data, long n, long k, type *result, long *ind);

template<class type>
void kgreatest(type *data, long n, long k, type *result, long *ind);

template void kleast<float>(float *data, long n, long k, float *result);
template void kgreatest<float>(float *data, long n, long k, float *result);

