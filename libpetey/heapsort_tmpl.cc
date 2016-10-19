
namespace libpetey {

template <class dt>
void heapsort_inplace(dt *data, long n) {

//***********************************************************************
//
// purpose:	Performs a heap-sort on an array of variables.  Can be used
//		to sort user defined data types (obviously) by defining '>'
//		(greater than) and '<' (less than) operators for the desired
//		class.
//
//		This version performs an in-place sort, requiring no extra
//		storage.
//
// usage:	heapsort_inplace(data, n)
//
// input/output:	data: an array of variables.  Can be used with any
// 		user defined type, so long as the comparator operators are
// 		defined for it.
//
//		n: the number of elements in the array.
//
// written by:	Peter Mills (peteymills@hotmail.com)
//
// history:	Created 2003-2-23 based on a heapsort function in IDL which
//		in turn was loosely based on that in Numerical Recipes.
//
//************************************************************************

  long i, j, jold, k;
  dt temp;

  k=n/2;
  
  //build the heap:
  for (i=k; i>=0; i--) {
    jold=i;
    j=i*2+1;
    
    while (j < n) {
      if (j < n-1) {
        if (data[j] < data[j+1]) j=j+1;
      }
      if (data[jold] < data[j]) {
        temp=data[jold];
        data[jold]=data[j];
        data[j]=temp;
      } else {
        break;
      }
      jold=j;
      j=j*2+1;
    }
  }
  
  //pull each element off the heap in turn:
  for (i=n-1; i>=1; i--) {
    temp=data[i];
    data[i]=data[0];
    data[0]=temp;
    
    jold=0;
    j=1;
    while (j < i) {
      if (j < i-1) {
        if (data[j] < data[j+1]) j=j+1;
      }
      if (data[jold] < data[j]) {
        temp=data[jold];
        data[jold]=data[j];
        data[j]=temp;
      } else {
        break;
      }
      jold=j;
      j=j*2+1;
    }
  }
  
}
  
template <class dt>
void heapsort(dt *data, long *ind, long n) {

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
// usage:	heapsort(data, ind, n)
//
// input/output:	data: an array of variables. 
//
//		ind: set of indices that sort the array
//
//		n: the number of elements in the array.
//
// written by:	Peter Mills (peteymills@hotmail.com)
//
// history:	Created 2003-2-25 based on a heapsort function in IDL which
//		in turn was loosely based on that in Numerical Recipes.
//		2004-1-21 PM: created new version that takes indices as argument
//				instead of returning them.
//
//************************************************************************

  long i, j, jold, k;
  long temp;

  for (i=0;i<n;i++) ind[i]=i;

  k=n/2;

  //build the heap:
  for (i=k; i>=0; i--) {
    jold=i;
    j=i*2+1;

    while (j < n) {
      if (j < n-1) {
        if (data[ind[j]] < data[ind[j+1]]) j=j+1;
      }
      if (data[ind[jold]] < data[ind[j]]) {
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
        if (data[ind[j]] < data[ind[j+1]]) j=j+1;
      }
      if (data[ind[jold]] < data[ind[j]]) {
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

}


template <class dt>
long * heapsort(dt *data, long n) {
  long * ind;
  ind=new long[n];
  heapsort(data, ind, n);
  return ind;
}


//because we're too lazy to write a heapsort routine that sorts in descending order:
template <class dt>
void reverse(dt *data, long n) {
  long i, j, k;
  dt temp;

  k=n/2;
  for (i=0;i<k;i++) {
    j=n-i-1;
    temp=data[i];
    data[i]=data[j];
    data[j]=temp;
  }
}

//maps a vector to a set of new positions given by a vector of indices:
template <class dt>
dt * map_vector(dt * vector, long * indices, long n) {
  long i;
  dt *result;

  result=new dt[n];
  for (i=0;i<n;i++) result[i]=vector[indices[i]];

  return result;
}

//maps a vector to a set of new positions given by a vector of indices:
template <class dt>
void map_vector_inplace(dt * vector, long * indices, long n) {
  long i;
  dt *result;

  result=new dt[n];
  for (i=0;i<n;i++) result[i]=vector[indices[i]];
  for (i=0;i<n;i++) vector[i]=result[i];
  delete [] result;

}

} //end namespace libpetey

