
namespace libpetey {

void heapsort_inplace(void **data, long n, void *func) {

//***********************************************************************
//
// purpose:	Performs a heap-sort on an array of variables.  Can be used
//		to sort user defined data types by passing the comparison
//              function in the third argument.
//
//		This version performs an in-place sort, requiring no extra
//		storage.
//
// usage:	heapsort_inplace(data, n, )
//
// input/output:	data: an array of variables. 
//
//		n: the number of elements in the array.
//
//		func: comparison function; should be of the form, 
//			int comp(type *el1, type *el2);
//			- returns -1 for lt, 0 for eq and 1 for gt
//
// written by:	Peter Mills (peteymills@hotmail.com)
//
// history:	Created 2003-2-23 based on a heapsort function in IDL which
//		in turn was loosely based on that in Numerical Recipes.
//
//************************************************************************

  long i, j, jold, k;
  void *temp;

  int (*comp) (void *, void *)=(int (*) (void *, void *)) func;

  k=n/2;
  
  //build the heap:
  for (i=k; i>=0; i--) {
    jold=i;
    j=i*2+1;
    
    while (j < n) {
      if (j < n-1) {
        if ((*comp)(data[j], data[j+1])<0) j=j+1;
      }
      if ((*comp)(data[jold], data[j])<1) {;
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
        if ((*comp)(data[j], data[j+1])<0) j=j+1;
      }
      if ((*comp)(data[jold], data[j]) < 1) {
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
  
void heapsort(void **data, long *ind, long n, void *func) {

//***********************************************************************
//
// purpose:	Performs a heap-sort on an array of variables.  Can be used
//		to sort user defined data types by passing the comparison
//              function in the third argument.
//
//		This version an array of longword integers giving the indices
//		of the sorted array, while leaving the original array untouched.
//
// usage:	ind=heapsort(data, n)
//
// input/output:	data: an array of variables.
//
//		n: the number of elements in the array.
//
//		func: comparison function; should be of the form, 
//			int comp(type *el1, type *el2);
//			- returns -1 for lt, 0 for eq and 1 for gt
//
// written by:	Peter Mills (peteymills@hotmail.com)
//
// history:	Created 2003-2-25 based on a heapsort function in IDL which
//		in turn was loosely based on that in Numerical Recipes.
//
//************************************************************************

  long i, j, jold, k;
  long temp;

  int (*comp) (void *, void *)=(int (*) (void *, void *)) func;

  for (i=0;i<n;i++) ind[i]=i;

  k=n/2;

  //build the heap:
  for (i=k; i>=0; i--) {
    jold=i;
    j=i*2+1;

    while (j < n) {
      if (j < n-1) {
        if ((*comp)(data[ind[j]], data[ind[j+1]]) < 0) j=j+1;
      }
      if ((*comp) (data[ind[jold]], data[ind[j]]) < 0) {
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
        if ((*comp) (data[ind[j]], data[ind[j+1]]) < 0) j=j+1;
      }
      if ((* comp) (data[ind[jold]], data[ind[j]]) < 0) {
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

long * heapsort(void **data, long n, void *func) {
  long * ind;
  ind=new long[n];
  heapsort(data, ind, n, func);
  return ind;
}

//uses a binary search to search an ordered list:
long bin_search (void **list, long n, void *value, void *func, long last_ind) {

//+****************************************************************************
//		BIN_SEARCH
//*****************************************************************************
//
// usage:	index=bin_search(list, n, value, func, last_ind)
//
// purpose:	Performs a binary search on an ordered list of values.  Returns
//		the subscript of the smallest array element closest.
//
// parameters:	list	An array of values sorted in ascending order.
//
//		n	Number of elements in the list.
//
//		value	The value to search for.
//
//		func    Comparison function: returns -1 for lt, 0 for eq, 1 for gt
//
//		last_ind	Last value found.  Set to -1 to tell the routine to bracket
//			the entire list right from the start.
//
// author:	Peter Mills (peter.mills@nrl.navy.mil)
//
// history:	First formally documented 2001-12-07 while making some routine
//		improvements.
//		-2002-02-11	PM added the last keyword
//
//		2003-2-25 PM: converted to a C++ template and cleaned up the code a bit
//		2003-3-09 PM: split into two routines;  changed behaviour so
//			that a value higher than the highest returns the
//			largest subscript instead of the number of elements.
//			
//		2003-12-13 PM: corrected bug that made it return a false
//			value if the search value was one step ahead of the last
//			Also made a small change to make it more efficient.
//
//		2004-1-27 PM: cleaned up the logic some more...
//
//		2015-6-23 PM: made this "old-style" version...
//
//-******************************************************************************

  long first, last, mid;		//indices to bracket list
  long step;
  double frac;
  int compres;

  int (*comp) (void *, void *)=(int (*) (void *, void *)) func;

  last=n-1;
  first=0;

  //check to see if the value is either first or last in the list:
  compres=(*comp)(value, list[last]);
  if (compres==0) {
    return last;
  } else if (compres>0) {
    return last;
  }
  compres=(*comp)(value, list[first]);
  if (compres==0) {
    return first;
  } else if (compres<0) {
    return -1;
  }

  //if last is set, search first in the vicinity of that
  //index, expanding outward from there:
  if (last_ind >= 0 && last_ind < n) {
    compres=(*comp)(value, list[last_ind]);
    if (compres==0) {
      return last_ind;
    } else if (compres<0) {
      first=last_ind;
      last=last_ind;
      step=-1;
      do {
        last=first;
        first=first+step;
        if (first<0) {
	  first=0;
	  break;
	}
	compres=(*comp)(value, list[first]);
        if (compres == 0) return first;
        step=step*2;
      } while (value < list[first]);
    } else {
      first=last_ind;
      last=last_ind;
      step=1;
      do {
        first=last;
        last=last+step;
        if (last >= n) {
	  last=n-1;
	  break;
	}
	compres=(*comp)(value, list[last]);
        if (compres == 0) return last;
        step=step*2;
      } while (compres > 0);
    }
  }


  //the actual binary search (not very big is it??)
  while (last-first > 1) {
    mid=(last+first)/2;
//    print, time_string(value), time_string list[mid]), comp
    compres=(*comp)(value, list[mid]);
    if (compres == 0) return mid;
    if (compres > 0) {
      first=mid;
    } else {
      last=mid;
    }
  }

  return first;

}

} //end namespace libpetey

