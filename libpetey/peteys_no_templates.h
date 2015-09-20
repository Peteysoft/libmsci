#ifndef NO_TMPL_LIB_INCLUDED
#define NO_TMPL_LIB_INCLUDED

namespace libpetey {

  //templates the old-fashioned way:
  void heapsort_inplace(void **data, 		//list of data elements
		  long n, 			//number of elements
		  void *func);			//comparison function
  void heapsort(void **data, long *ind, long n, void *func);
  long *heapsort(void **data, long n, void *func);
  long bin_search(void **data, 			//sorted list
		  long n, 			//number of elements
		  void *val, 			//value to find
		  void *func, 			//comparision function
		  long lastind=-1);		//last element found

}

#endif
