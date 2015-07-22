#include <stdint.h>
#include <assert.h>

//#include "string_operators.h"
#include "string_petey.h"
#include "symbol_table.h"
#include "peteys_tmpl_lib.h"

namespace libpetey {

template <class sym_t>
long * symbol_table<sym_t>::collect_garbage(long &ndup) {
  update();
  long *result;
  long *ind;
  long k;

  ndup=n-nunq;
  result=new long[ndup];
  
  ind=heapsort(sym, n);
  k=0;
  for (long i=1; i<n; i++) {
    if (sym[ind[i-1]]==sym[ind[i]]) {
      if (ind[i-1] > ind[i]) {
        result[k]=ind[i];
        continue; 
      }
    } else {
      result[k]=ind[i-1];
      k++;
    }
  }
  assert(k==ndup);
  delete [] ind;

  return result;
}

template <class sym_t>
void symbol_table<sym_t>::update() {
  long k;
  if (update_flag) return;
  if (sind!=NULL) delete [] sind;
  sind=heapsort(sym, n);
  //search for duplicates (no garbage collection yet):
  nunq=n;
  update_flag=1;
  //return;

  k=0;
  for (long i=1; i<n; i++) {
    if (sym[sind[i-1]]==sym[sind[i]]) {
      if (sind[i-1] > sind[i]) {
        //k++;
        continue; 
      }
    } else {
      k++;
    }
    sind[k]=sind[i];
  }
  nunq=k+1;
}
      
template <class sym_t>
symbol_table<sym_t>::symbol_table() {
  array_size=1;
  n=0;
  nunq=0;
  sym=new sym_t[array_size];
  sind=NULL;
  update_flag=1;
}

template <class sym_t>
symbol_table<sym_t>::~symbol_table() {
  delete [] sym;
  if (sind!=NULL) delete [] sind;
}

template <class sym_t>
long symbol_table<sym_t>::add(sym_t name) {
  sym_t *sym_new;

  //if new symbol doesn't fit, increase size of array:
  if (n+1>=array_size) {
    array_size=array_size*2;
    sym_new=new sym_t[array_size];
    for (long i=0; i<n; i++) {
      sym_new[i]=sym[i];
    }
    delete [] sym;
    sym=sym_new;
  }

  //pretty basic, just paste the value onto the end of the array:
  sym[n]=name;
  n++;
  update_flag=0;

  return n-1;
}

//looks up a symbol in the table and returns a unique id:
template <class sym_t>
long symbol_table<sym_t>::lookup(sym_t name) {
  long low, mid, high;

  update();
  if (n==0) return -1;

  low=0;
  high=nunq-1;

  //binary search:
  if (sym[sind[low]]==name) {
    return sind[low];
  }

  if (sym[sind[low]] > name) {
    return -1;
  }

  if (sym[sind[high]]==name) {
    return sind[high];
  }

  if (sym[sind[high]] < name) {
    return -1;
  }

  while (high-low > 1) {
    mid=(low+high)/2;
    if (sym[sind[mid]] < name) {
      low=mid;
    } else if (sym[sind[mid]] > name) {
      high=mid;
    } else {
      return sind[mid];
    }
  }
  return -1;

}
      
template <class sym_t>
sym_t symbol_table<sym_t>::get(long id) {
  return sym[id];
}

template <class sym_t>
long symbol_table<sym_t>::getid(long ind) {
  update();
  return sind[ind];
}

template <class sym_t>
sym_t symbol_table<sym_t>::let(long ind) {
  update();
  return sym[sind[ind]];
}

template <class sym_t>
long symbol_table<sym_t>::entries(int uflag) {
  if (uflag) {
    update();
    return nunq; 
  }
  return n;
}

template <>
void symbol_table<string_petey>::print() {
  update();
  for (long i=0; i<n; i++) {
    printf("%ld ", i);
    sym[i].print();
  }
  for (long i=0; i<nunq; i++) {
    printf("%ld ", sind[i]);
    sym[sind[i]].print();
  }
}
  

//template class symbol_table<char *>;
template class symbol_table<string_petey>;
template class symbol_table<float>;
template class symbol_table<double>;
template class symbol_table<int32_t>;
template class symbol_table<int64_t>;
//template class symbol_table<char>;

} //end namespace libpetey

