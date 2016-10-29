#include <stdint.h>
#include <stdio.h>
#include <algorithm>

#include "randomize.h"
#include "quicksort.h"
#include "kselect.h"
#include "kextreme.h"

using namespace std;

namespace libpetey {

  //lets have a showdown:
  template <class type>
  kselect_base<type>::kselect_base(long k1) {
    k=k1;
  }

  template <class type>
  kselect_base<type>::kselect_base() {}

  template <class type>
  kselect_base<type>::~kselect_base() {}

  template <class type>
  kselect_naive<type>::kselect_naive(long k1) {
    select=new type[k1];
    this->k=k1;
    ncur=0;
  }

  template <class type>
  kselect_naive<type>::~kselect_naive() {
    delete [] select;
  }

  template <class type>
  long kselect_naive<type>::add(type val) {
    long gind=-1;
    if (ncur<this->k) {
      select[ncur]=val;
      ncur++;
    } else {
      for (long i=0; i<this->k; i++) {
        if (gind!=-1) {
          if (select[i]>select[gind]) gind=i;
	} else if (val<select[i]) {
	  gind=i;
        }
      }
      if (gind!=-1) select[gind]=val;
    }
    return ncur;
  }

  template <class type>
  void kselect_naive<type>::get(type *kleast) {
    for (long i=0; i<this->k; i++) kleast[i]=select[i];
  }

  template <class type>
  kselect_tree<type>::kselect_tree(long k1) {
    this->k=k1;
  }

  template <class type>
  kselect_tree<type>::~kselect_tree() {
  }

  template <class type>
  long kselect_tree<type>::add(type val) {
    long n;
    n=data.add(val);
    if (n>this->k) {
      n=data.delete_greatest();
    }
    return n;
  }

  template <class type>
  void kselect_tree<type>::get(type *kleast) {
    data.decompose(kleast, this->k);
  }

  template <class type>
  kselect_heap<type>::kselect_heap(long k1) {
    this->k=k1;
    data=new vector<type>(0);
  }

  template <class type>
  kselect_heap<type>::~kselect_heap() {
    delete data;
  }

  template <class type>
  long kselect_heap<type>::add(type val) {
    long n;
    n=data->size();
    data->push_back(val);
    if (n>this->k) {
      (*data)[n]=val;
      push_heap(data->begin(), data->end());
      pop_heap(data->begin(), data->end());
      data->pop_back();
    } else if (n==this->k) {
      make_heap(data->begin(), data->end());
      data->push_back(val);
    }
    return n;
  }

  template <class type>
  void kselect_heap<type>::get(type *kleast) {
    sort_heap(data->begin(), data->end());
    for (long i=0; i<this->k; i++) kleast[i]=(*data)[i];
  }

  template <class type>
  kselect_quick<type>::kselect_quick(long k1) {
    this->k=k1;
    data=new vector<type>(0);
    ncur=0;
  }

  template <class type>
  kselect_quick<type>::kselect_quick(long k1, long n) {
    this->k=k1;
    //since I can't figure out how to (easily) assign to a vector, extending
    //it if necessaray, the n parameter is meaningless:
    data=new vector<type>(0);
    ncur=0;
  }

  template <class type>
  kselect_quick<type>::~kselect_quick() {
    delete data;
  }

  template <class type>
  long kselect_quick<type>::add(type val) {
    data->push_back(val);
    ncur++;
    //printf("kselect_quick: adding value; current list:\n");
    //data->print(stdout);
    return ncur;
  }

  template <class type>
  void kselect_quick<type>::get(type *kleast) {
    type *d2;
    long n;
    //either one will work, but the second is less complex:
    //d2=&((*data)[0]);
    d2=data->data();
    n=data->size();
    for (long i=0; i<n; i++) printf("%g ", d2[i]);
    printf("\n");
    kleast_quick(d2, n, this->k);
    for (long i=0; i<this->k; i++) kleast[i]=d2[i];
  }

  template <class type>
  kiselect_base<type>::kiselect_base(long k1) {
    k=k1;
    ncur=0;
  }

  template <class type>
  kiselect_base<type>::kiselect_base() {}

  template <class type>
  kiselect_base<type>::~kiselect_base() {}

  //trying to move towards more of a test-driven development:
  template <class type>
  int kiselect_base<type>::test(long n) {
    int err=0;
    type list[n];
    type kleast[k];
    long ind[k];
    long lind;			//index of largest of k-least
    int flag[n];
   
    //generate a list of random numbers and apply k-least algo to it: 
    for (long i=0; i<n; i++) {
      list[i]=ranu();
      add(list[i]);
    }
    get(kleast, ind);

    //find the largest of the k-least:
    lind=0;
    for (long i=1; i<k; i++) {
      if (kleast[i]>kleast[lind]) lind=i;
    }

    //set flags to exclude all k-least from the comparison:
    for (long i=0; i<n; i++) flag[i]=1;
    for (long i=0; i<k; i++) flag[ind[i]]=0;
    
    //largest of k-least must be smaller than all others in the list:
    for (long i=0; i<n; i++) {
      if (flag[i] && kleast[lind]>list[i]) {
        err=-1;
        break;
      }
      if (err!=0) break;
    }
    int err2=verify_kleast(list, n, kleast, k);
    if ((err2==0) && (err!=0) || (err2==1 && err!=-1)) {
      fprintf(stderr, "Two error states differ: %d %d\n", err, err2);
      err=1;
    }
    if (err!=0) {
      for (long i=0; i<n; i++) printf("%g ", list[i]);
      printf("\n");
      for (long i=0; i<k; i++) printf("%g ", kleast[i]);
      printf("\n");
    }
    return err;
  }

  template <class type>
  kiselect_tree<type>::kiselect_tree(long k1) {
    this->k=k1;
    this->ncur=0;
  }

  template <class type>
  kiselect_tree<type>::~kiselect_tree() {
  }

  template <class type>
  kiselect_naive<type>::kiselect_naive(long k1) {
    select=new type[k1];
    ind=new long [k1];
    this->k=k1;
    this->ncur=0;
  }

  template <class type>
  kiselect_naive<type>::~kiselect_naive() {
    delete [] select;
    delete [] ind;
  }

  template <class type>
  long kiselect_naive<type>::add(type val) {
    long gind=-1;
    if (this->ncur<this->k) {
      select[this->ncur]=val;
      ind[this->ncur]=this->ncur;
    } else {
      for (long i=0; i<this->k; i++) {
        if (gind!=-1) {
          if (select[i]>select[gind]) gind=i;
	} else if (val<select[i]) {
	  gind=i;
        }
      }
      if (gind!=-1) {
        select[gind]=val;
        ind[gind]=this->ncur;
      }
    }
    this->ncur++;
    return this->ncur;
  }

  template <class type>
  long kiselect_naive<type>::add(type val, long i) {
    if (this->ncur<this->k) {
      select[this->ncur]=val;
      ind[this->ncur]=i;
    } else {
      for (long j=0; j<this->k; j++) {
        if (val<select[i]) {
          select[j]=val;
          ind[j]=i;
          break;
        }
      }
    }
    this->ncur++;
    return this->ncur;
  }

  template <class type>
  void kiselect_naive<type>::get(type *kleast, long *j) {
    for (long i=0; i<this->k; i++) {
      kleast[i]=select[i];
      j[i]=ind[i];
    }
  }

  template <class type>
  long kiselect_tree<type>::add(type val) {
    long n;
    n=data.add(val, this->ncur);
    this->ncur++;
    if (n>this->k) {
      n=data.delete_greatest();
    }
    return n;
  }

  template <class type>
  long kiselect_tree<type>::add(type val, long ind) {
    long n;
    n=data.add(val, ind);
    this->ncur++;
    if (n>this->k) {
      n=data.delete_greatest();
    }
    return n;
  }

  template <class type>
  void kiselect_tree<type>::get(type *kleast, long *ind) {
    data.decompose(kleast, ind, this->k);
  }

  template <class type>
  int heap_ind_el_comp(const heap_ind_el<type> &v1, heap_ind_el<type> &v2) {
    if (v1.data<v2.data) return 1;
    return 0;
  }

  template <class type>
  kiselect_heap<type>::kiselect_heap(long k1) {
    this->k=k1;
    data=new vector<heap_ind_el<type> >(0);
    this->ncur=0;
  }

  template <class type>
  kiselect_heap<type>::~kiselect_heap() {
    delete data;
  }

  template <class type>
  long kiselect_heap<type>::add(type val) {
    long n;
    heap_ind_el<type> v1;
    v1.data=val;
    v1.ind=this->ncur;
    this->ncur++;
    data->push_back(v1);
    n=data->size();
    if (n>this->k) {
      push_heap(data->begin(), data->end(), &heap_ind_el_comp<type>);
      pop_heap(data->begin(), data->end(), &heap_ind_el_comp<type>);
      data->pop_back();
    } else if (n==this->k) {
      make_heap(data->begin(), data->end(), &heap_ind_el_comp<type>);
    }
    return n;
  }

  template <class type>
  long kiselect_heap<type>::add(type val, long ind) {
    long n;
    heap_ind_el<type> v1;
    v1.data=val;
    v1.ind=ind;
    this->ncur++;
    data->push_back(v1);
    n=data->size();
    if (n>this->k) {
      push_heap(data->begin(), data->end(), &heap_ind_el_comp<type>);
      pop_heap(data->begin(), data->end(), &heap_ind_el_comp<type>);
      data->pop_back();
    } else if (n==this->k) {
      make_heap(data->begin(), data->end(), &heap_ind_el_comp<type>);
    }
    return n;
  }

  template <class type>
  void kiselect_heap<type>::get(type *kleast, long *ind) {
    sort_heap(data->begin(), data->end(), &heap_ind_el_comp<type>);
    for (long i=0; i<this->k; i++) {
      kleast[i]=(*data)[i].data;
      ind[i]=(*data)[i].ind;
    }
  }

  template <class type>
  kiselect_quick<type>::kiselect_quick(long k1) {
    this->k=k1;
    data=new vector<type>(0);
    ind=new vector<long>(0);
    this->ncur=0;
  }

  template <class type>
  kiselect_quick<type>::kiselect_quick(long k1, long n) {
    this->k=k1;
    data=new vector<type>(0);
    ind=new vector<long>(0);
    this->ncur=0;
  }

  template <class type>
  kiselect_quick<type>::~kiselect_quick() {
    delete data;
    delete ind;
  }

  template <class type>
  long kiselect_quick<type>::add(type val) {
    data->push_back(val);
    ind->push_back(this->ncur);
    this->ncur++;
    return data->size();
  }

  template <class type>
  long kiselect_quick<type>::add(type val, long ind2) {
    data->push_back(val);
    ind->push_back(ind2);
    this->ncur++;
    return data->size();
  }

  template <class type>
  void kiselect_quick<type>::get(type *kleast, long *ind2) {
    long n;
    type *d2;
    long *ind1;
    n=data->size();
    d2=data->data();
    ind1=ind->data();
    n=data->size();
    kleast_quick(d2, n, this->k, ind1, 0, n-1);
    for (long i=0; i<this->k; i++) {
      //kleast[i]=d2[i];
      kleast[i]=d2[ind1[i]];
      ind2[i]=ind1[i];
    }
  }

  template class kselect_naive<float>;
  template class kselect_tree<float>;
  template class kselect_heap<float>;
  template class kselect_quick<float>;

  template class kselect_naive<double>;
  template class kselect_tree<double>;
  template class kselect_heap<double>;
  template class kselect_quick<double>;

  template class kiselect_base<float>;
  template class kiselect_naive<float>;
  template class kiselect_tree<float>;
  template class kiselect_heap<float>;
  template class kiselect_quick<float>;

  template class kiselect_base<double>;
  template class kiselect_naive<double>;
  template class kiselect_tree<double>;
  template class kiselect_heap<double>;
  template class kiselect_quick<double>;

} //end namespace libpetey
