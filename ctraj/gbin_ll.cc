#include <assert.h>

#include "gbin_ll.h"
#include "azeq_binsub.h"
#include "linked.cc"

#define MAX_NEDGE 20

template <class gbin_index, class coord_t, class binsub_t>
gbin_ll<gbin_index, coord_t, binsub_t>::gbin_ll() {
  this->default_init(NULL);
  bininit=new linked_list<global_bin<gbin_index, coord_t> *> *[this->nbin];
  for (int32_t i=0; i<this->nbin; i++) {
    bininit[i]=new linked_list<global_bin<gbin_index, coord_t> *>();
  }
  update_flag=1;
}

template <class gbin_index, class coord_t, class binsub_t>
gbin_ll<gbin_index, coord_t, binsub_t>::gbin_ll(void *params) {
  this->default_init(params);
  bininit=new linked_list<global_bin<gbin_index, coord_t> *> *[this->nbin];
  for (int32_t i=0; i<this->nbin; i++) {
    bininit[i]=new linked_list<global_bin<gbin_index, coord_t> *>();
  }
  update_flag=1;
}

template <class gbin_index, class coord_t, class binsub_t>
gbin_ll<gbin_index, coord_t, binsub_t>::~gbin_ll() {
  this->clear();
  for (int32_t i=0; i<this->nbin; i++) {
    delete bininit[i];
  }
  delete [] bininit;
}

template <class gbin_index, class coord_t, class binsub_t>
int gbin_ll<gbin_index, coord_t, binsub_t>::update() {
  global_bin<gbin_index, coord_t> **newbin;
  long nnew;
  global_bin<gbin_index, coord_t> **newbina;

  if (update_flag) return 0;
  //transfer data from linked list to array:
  for (int32_t i=0; i<this->nbin; i++) {
    newbin=bininit[i]->make_array(nnew);
    if (nnew==0) continue;
    //printf("%d elements found in bin %d\n", nnew, i);
    if (this->nel[i]==0) {
      this->nel[i]=nnew;
      this->bins[i]=newbin;
      continue;
    }
    newbina=new global_bin<gbin_index, coord_t> *[nnew+this->nel[i]];
    for (int32_t j=0; j<this->nel[i]; j++) newbina[j]=this->bins[i][j];
    for (int32_t j=0; j<nnew; j++) newbina[this->nel[i]+j]=newbin[j];
    delete [] this->bins[i];
    delete [] newbin;
    this->bins[i]=newbina;
    this->nel[i]+=nnew;
    bininit[i]->reset();
  }
  update_flag=1;
  return 0;		//should return something useful...
}

template <class gbin_index, class coord_t, class binsub_t>
int gbin_ll<gbin_index, coord_t, binsub_t>::add_el(global_bin<gbin_index, coord_t> *newel) {
  int32_t sub;

  update_flag=0;
  sub=this->binsub->getsub(newel->coords);
  bininit[sub]->add(newel);

  return 0;		//ditto
}


template class gbin_ll<int64_t,lonlat_coord,azeq_binsub>;

