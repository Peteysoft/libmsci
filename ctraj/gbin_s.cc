#include <math.h>
#include <assert.h>

#include "peteys_tmpl_lib.h"
#include "randomize.h"
#include "gbin_s.h"
#include "azeq_binsub.h"
#include "linked.cc"

#define MAX_NEDGE 20

template <class gbin_index, class coord_t, class binsub_t>
gbin_s<gbin_index, coord_t, binsub_t>::gbin_s() {
  this->default_init(NULL);
}

template <class gbin_index, class coord_t, class binsub_t>
gbin_s<gbin_index, coord_t, binsub_t>::gbin_s(void *params) {
  this->default_init(params);
}

template <class gbin_index, class coord_t, class binsub_t>
gbin_s<gbin_index, coord_t, binsub_t>::~gbin_s() {
  this->clear();
}

template <class gbin_index, class coord_t, class binsub_t>
gbin_s<gbin_index, coord_t, binsub_t>::gbin_s(coord_t *v, int32_t n, int32_t nsamp) {
  float d[nsamp];
  float dcalc, dmin;
  float dave, dstd;
  int32_t k1, k2;
  int32_t binpside;

  //figure out the optimal bin size:
  //(interesting, in order to figure out the optimal bin size,
  //we have to perform a nearest-neighbours search...)
  //oh well, lets run with it...
  for (int32_t i=0; i<nsamp; i++) {
    k1=n*ranu();
    dmin=-1;
    for (int32_t j=0; j<n; j++) { 
      //printf("  %d\n", j);
      if (j==k1) continue;
      dcalc=this->binsub->metric(v[k1], v[j]);
      if (dmin==-1 || dcalc<dmin) {
        dmin=dcalc;
      }
    }
    d[i]=sqrt(dmin);
    //printf("%d; %f\n", i, d[i]);
  }

  dave=0;
  for (int32_t i=0; i<nsamp; i++) {
    dave+=d[i];
  }
  dave=dave/nsamp;

  dstd=0;
  for (int32_t i=0; i<nsamp; i++) {
    dstd+=(d[i]-dave)*(d[i]-dave);
  }
  dstd=sqrt(dstd/(nsamp-1));

  binpside=20000/(dave+2*dstd);
  printf("Using a bin size of %f; %d bins per side\n", dave+2*dstd, binpside);

  this->binsub=new binsub_t(&binpside);
  this->nbin=this->binsub->nbin();
  this->bins=new global_bin<gbin_index, coord_t> **[this->nbin];
  this->nel=new int32_t[this->nbin];
  for (int32_t i=0; i<this->nbin; i++) {
    this->bins[i]=NULL;
    this->nel[i]=0;
  }

  init(v, n); 
}

template <class gbin_index, class coord_t, class binsub_t>
void gbin_s<gbin_index, coord_t, binsub_t>::init(coord_t *v, int32_t n) {
  int32_t sub[n];
  long *sind;
  int32_t k;
  global_bin<gbin_index, coord_t> *newel;

  this->empty();

  for (int32_t i=0; i<n; i++) sub[i]=this->binsub->getsub(v[i]);
  sind=heapsort(sub, n);
  map_vector_inplace(sub, sind, n);

  k=0;
  for (int32_t i=1; i<n; i++) {
    if (sub[i]!=sub[i-1]) {
      this->nel[sub[i-1]]=i-k;
      this->bins[sub[i-1]]=new global_bin<gbin_index, coord_t> *[i-k];
      for (int32_t j=0; j<this->nel[sub[i-1]]; j++) {
        newel=new global_bin<gbin_index, coord_t>;
        newel->coords=v[sind[k+j]];
        newel->ind=sind[k+j];
        this->bins[sub[i-1]][j]=newel;
      }
      k=i;
    }
    this->nel[sub[i-1]]=n-k+1;
    this->bins[sub[i-1]]=new global_bin<gbin_index, coord_t> *[n-k+1];
    for (int32_t j=0; j<this->nel[sub[i-1]]; j++) {
      newel=new global_bin<gbin_index, coord_t>;
      newel->coords=v[sind[k+j]];
      newel->ind=sind[k+j];
      this->bins[sub[i-1]][j]=newel;
    }
  }

}


template class gbin_s<int64_t,lonlat_coord,azeq_binsub>;

