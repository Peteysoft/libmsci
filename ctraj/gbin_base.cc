#include <math.h>
#include <assert.h>

#include "gbin_base.h"
#include "azeq_binsub.h"
#include "linked.cc"

#define MAX_NEDGE 20

template <class gbin_index, class coord_t, class binsub_t>
gbin_base<gbin_index, coord_t, binsub_t>::gbin_base() {
  binsub=NULL;
  nbin=0;
  bins=NULL;
  nel=NULL;
}

template <class gbin_index, class coord_t, class binsub_t>
gbin_base<gbin_index, coord_t, binsub_t>::~gbin_base() {

}

template <class gbin_index, class coord_t, class binsub_t>
void gbin_base<gbin_index, coord_t, binsub_t>::default_init(void *params) {
  if (params==NULL) binsub=new binsub_t(); else binsub=new binsub_t(params);
  this->nbin=this->binsub->nbin();
  this->bins=new global_bin<gbin_index, coord_t> **[this->nbin];
  this->nel=new int32_t[this->nbin];
  for (int32_t i=0; i<this->nbin; i++) {
    this->bins[i]=NULL;
    this->nel[i]=0;
  }
}

template <class gbin_index, class coord_t, class binsub_t>
void gbin_base<gbin_index, coord_t, binsub_t>::empty() {
  for (int32_t i=0; i<this->nbin; i++) {
    for (int32_t j=0; j<this->nel[i]; j++) delete this->bins[i][j];
    if (this->bins[i]!=NULL) delete [] this->bins[i];
    bins[i]=NULL;
    nel[i]=0;
  }
}

template <class gbin_index, class coord_t, class binsub_t>
void gbin_base<gbin_index, coord_t, binsub_t>::clear() {
  empty();
  if (bins!=NULL) delete [] this->bins;
  if (nel!=NULL) delete [] this->nel;
  if (binsub!=NULL) delete this->binsub;
}

template <class gbin_index, class coord_t, class binsub_t>
float *gbin_base<gbin_index, coord_t, binsub_t>::calc_distances(coord_t &coords, int32_t sub) {
  float *d;

  if (nel[sub]==0) return NULL;
  d=new float[nel[sub]];
  for (int32_t i=0; i<nel[sub]; i++) {
    d[i]=binsub->metric(coords, bins[sub][i]->coords);
  }
  return d;
}

template <class gbin_index, class coord_t, class binsub_t>
float gbin_base<gbin_index, coord_t, binsub_t>::nn1(coord_t &coords, global_bin<gbin_index, coord_t> *&nn) {
  float *d1, dmin, dedge;
  int32_t sub0, sub1;
  int32_t i0;

  dmin=-1;
  sub0=binsub->getsub(coords);
  dedge=binsub->min_edge_distance(coords);
  d1=calc_distances(coords, sub0);

  if (nel[sub0] > 0) {
    dmin=d1[0];
    nn=bins[sub0][0];
    for (int32_t i=1; i<nel[sub0]; i++) {
      if (d1[i] < dmin) {
        dmin=d1[i];
        nn=bins[sub0][i];
      }
    }
    delete [] d1;
  }

  if (dmin < dedge && dmin > 0) return dmin;

  for (int32_t j=0; j<nbin; j++) {
    if (j==sub0) continue;
    d1=calc_distances(coords, j);
    if (d1 == NULL) continue;

    i0=0;
    if (dmin == -1) {
      dmin=d1[0];
      nn=bins[j][0];
      i0=1;
    }
    for (int32_t i=i0; i<nel[j]; i++) {
      if (d1[i] < dmin) {
        dmin=d1[i];
        nn=bins[j][i];
        //printf("dmin=%d; ind=%ld; (%f, %f)\n", dmin, nn->ind, nn->coords.lon, nn->coords.lat);
      }
    }
    delete [] d1;
  }
  return dmin;
}

template <class gbin_index, class coord_t, class binsub_t>
float gbin_base<gbin_index, coord_t, binsub_t>::nn(coord_t &coords, global_bin<gbin_index, coord_t> *&nn) {
  float *d1, d2[MAX_NEDGE], dmin;
  int32_t nedge;
  int32_t sub0, sub1;
  int32_t sub[MAX_NEDGE];
  int32_t k;
  int32_t i0;
  int contflag;

  dmin=-1;
  sub0=binsub->getsub(coords);
  d1=calc_distances(coords, sub0);

  if (nel[sub0] > 0) {
    dmin=d1[0];
    nn=bins[sub0][0];
    for (int32_t i=1; i<nel[sub0]; i++) {
      if (d1[i] < dmin) {
        dmin=d1[i];
        nn=bins[sub0][i];
      }
    }
    delete [] d1;
  }

  nedge=binsub->edge_distances(coords, sub, d2);

/*
  for (int32_t i=0; i<nedge; i++) {
    int32_t xind, yind;
    short hemi;
    binsub->convertsub(sub[i], xind, yind, hemi);
    printf("%d %d %d\n", xind, yind, hemi);
    printf("%f %d %d\n", sqrt(d2[i]), sub[i], nel[sub[i]]);
  }
  printf("%f\n", sqrt(d2[nedge]));
*/

  if (dmin < d2[0] && dmin > 0) return dmin;

  k=0;
  while ((dmin<0 || dmin > d2[k]) && k<nedge) {
    d1=calc_distances(coords, sub[k]);
    if (d1 == NULL) {
      k++;
      continue;
    }
    i0=0;
    if (dmin == -1) {
      dmin=d1[0];
      nn=bins[sub[k]][0];
      i0=1;
    }
    for (int32_t i=i0; i<nel[sub[k]]; i++) {
      if (d1[i] < dmin) {
        dmin=d1[i];
        nn=bins[sub[k]][i];
      }
    }
    delete [] d1;
    k++;
  }

  if (dmin < d2[nedge] && dmin > 0) return dmin;

  assert(k==nedge);

  fprintf(stderr, "global_binning: NN not found in adjactent bins\n");
  fprintf(stderr, "	searching all entries (%f %f)--%d\n", dmin, d2[nedge-1], k);
  for (int32_t j=0; j<nbin; j++) {
    //check to see if we've already calculated the distances:
    if (j==sub0) continue;
    contflag=0;
    for (int32_t k1=0; k1<k; k1++) {
      if (j==sub[k1]) {
        contflag=1;
        break;
      }
    }
    if (contflag) continue;

    d1=calc_distances(coords, j);
    if (d1 == NULL) continue;

    i0=0;
    if (dmin == -1) {
      dmin=d1[0];
      nn=bins[j][0];
      i0=1;
    }
    for (int32_t i=i0; i<nel[j]; i++) {
      if (d1[i] < dmin) {
        dmin=d1[i];
        nn=bins[j][i];
        //printf("dmin=%d; ind=%ld; (%f, %f)\n", dmin, nn->ind, nn->coords.lon, nn->coords.lat);
      }
    }
    delete [] d1;
  }

  return dmin;
}


template class gbin_base<int64_t,lonlat_coord,azeq_binsub>;

