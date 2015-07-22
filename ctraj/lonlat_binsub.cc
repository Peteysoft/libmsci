#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "lonlat_binsub.h"

#include "global_metric.h"
#include "peteys_tmpl_lib.h"

lonlat_binsub::lonlat_binsub() {
  nlat=17;
  offset=5;

  nlon=new int32_t[nlat+2];
  nlon[0]=1;
  nlon[1]=9;
  nlon[2]=15;
  nlon[3]=21;
  nlon[4]=26;
  nlon[5]=30;
  nlon[6]=33;
  nlon[7]=35;
  nlon[8]=36;
  nlon[9]=36;
  nlon[10]=36;
  nlon[11]=35;
  nlon[12]=33;
  nlon[13]=30;
  nlon[14]=26;
  nlon[15]=21;
  nlon[16]=15;
  nlon[17]=9;
  nlon[18]=1;

  cumulate_nlon();
}

lonlat_binsub::lonlat_binsub(void *params) {
  FILE *fs;
  float lat;
  float checklat;
  long err=0;
  char *fname;

  fname=(char *) params[0];

  fs=fopen(fname, "r");
  fscanf(fs, "%ld", &nlat);
  nlon=new int32_t[nlat+2];
  fscanf(fs, "%f %ld", &lat, nlon+1);
  offset=90+lat;
  for (long i=1; i<nlat; i++) {
    checklat=i*(180-offset*2)/nlat-90+offset;
    fscanf(fs, "%f %d", &lat, nlon+i+1);
    if (lat!=checklat) {
      fprintf(stderr, "Scanned and calculated latitudes do not agree: %f vs. %f\n", lat, checklat);
      err=1;
    }
  }

  nlon[0]=1;
  nlon[nlat+1]=1;

  fclose(fs);

  cumulate_nlon();
  //return err;

}

lonlat_binsub::~lonlat_binsub() {
  delete [] nlon;
  delete [] nlon_cum;
}

void lonlat_binsub::cumulate_nlon() {
  nlon_cum=new int32_t[nlat+3];
  nlon_cum[0]=0;
  for (int32_t i=1; i<nlat+3; i++) {
    nlon_cum[i]=nlon_cum[i-1]+nlon[i-1];
  }
}

void lonlat_binsub::getsub(float lon, float lat, int32_t &lonind, int32_t &latind) {
  latind=(lat+90.-offset)/(180.-2*offset)*nlat+1;
  lonind=(lon/360.*nlon[latind]);
}

int32_t lonlat_binsub::convert_sub(int32_t lonind, int32_t latind) {
  return nlon_cum[latind]+lonind;
}

int32_t lonlat_binsub::getsub(const lonlat_coord &coord) {
  int32_t lonind, latind;
  int32_t sub;
  getsub(coord.lon, coord.lat, lonind, latind);
  sub=convert_sub(lonind, latind);
  //printf("nbins=%d; binsub=%d\n", nlon_cum[nlat+2], sub);
  return sub;
}

float lonlat_binsub::metric(const lonlat_coord &c1, const lonlat_coord &c2) {
  return sdist(c1.lon, c1.lat, c2.lon, c2.lat);
}

int32_t lonlat_binsub::edge_distances(const lonlat_coord &coord, int32_t *sub, float *d) {

  int32_t lonind0, latind0;
  int32_t lonind, latind;
  int32_t n;
  int32_t l1, l2;
  long *ind;
  float lon, lat;

  n=0;
  if (nlon[latind0] > 1) {
    //bin to right:
    d[n]=sdist(coord.lon, coord.lat, 360*(lonind0+1)/nlon[latind0], coord.lat);
    lonind=lonind0+1;
    if (lonind0 >= nlon[latind0]) lonind0=0;
    latind=latind0;
    sub[n]=convert_sub(lonind, latind);
    printf("(%d, %d): %f\n", lonind, latind, d[n]);
    n++;

    if (nlon[latind0] > 2) {
      //bin to left:
      getsub(coord.lon, coord.lat, lonind0, latind0);
      d[n]=sdist(coord.lon, coord.lat, 360*lonind0/nlon[latind0], coord.lat);
      lonind=lonind0-1;
      if (lonind0 < 0) lonind0=nlon[latind0]-1;
      latind=latind0;
      sub[n]=convert_sub(lonind, latind);
      printf("(%d, %d): %f\n", lonind, latind, d[n]);
      n++;
    }
  }

  if (latind0 > 0) {
    //bin directly below:
    lat=lat_from_sub(latind0);
    d[n]=sdist(coord.lon, coord.lat, coord.lon, lat);
    lat=lat_from_sub(latind0-0.5);
    getsub(coord.lon, lat, lonind, latind);
    sub[n]=convert_sub(lonind, latind);
    printf("(%d, %d): %f\n", lonind, latind, d[n]);
    n++;

    //diagonally adjacent bins:
    l1=ceil(1.*nlon[latind0-1]*(lonind0)/(nlon[latind0]));
    l2=1.*nlon[latind0-1]*(lonind0+1)/nlon[latind0];
    printf("l1=%d; l2=%d\n", l1, l2);
    for (int32_t i=l1; i<=l2; i++) {
      lon=360.*i/nlon[latind0-1];
      lat=lat_from_sub(latind0);
      d[n]=sdist(coord.lon, coord.lat, lon, lat);
      if (coord.lon > 360*i/nlon[latind0-1]) {
        lonind=i-1;
      } else {
        lonind=i;
      }
      latind=latind0-1;
      sub[n]=convert_sub(lonind, latind);
      printf("(%d, %d): %f\n", lonind, latind, d[n]);
      n++;
    }
  }

  if (latind0 <= nlat) {
    //bin directly above:
    lat=lat_from_sub(latind0+1);
    d[n]=sdist(coord.lon, coord.lat, coord.lon, lat);
    lat=lat_from_sub(latind0+1.5);
    getsub(coord.lon, lat, lonind, latind);
    sub[n]=convert_sub(lonind, latind);
    printf("(%d, %d): %f\n", lonind, latind, d[n]);
    n++;

    //diagonally adjacent bins:
    l1=ceil(1.*nlon[latind0+1]*lonind0/nlon[latind0]);
    l2=1.*nlon[latind0+1]*(lonind0+1)/nlon[latind0];
    printf("l1=%d; l2=%d\n", l1, l2);
    for (int32_t i=l1; i<=l2; i++) {
      lon=360.*i/nlon[latind0+1];
      lat=lat_from_sub(latind0+1);
      d[n]=sdist(coord.lon, coord.lat, lon, lat);
      if (coord.lon > 360*i/nlon[latind0+1]) {
        lonind=i-1;
      } else {
        lonind=i;
      }
      latind=latind0+1;
      sub[n]=convert_sub(lonind, latind);
      printf("(%d, %d): %f\n", lonind, latind, d[n]);
      n++;
    }
  }

  //sort the distances:
  ind=new long[n];
  heapsort(d, ind, n);
  map_vector_inplace(d, ind, n);
  map_vector_inplace(sub, ind, n);
  delete [] ind;

  return n;
}

int32_t lonlat_binsub::nbin() {
  return nlon_cum[nlat+2];
}

