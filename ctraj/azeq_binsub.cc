#include <stdio.h>
#include <stdint.h>
#include <math.h>

#include "azeq_binsub.h"
#include "peteys_tmpl_lib.h"
#include "tcoord_defs.h"

#include "global_metric.h"

azeq_binsub::azeq_binsub() {
  nx=20;
  ny=20;
  sidelength=20000.;
  xoffset=-10000.;
  yoffset=-10000.;
}

azeq_binsub::azeq_binsub(void *params) {
  nx=*(long *) params;
  ny=nx;
  sidelength=20000.;
  xoffset=-10000.;
  yoffset=-10000.;
}

azeq_binsub::~azeq_binsub() {
}

int32_t azeq_binsub::nbin() {
  return 2*nx*ny;
}

int32_t azeq_binsub::convertsub(int32_t xind, int32_t yind, short hemi) {
  return yind*nx+xind+(hemi+1)*nx*ny/2;
}

void azeq_binsub::convertsub(int32_t sub, int32_t &xind, int32_t &yind, short &hemi) {
  if (sub < nx*ny) hemi=-1; else hemi=1;
  sub=sub % (nx*ny);
  xind=sub % nx;
  yind=sub/nx;
}

void azeq_binsub::getsub(float x, float y, int32_t &xind, int32_t &yind) {
  xind=nx*(x-xoffset)/sidelength;
  yind=ny*(y-yoffset)/sidelength;
}

float azeq_binsub::xsub2x(float xind) {
  return xoffset+sidelength*xind/nx;
}

float azeq_binsub::ysub2y(float yind) {
  return yoffset+sidelength*yind/ny;
}

float azeq_binsub::metric(lonlat_coord &c1, lonlat_coord &c2) {
  return sdist (c1.lon, c1.lat, c2.lon, c2.lat);
}

int32_t azeq_binsub::getsub(lonlat_coord &coord) {
  float x, y;
  short hemi;
  int32_t xind, yind;

  hemi=coord.lat/fabs(coord.lat);
  tcoord2_2lonlat(x, y, -1, hemi, coord.lon, coord.lat);

  getsub(x, y, xind, yind);

  return convertsub(xind, yind, hemi);
}

float azeq_binsub::min_edge_distance(lonlat_coord &coord) {

  float x0, y0;
  float x, y;
  float d, dthresh;
  int32_t xind0, yind0;

  //right:
  x=xsub2x(xind0+1);
  y=y0;
  d=tcoord_ds2(x0, y0, x, y);
  dthresh=d;

  //left:
  x=xsub2x(xind0);
  y=y0;
  d=tcoord_ds2(x0, y0, x, y);
  if (d < dthresh) dthresh=d;

  //up:
  x=x0;
  y=ysub2y(yind0+1);
  d=tcoord_ds2(x0, y0, x, y);
  if (d < dthresh) dthresh=d;

  //down:
  x=x0;
  y=ysub2y(yind0);
  d=tcoord_ds2(x0, y0, x, y);
  if (d < dthresh) dthresh=d;

  return dthresh;
}

int32_t azeq_binsub::edge_distances(lonlat_coord &coord,
		int32_t *sub, float *d) {
  float x0, y0;
  float x, y;
  int32_t n, n1;
  int32_t xind0, yind0;
  int32_t xind1, yind1;
  int32_t xind[20], yind[20];
  short hemi;
  float dthresh;	//as long as dmin is below this, it's contained
			//within the cluster of neighbouring bins
  int32_t ithresh;	//as long as dmin is below this, it's contained

  float dx, dy;

  long *ind;

  hemi=coord.lat/fabs(coord.lat);
  tcoord2_2lonlat(x0, y0, -1, hemi, coord.lon, coord.lat);

  getsub(x0, y0, xind0, yind0);
  //printf("%d %d %d %d\n", convertsub(xind0, yind0, hemi), xind0, yind0, hemi);
  n=0;

  //not exact, but then again, neither are our distance calculations:
  //(keeps the code simple...)
  if (xind0==nx-1) xind1=xind0-1; else xind1=xind0;
  if (yind0==ny-1) yind1=yind0-1; else yind1=yind0;
  dx=sqrt(tcoord_ds2(xsub2x(xind1), y0, xsub2x(xind1+1), y0));
  dy=sqrt(tcoord_ds2(x0, ysub2y(yind1), x0, ysub2y(yind1+1)));
  printf("dx=%f; dy=%f\n", dx, dy);

  //right:
  xind[n]=xind0+1;
  yind[n]=yind0;

  dthresh=2*sqrt(dx*dx+dy*dy);

  if (xind[n] < nx) {
    x=xsub2x(xind[n]);
    y=y0;
    d[n]=tcoord_ds2(x0, y0, x, y);
    sub[n]=convertsub(xind[n], yind[n], hemi);
    dthresh=sqrt(d[n])+dx;
    n++;
  }

  //left:
  xind[n]=xind0-1;
  yind[n]=yind0;

  if (xind[n]>=0) {
    x=xsub2x(xind0);
    y=y0;
    d[n]=tcoord_ds2(x0, y0, x, y);
    sub[n]=convertsub(xind[n], yind[n], hemi);
    if (sqrt(d[n])+dx < dthresh) dthresh=sqrt(d[n])+dx;
    n++;
  }

  //up:
  xind[n]=xind0;
  yind[n]=yind0+1;

  if (yind[n] < ny) {
    x=x0;
    y=ysub2y(yind[n]);
    d[n]=tcoord_ds2(x0, y0, x, y);
    sub[n]=convertsub(xind[n], yind[n], hemi);
    if (sqrt(d[n])+dy < dthresh) dthresh=sqrt(d[n])+dy;
    n++;
  }

  //down:
  xind[n]=xind0;
  yind[n]=yind0-1;

  if (yind[n]>=0) {
    x=x0;
    y=ysub2y(yind0);
    d[n]=tcoord_ds2(x0, y0, x, y);
    sub[n]=convertsub(xind[n], yind[n], hemi);
    if (sqrt(d[n])+dy < dthresh) dthresh=sqrt(d[n])+dy;
    n++;
  }

  //diagonals:
  xind[n]=xind0+1;
  yind[n]=yind0+1;

  if (xind[n] < nx && yind[n] < ny) {
    x=xsub2x(xind[n]);
    y=ysub2y(yind[n]);
    d[n]=tcoord_ds2(x0, y0, x, y);
    sub[n]=convertsub(xind[n], yind[n], hemi);
    n++;
  }

  //diagonals:
  xind[n]=xind0-1;
  yind[n]=yind0+1;

  if (xind[n] >=0 && yind[n] < ny) {
    x=xsub2x(xind0);
    y=ysub2y(yind[n]);
    d[n]=tcoord_ds2(x0, y0, x, y);
    sub[n]=convertsub(xind[n], yind[n], hemi);
    n++;
  }

  //diagonals:
  xind[n]=xind0+1;
  yind[n]=yind0-1;

  if (xind[n] < nx && yind[n] >= 0) {
    x=xsub2x(xind[n]);
    y=ysub2y(yind0);
    d[n]=tcoord_ds2(x0, y0, x, y);
    sub[n]=convertsub(xind[n], yind[n], hemi);
    n++;
  }

  //diagonals:
  xind[n]=xind0-1;
  yind[n]=yind0-1;

  if (xind[n] >=0 && yind[n] >= 0) {
    x=xsub2x(xind0);
    y=ysub2y(yind0);
    d[n]=tcoord_ds2(x0, y0, x, y);
    sub[n]=convertsub(xind[n], yind[n], hemi);
    n++;
  }

  //if the bin is close to the equator, we have to check
  //the bins on the other side of the equator:
  //(bit wastefull--most of the bins will be empty...)
  xind1=xind0+x0/fabs(x0);
  yind1=yind0+y0/fabs(y0);

  x=xsub2x(xind1);
  y=ysub2y(yind1);

  if ((x*x+y*y) > 1.e8) {
    hemi=-hemi;
    tcoord2_2lonlat(x0, y0, -1, hemi, coord.lon, coord.lat);
    getsub(x0, y0, xind0, yind0);

    //right:
    xind[n]=xind0+1;
    yind[n]=yind0;

    if (xind[n] < nx && xind[n] >= 0 && yind[n] < ny && yind[n] >= 0) {
      x=xsub2x(xind[n]);
      y=y0;
      d[n]=tcoord_ds2(x0, y0, x, y);
      sub[n]=convertsub(xind[n], yind[n], hemi);
      if (sqrt(d[n])+dx < dthresh) dthresh=sqrt(d[n])+dx;
      n++;
    }

    //left:
    xind[n]=xind0-1;
    yind[n]=yind0;

    if (xind[n] < nx && xind[n] >= 0 && yind[n] < ny && yind[n] >= 0) {
      x=xsub2x(xind0);
      y=y0;
      d[n]=tcoord_ds2(x0, y0, x, y);
      sub[n]=convertsub(xind[n], yind[n], hemi);
      if (sqrt(d[n])+dx < dthresh) dthresh=sqrt(d[n])+dx;
      n++;
    }

    //up:
    xind[n]=xind0;
    yind[n]=yind0+1;

    if (xind[n] < nx && xind[n] >= 0 && yind[n] < ny && yind[n] >= 0) {
      x=x0;
      y=ysub2y(yind[n]);
      d[n]=tcoord_ds2(x0, y0, x, y);
      sub[n]=convertsub(xind[n], yind[n], hemi);
      if (sqrt(d[n])+dy < dthresh) dthresh=sqrt(d[n])+dy;
      n++;
    }

    //down:
    xind[n]=xind0;
    yind[n]=yind0-1;

    if (xind[n] < nx && xind[n] >= 0 && yind[n] < ny && yind[n] >= 0) {
      x=x0;
      y=ysub2y(yind0);
      d[n]=tcoord_ds2(x0, y0, x, y);
      sub[n]=convertsub(xind[n], yind[n], hemi);
      if (sqrt(d[n])+dy < dthresh) dthresh=sqrt(d[n])+dy;
      n++;
    }

    //diagonals:
    xind[n]=xind0+1;
    yind[n]=yind0+1;

    if (xind[n] < nx && xind[n] >= 0 && yind[n] < ny && yind[n] >= 0) {
      x=xsub2x(xind[n]);
      y=ysub2y(yind[n]);
      d[n]=tcoord_ds2(x0, y0, x, y);
      sub[n]=convertsub(xind[n], yind[n], hemi);
      n++;
    }

    //diagonals:
    xind[n]=xind0-1;
    yind[n]=yind0+1;

    if (xind[n] < nx && xind[n] >= 0 && yind[n] < ny && yind[n] >= 0) {
      x=xsub2x(xind0);
      y=ysub2y(yind[n]);
      d[n]=tcoord_ds2(x0, y0, x, y);
      sub[n]=convertsub(xind[n], yind[n], hemi);
      n++;
    }

    //diagonals:
    xind[n]=xind0+1;
    yind[n]=yind0-1;

    if (xind[n] < nx && xind[n] >= 0 && yind[n] < ny && yind[n] >= 0) {
      x=xsub2x(xind[n]);
      y=ysub2y(yind0);
      d[n]=tcoord_ds2(x0, y0, x, y);
      sub[n]=convertsub(xind[n], yind[n], hemi);
      n++;
    }

    //diagonals:
    xind[n]=xind0-1;
    yind[n]=yind0-1;

    if (xind[n] < nx && xind[n] >= 0 && yind[n] < ny && yind[n] >= 0) {
      x=xsub2x(xind0);
      y=ysub2y(yind0);
      d[n]=tcoord_ds2(x0, y0, x, y);
      sub[n]=convertsub(xind[n], yind[n], hemi);
      n++;
    }
  }

  //sort the distances:
  ind=new long[n];
  heapsort(d, ind, n);
  map_vector_inplace(d, ind, n);
  map_vector_inplace(sub, ind, n);

  //print out all the data in one big schwab:
/*
  printf("%d %d %d %d\n", convertsub(xind0, yind0, hemi), xind0, yind0, hemi);
  for (int32_t i=0; i<n; i++) {
    printf("%d %d %d  %f\n", sub[i], xind[ind[i]], yind[ind[i]], sqrt(d[i]));
  }
  printf("%f\n", dthresh); 
*/

  delete [] ind;

  n1=n;
  dthresh=dthresh*dthresh;
  for (int32_t i=0; i<n; i++) {
    if (d[i] > dthresh) {
      n1=i;
      break;
    }
  }

  d[n1]=dthresh;

  return n1;
}

