#include <stdio.h>
#include <math.h>
#include <assert.h>

#include "time_class.h"
#include "randomize.h"
#include "linked.h"
#include "bit_array.h"

#include "meas_data.h"
#include "tcoord_defs.h"
#include "global_metric.h"

#define MAXLL 500

using namespace libpetey;

namespace ctraj {

void interpol_tween_coords(float lon1, float lat1, 
		float lon2, float lat2,
		float *frac, int n, float *lon, float *lat) {
  if (fabs(lat1) > 87.5 || fabs(lat2) > 87.5) {
    float x1, y1;
    float x2, y2;
    float x, y;
    int hemi=fabs(lat1)/lat1;
    tcoord2_2lonlat(x1, y1, -1, hemi, lon1, lat1);
    tcoord2_2lonlat(x2, y2, -1, hemi, lon2, lat2);
    for (int i=0; i<n; i++) {
      x=x1+(x2-x1)*frac[i];
      y=y1+(y2-y1)*frac[i];
      tcoord2_2lonlat(x, y, 1, hemi, lon[i], lat[i]);
    }
  } else {
    float dlon;
    if (lon1-lon2 > 180) {
      dlon=(lon2-lon1-360);
    } else if (lon1-lon2 < -180) {
      dlon=(lon2-lon1+360);
    } else {
      dlon=(lon2-lon1);
    }
    for (int i=0; i<n; i++) {
      lon[i]=lon1+dlon*frac[i];
      if (lon[i] < 0) lon[i]=360+lon[i];
      lon[i]=lon[i]-360*int (lon[i]/360);
      lat[i]=lat1+(lat2-lat1)*frac[i];
    }
  }
}

int32_t read_bev_file(char *filename, int32_t ind, time_class &t, float *&x, float *&y) {
  FILE *fs;
  int32_t magic;
  int32_t n;
  int32_t npt;
  short yr, mon, dy, hr, min;
  float sec;

  fs=fopen(filename, "r");

  fread(&magic, sizeof(magic), 1, fs);
  assert(magic == MAGIC);
  fread(&n, sizeof(n), 1, fs);

  for (int32_t i=0; i<ind; i++) {
    fread(&yr, sizeof(yr), 1, fs);
    fread(&mon, sizeof(mon), 1, fs);
    fread(&dy, sizeof(dy), 1, fs);
    fread(&hr, sizeof(hr), 1, fs);
    fread(&min, sizeof(min), 1, fs);
    fread(&sec, sizeof(sec), 1, fs);
    fread(&npt, sizeof(npt), 1, fs);
    fseek(fs, 8*npt, SEEK_CUR);
  }

  fread(&yr, sizeof(yr), 1, fs);
  fread(&mon, sizeof(mon), 1, fs);
  fread(&dy, sizeof(dy), 1, fs);
  fread(&hr, sizeof(hr), 1, fs);
  fread(&min, sizeof(min), 1, fs);
  fread(&sec, sizeof(sec), 1, fs);
  fread(&npt, sizeof(npt), 1, fs);

  x=new float[npt];
  y=new float[npt];

  fread(x, sizeof(float), npt, fs);
  fread(y, sizeof(float), npt, fs);

  fclose(fs);

  return npt;
}

int32_t bev_index(char *filename, time_class *&t, int32_t *&npt) {
  FILE *fs;
  int32_t magic;
  int32_t n;
  short yr, mon, dy, hr, min;
  float sec;

  fs=fopen(filename, "r");

  fread(&magic, sizeof(magic), 1, fs);
  assert(magic == MAGIC);
  fread(&n, sizeof(n), 1, fs);

  t=new time_class[n];
  npt=new int32_t[n];

  for (int32_t i=0; i<n; i++) {
    fread(&yr, sizeof(yr), 1, fs);
    fread(&mon, sizeof(mon), 1, fs);
    fread(&dy, sizeof(dy), 1, fs);
    fread(&hr, sizeof(hr), 1, fs);
    fread(&min, sizeof(min), 1, fs);
    fread(&sec, sizeof(sec), 1, fs);
    fread(npt+i, sizeof(int32_t), 1, fs);
    t[i].init(yr, mon, dy, hr, min, sec);
    fseek(fs, 8*npt[i], SEEK_CUR);
  }
  return n;
}

long read_ascii_contour(FILE *fs, float *&x, float *&y) {
  linked_list<float> xdata;
  linked_list<float> ydata;
  float x1, y1;
  char *check;
  int ngot;
  long npt;
  char line[MAXLL];

  while (feof(fs) == 0) {
    check=fgets(line, MAXLL, fs);
    if (check==NULL) break;
    ngot=sscanf(line, "%f %f", &x1, &y1);
    if (ngot<2) break;
    xdata.add(x1);
    ydata.add(y1);
  }

  x=xdata.make_array(npt);
  y=ydata.make_array(npt);

  return npt;
}

//determines whether or not two points cross
//a boundary
int crosses_boundary(float lon1, float lat1, float lon2, float lat2, float *lon, float *lat, long n, int wrap_flag) {

  float minlon_v;
  float maxlon_v;
  float minlon;
  float maxlon;

  float lon_i, lon_ip1, lat_i, lat_ip1;
  float m1, b1, m2, b2;
  float xint;

  if (lon1 < lon2) {
    minlon_v=lon1;
    maxlon_v=lon2;
  } else {
    minlon_v=lon2;
    maxlon_v=lon1;
  }

  for (long i=0; i<n; i++) {
    lon_i=lon[i];
    lon_ip1=lon[i+1];
    lat_i=lat[i];
    lat_ip1=lat[i+1];

    if (fabs(lon_i-lon_ip1) > 180) {
      if (lon_i < lon_ip1) {
        lon_i=lon_i+360;
      } else {
        lon_ip1=lon_ip1+360;
      }
    }
    
    m1=(lat_i-lat_ip1)/(lon_i-lon_ip1);
    b1=lat_i-m1*lon_i;
    m2=(lat1-lat2)/(lon1-lon2);
    b2=lat1-m2*lon1;
    
    xint=(b2-b1)/(m1-m2);

    if (lon_i < lon_ip1) {
      minlon=lon_i;
      maxlon=lon_ip1;
    } else {
      minlon=lon_ip1;
      maxlon=lon_i;
    }

    if (xint > minlon && xint > minlon_v && xint < maxlon && xint < maxlon_v) {
      return 1;
    }
  }

  if (wrap_flag==0) return 0;
  
  lon_i=lon[n-1];
  lon_ip1=lon[0];
  lat_i=lat[n-1];
  lat_ip1=lat[0];

  if (fabs(lon_i-lon_ip1) > 180) {
    if (lon_i < lon_ip1) {
      lon_i=lon_i+360;
    } else {
      lon_ip1=lon_ip1+360;
    }
  }
    
  m1=(lat_i-lat_ip1)/(lon_i-lon_ip1);
  b1=lat_i-m1*lon_i;
  m2=(lat1-lat2)/(lon1-lon2);
  b2=lat1-m2*lon2;
    
  xint=(b2-b1)/(m1-m2);

  if (lon_i < lon_ip1) {
    minlon=lon_i;
    maxlon=lon_ip1;
  } else {
    minlon=lon_ip1;
    maxlon=lon_i;
  }

  if (xint > minlon && xint > minlon_v && xint < maxlon && xint < maxlon_v) {
    return 1;
  }
  return 0;

}


//returns the uncertainty fraction of a boundary, given an error value:
float uncertainty_fraction(float *lon, float *lat, long npt, float epsilon, int hemi,  
		int min_uncertain, int max_total, int wrap_flag) {

  time_class t1, t2;
  meas_data *data;
  float angle, late, lone;
  
  float unf;
  int uncertain, total;

  uncertain=0;
  total=0;
  while (uncertain <= min_uncertain && total <= max_total) {
//      print, total, uncertain
    data=generate_random_global(t1, t2, 1, hemi);
    angle=2*M_PI*ranu();
    late=epsilon*sin(angle)/KMPERDEG+data[0].lat;
    lone=epsilon*cos(angle)/KMPERDEG/cos(M_PI*(late+data[0].lat)/360)+data[0].lon;
    if (lone < 0) {
      lone=lone+360;
      data[0].lon=data[0].lon+360;
    }
    if (crosses_boundary(data[0].lon, data[0].lat, lone, late, lon, lat, npt, wrap_flag)) uncertain=uncertain+1;
    total=total+1;
    delete [] data;
  }
  unf=1.0*uncertain/total;
  
  return unf;
  
}

//counts the number of boxes in an az. eq. set of bins:
long count_boxes_global(float *lon,
			float *lat, 
			long npt,
			float eps,		//box size
			int hemi,
			int wrapflag)
{

  //int *bin;
  char *bin2;			//for validation
  float x, y;			//coords in az-eq system
  int xind, yind;		//bin indices
  long nbox, nbox1, nbox2;	//for validation
  int hemi2;			//hemisphere used for coord-trans
  int hemiold;
  long n=ceil(20000/eps);
  		//number of bins per side (n x n)
  float ds;			//distance between adjacent points
  float xold, yold;		//previous point
  int noold=1;			//previous point not valid (must be calculated)
  float xtween, ytween;		//interpolated points
  int ntween;			//number of points to interpolate
  int i0, i1;
  int baseind;
  bit_array *bin;		//bins for counting boxes

  if (hemi==0) bin=new bit_array(2*n*n, 0); else bin=new bit_array(n*n, 0);

  if (wrapflag) {
    i0=npt-1;
    i1=0;
  } else {
    i0=0;
    i1=1;
  }

  if (hemi==0) {
    hemi2=fabs(lat[i0])/lat[i0];
    baseind=(hemi2+1)/2;
  } else {
    hemi2=abs(hemi)/hemi;
    baseind=0;
    if (hemi<0 && lat[i0] > 0) {
      goto skip;
    } else if (hemi>0 && lat[i0] < 0) {
      goto skip;
    }
  }
      
    tcoord2_2lonlat(xold, yold, -1, hemi2, lon[i0], lat[i0]);
    xind=n*(xold+10000)/(20000);
    yind=n*(yold+10000)/(20000);
    bin->on(baseind*n*n+xind*n+yind);
    noold=0;
  skip:
  hemiold=hemi2;

  for (long i=i1; i<npt; i++) {
    if (hemi==0) {
      hemi2=fabs(lat[i])/lat[i];
      baseind=(hemi2+1)/2;
    } else if (hemi<0 && lat[i] > 0) {
      noold=1;
      continue;
    } else if (hemi>0 && lat[i] < 0) {
      noold=1;
      continue;
    }
      
    tcoord2_2lonlat(x, y, -1, hemi2, lon[i], lat[i]);
    xind=n*(x+10000)/(20000);
    yind=n*(y+10000)/(20000);
    //bin[(hemi2+1)*n*n/2+xind*n+yind]++;
    //bin2[(hemi2+1)*n*n/2+xind*n+yind]=1;
    bin->on(baseind*n*n+xind*n+yind);

    //boxes won't cover every single piece of the curve, but the boxes
    //should nonetheless be contiguous, which as all we need...
    if (noold) {
      tcoord2_2lonlat(xold, yold, -1, hemi2, lon[i-1], lat[i-1]);
    } else if (hemi2 != hemiold) {
      //if the line segment is right on the equator, flip everything
      //around and do it twice...
      tcoord_N2S(x, y);
      baseind=1-baseind;
      for (int j=1; j<ntween; j++) {
        xtween=xold+(x-xold)*j/ntween;
        ytween=yold+(y-yold)*j/ntween;
        //damn, this is really crude...
        if (xtween*xtween+ytween*ytween > 1e8) continue;
        xind=n*(xtween+10000)/(20000);
        yind=n*(ytween+10000)/(20000);
        bin->on(baseind*n*n+xind*n+yind);
      }
      tcoord_N2S(x, y);
      tcoord_N2S(xold, yold);
      baseind=1-baseind;
    }
    ds=sqrt(tcoord_ds2(x, y, xold, yold));
    ntween=n/eps+1;
    for (int j=1; j<ntween; j++) {
      xtween=xold+(x-xold)*j/ntween;
      ytween=yold+(y-yold)*j/ntween;
      //damn, this is really crude...
      if (xtween*xtween+ytween*ytween > 1e8) continue;
      xind=n*(xtween+10000)/(20000);
      yind=n*(ytween+10000)/(20000);
      bin->on(baseind*n*n+xind*n+yind);
    }

    xold=x;
    yold=y;
    hemiold=hemi2;
    noold=0;
  }

  if (hemi==0) baseind=2; else baseind=1;
  nbox=bin->nnonzero(baseind*n*n-1);
  delete bin;
  return nbox;

  nbox=0;
  nbox1=0;
  for (long i=0; i<2*n*n; i++) {
    printf("%d %d\n", (*bin)[i], bin2[i]);
    //printf("%d", bin2[i]);
    nbox+=bin2[i];
    nbox1+=(*bin)[i];
  }
  printf("\n");
  //bin.print();
  delete [] bin2;
  nbox2=bin->nnonzero(2*n*n-1);
  printf("nbox=%d %d %d\n", nbox, nbox1, nbox2);

  return nbox;
}

long measure_contour(float *lon,
                        float *lat,
                        long npt,
                        float ruler,
                        int wrap_flag)
{
  float d1, d2, ds;
  float frac;
  long nseg, nsegnew;
  float lonc, latc;
  float lon1, lon2;

  long j, k;

  nseg=0;
  d1=0;
  lonc=lon[0];
  latc=lat[0];

  for (long i=1; i<npt; i++) {
    d2=sqrt(sdist(lonc, latc, lon[i], lat[i]));
    if (d2 >= ruler) {
      k=i-1;
      ds=sqrt(sdist(lon[k], lat[k], lon[i], lat[i]));
      frac=(ruler-d1)/(d2-d1);
      //printf("frac=%g\n", frac);
      nsegnew=ds*(1-frac)/ruler;
      frac=frac+ruler*nsegnew/ds;
      interpol_tween_coords(lon[k], lat[k], lon[i], lat[i], &frac, 1, &lonc, &latc);
      nseg=nseg+nsegnew+1;
      d1=ds*(1-frac);
    } else {
      d1=d2;
    }
  }

  if (wrap_flag) {
    for (long i=0; i<npt; i++) {
      d2=sqrt(sdist(lonc, latc, lon[i], lat[i]));
      if (d2 >= ruler) {
        k=i-1;
        if (k<0) k=npt+k;
        ds=sqrt(sdist(lon[k], lat[k], lon[i], lat[i]));
        frac=(ruler-d1)/(d2-d1);
        nsegnew=ds*(1-frac)/ruler;
        nseg+=nsegnew+1;
        if (i==0) nseg++;
        break;
      }
      d1=d2;
    }
  }
      
  return nseg;
}

long contour_interpolate(float *lon, float *lat, long n, float ds, float *&lonnew, float *&latnew, int wrap_flag) {
  int ntween[n];
  float ds0;
  long nnew;

  for (long i=1; i<n; i++) {
    ds0=sdist(lon[i-1], lat[i-1], lon[i], lat[i]);
    ntween[i]=ds0/ds;
    if (ds0-ds*ntween[i] != 0) ntween[i]++;
    nnew+=ntween[i];
  }

  if (wrap_flag) {
    ds0=sdist(lon[n-1], lat[n-1], lon[0], lat[0]);
    ntween[0]=ds0/ds;
    if (ds0-ds*ntween[0] != 0) ntween[0]++;
    nnew+=ntween[0];
  }
   
  lonnew=new float[nnew];
  latnew=new float[nnew];

  for (long i=1; i<n; i++) {
  }

}

/*
long measure_contour(float *lon,
			float *lat,
			long npt,
			float ruler,
			int wrap_flag)
{
  float d1, d2, ds;
  float frac;
  long nseg, nsegnew;
  float lonc, latc;
  float lon1, lon2;

  long i,j, k;
  int wrap=0;

  nseg=0;
  d2=0;
  i=0;
  lonc=lon[0];
  latc=lat[0];
  do {
    while (d2 < ruler) {
      j++;
      d1=d2;
      if (j>=npt) {
        if (wrap_flag==0) return nseg;
        j=0;
        wrap=1;
      }
      d2=sqrt(sdist(lonc, latc, lon[j], lat[j]));
    }
    if (j==0) k=npt-1; else k=j-1;
    ds=sqrt(sdist(lon[k], lat[k], lon[j], lat[j]));
    frac=(ruler-d1)/(d2-d1);
    nsegnew=ds*(1-frac)/ruler;
    frac+=ruler*nsegnew/ds;
    d1=0;
    d2=ds*(1-frac);
    lon1=lon[k]; lon2=lon[j];
    if (fabs(lon1-lon2) > 360) {
      if (lon2 > 90) lon2=lon2-360; else lon2=lon2+360;
    }
    lonc=lon1+(lon2-lon1)*(frac);
    latc=lat[k]+(lat[j]-lat[k])*(frac);
    //printf("frac=%f\n", frac);
    nseg+=nsegnew+1;
    //printf("%d: (%f %f); %d; %d\n", i, lonc, latc, nsegnew, nseg);
  } while (wrap==0);

  return nseg;
}
*/

} //end namespace ctraj

