#include <assert.h>
#include <math.h>

#include <stdint.h>

#include "rk_dumb_ts.h"

#include "gsl/gsl_interp.h"
#include "gsl/gsl_spline.h"

#include "time_class.h"

#include "ctraj_boundary_list.h"

#include "ctraj_defaults.h"

using namespace libpetey;

namespace ctraj {

//default constructor:
template <class real, class vreal>
ctraj_boundary_list<real, vreal>::ctraj_boundary_list() {

  x=NULL;
  y=NULL;
  domain=NULL;
  n=0;

  tind=0;
  min_spac=0;
  max_spac=0;
  thresh_arc=0;
  wrap_flag=0;
}

template <class real, class vreal>
void ctraj_boundary_list<real, vreal>::sett(double t) {
  tind=t;
}

template <class real, class vreal>
double ctraj_boundary_list<real, vreal>::gett() {
  return tind;
}

template <class real, class vreal>
ctraj_boundary_list<real, vreal>::ctraj_boundary_list(ctraj_vfield_base<vreal> *v, metric_base<vreal> *m, real tarc, real mins, real maxs, int df) {
  vfield=v;
  metric=m;

  x=NULL;
  y=NULL;
  domain=NULL;
  n=0;

  thresh_arc=tarc;
  min_spac=mins;
  max_spac=maxs;

  wrap_flag=1;

  dflag=df;
/*
  for (ind_type i=0; i<tgrid->nel(); i++) {
    tgrid->get(tel, i);
    //tel.write_string(tstr);
    //printf("%s\n", tstr);
  }
  */

}

template <class real, class vreal>
long ctraj_boundary_list<real, vreal>::init_circle(real x0, real y0, real r) {
  real ang;

  n=(long) (2*M_PI/thresh_arc);
  x=new real[n+wrap_flag];
  y=new real[n+wrap_flag];
  domain=new int[n+wrap_flag];

  for (long i=0; i<n+wrap_flag; i++) {
    ang=i*thresh_arc;
    //no attempt to apply metric corrections:
    x[i]=x0+r*cos(ang);
    y[i]=y0+r*sin(ang);
    domain[i]=1;
  }

  return n;
}

template <class real, class vreal>
void ctraj_boundary_list<real, vreal>::wrap_on() {
  real *xnew;
  real *ynew;
  int32_t *hnew;

  if (!wrap_flag) {
    wrap_flag=1;
    if (n != 0) {
      xnew=new real[n+1];
      ynew=new real[n+1];
      hnew=new int[n+1];
      for (long i=0; i<n; i++) {
        xnew[i]=x[i];
        ynew[i]=y[i];
	hnew[i]=domain[i];
      }
      xnew[n]=x[0];
      ynew[n]=y[0];
      hnew[n]=domain[0];
      delete [] x;
      delete [] y;
      delete [] domain;
      x=xnew;
      y=ynew;
      domain=hnew;
    }
  }
}

template <class real, class vreal>
void ctraj_boundary_list<real, vreal>::wrap_off() {
  wrap_flag=0;
}

template <class real, class vreal>
ctraj_boundary_list<real, vreal>::~ctraj_boundary_list() {
  if (x != NULL) delete [] x;
  if (y != NULL) delete [] y;
  if (domain != NULL) delete [] domain;

}

template <class real, class vreal>
double ctraj_boundary_list<real, vreal>::advance(double tstep, int32_t nrk) {
  int32_t hold;
  vreal x0[2];
  time_class t1;
  char tstr[30];
  void *param[2];
  vreal **result;

  //allocate intermediate values:
  result=new vreal*[n];
  result[0]=new vreal[n*2];
  for (long i=1; i<n; i++) result[i]=result[0]+i*2;

  //pass the grids to the global variables:
  param[0]=vfield;

  for (long i=0; i<n; i++) {
    x0[0]=x[i];
    x0[1]=y[i];
    domain[i]=vfield->fix(domain[i], tind, x0);
    param[1]=domain+i;
    rk_dumb(tind, (vreal *) x0, 2L, tstep, nrk, result, (void *) param, &ctraj_deriv<vreal>);
    x[i]=result[nrk][0];
    y[i]=result[nrk][1];
  }

  tind+=tstep*nrk;

  delete [] result[0];
  delete [] result;

  return tind;
}

template <class real, class vreal>
long ctraj_boundary_list<real, vreal>::fix() {
  long nn=n+wrap_flag;
  real s[nn];		//the parametric distance along the curve
  real dx2;		//2nd derivative of x wrt s
  real dy2;		//2nd derivative of y wrt s
  real dels;		//the distance between a pair of points
  real dstotal;
  real dangle;		//the radians of arc (^2) between a pair of points
  real curv2;		//the curvature squared
  int nnew[nn-1];	//number of new points to insert between each pair of points
  long totaln, dn;	//final total number of points
  real * xnew;
  real * ynew;
  int * hnew;
  long j, k;
  long offset;		//for removing redundant nodes
  int32_t refd;		//reference domain

  vreal mcoef[2];
  real gcurv;

  //we do the calculations with this:
  vreal x1[2], x2[2];

  gsl_interp *xinterp;
  gsl_interp *yinterp;
  gsl_interp_accel *xaccel;
  gsl_interp_accel *yaccel;

  gsl_interp_type *spltype;

  if (wrap_flag) {
    x[n]=x[0];
    y[n]=y[0];
    domain[n]=domain[0];
  }

  s[0]=0;
  //we work in the N. hemisphere:
  if (domain[0] == 0) {
    x1[0]=x[0];
    x1[1]=y[0];
    refd=vfield->reference(domain[0], x1);
    x[0]=x1[0];
    y[0]=x1[1];
  }
  offset=0;
  //printf("nn=%d\n", nn);
  for (long i=1; i<nn; i++) {
    if (offset != 0) {
      //because of any redundant nodes, all the nodes must be shifted:
      j=i-offset;
      x[j]=x[i];
      y[j]=y[i];
      domain[j]=domain[i];
    } else {
      j=i;
    }
    x2[0]=x[j];
    x2[1]=y[j];
    vfield->reference(domain[0], x2);
    x[j]=x2[0];
    y[j]=x2[1];
    x1[0]=x[j-1];
    x1[1]=y[j-1];

    dels=metric->ds2(x1, x2);
    //printf("%f\n", dels);
    assert(dels>=0);
    s[j]=s[j-1]+sqrt(dels);
    //printf("%g %g\n", 
    //		sqrt((x[j-1]-x[j])*(x[j-1]-x[j])+(y[j-1]-y[j])*(y[j-1]-y[j])), 
    //		sqrt(dels));
    //check for nodes made redundant by round-off:
    if (s[j]-s[j-1] <= 0) {
      printf("Warning: removing %ldth node made redundant by round-off\n", i);
      offset++;
    }
    //printf("%d (%d %g, %g): %g\n", i, domain[i], x[i], y[i], dels);
  }
  n-=offset;
  nn-=offset;

  //do a cubic spline interpolation on both variables:
  if (wrap_flag) {
    //closed contour implies periodic boundary conditions:
    spltype=(gsl_interp_type *) gsl_interp_cspline_periodic;
//    printf("dx0ds=%g; dy0ds=%g\n", dx0, dy0);
  } else {
    spltype=(gsl_interp_type *) gsl_interp_cspline;
  }

  //initialize our cubic spline fitters:
  xinterp=gsl_interp_alloc(spltype, nn);
  yinterp=gsl_interp_alloc(spltype, nn);

  xaccel=gsl_interp_accel_alloc();
  yaccel=gsl_interp_accel_alloc();

  //start the cubic spline interpolation:
  gsl_interp_init(xinterp, s, x, nn);
  gsl_interp_init(yinterp, s, y, nn);

  //calculate the angle traced out by each pair of points:
  //and use to figure out how many new points to add...
  totaln=1;
  dangle=0;
  dstotal=0;
//  printf("i, s, dx2ds2, dy2ds2, dangle, #new\n");
  for (long i=1; i<nn; i++) {
    dels=s[i]-s[i-1];
    dstotal+=dels;
    assert(dels>=0);
    //if the point spacing is less than the threshold, either
    //remove next point or add no new:
    if (dstotal < min_spac) {
      if (2*dstotal < min_spac) {
        nnew[i-1]=-1;		//remove point
      } else {
        nnew[i-1]=0;
      }
    } else {
      //get second order derivatives:
      gsl_interp_eval_deriv2_e(xinterp, s, x, (s[i-1]+s[i])/2, xaccel, &dx2);
      gsl_interp_eval_deriv2_e(xinterp, s, y, (s[i-1]+s[i])/2, yaccel, &dy2);
      x1[0]=(x[i-1]+x[i])/2;
      x1[1]=(y[i-1]+y[i])/2;

      //get properties of space:
      metric->mcoef2(x1, mcoef);
      //printf("mcoef=%g, %g\n", mcoef[0], mcoef[1]);
      gcurv=metric->gcurv(x1);

      //use that to calculate the sweep of arc:
      dangle+=sqrt(dx2*dx2*mcoef[0]+dy2*dy2*mcoef[1]+gcurv)*dels;
      //dangle+=sqrt(dx2*dx2+dy2*dy2+gcurv)*dels;

      //printf("r=%g\n", dels/dangle);

      //if the accumulated angle of arc is greater than the threshold,
      //add new points, otherwise, remove the next point...
      if (dangle*2 < thresh_arc) {
        nnew[i-1]=-1;
      } else {
        nnew[i-1]=(int) (dangle/thresh_arc-0.5);
        //dangle=dangle-(nnew[i-1]+0.5)*thresh_arc;
	dangle=0;
      }

      //to correct for instabilities caused by 
      //points being too widely separated:
      if (dstotal > max_spac) {
        long checkn=(long) (dstotal/max_spac);
        if (checkn > nnew[i-1]) nnew[i-1]=checkn;
      }
      dstotal=0;
    }
    totaln+=nnew[i-1]+1;
    //printf("%d, %g, %g, %g, %g, %d\n", i, s[i], dx2, dy2, dangle, nnew[i-1]);
  }

  xnew=new real[totaln];
  ynew=new real[totaln];

  j=1;			//index of new nodes
  xnew[0]=x[0];		//always keeps the first point
  ynew[0]=y[0];

  k=0;			//start interpolating from this node
  for (long i=1; i<nn; i++) {
    if (nnew[i-1] > 0) {
      dstotal=s[i]-s[k];
      dels=dstotal/(nnew[i-1]+1);
      if (2*dels < min_spac) {
        dn=nnew[i-1]-(int) (dstotal/min_spac*2);
        nnew[i-1]-=dn;
	totaln-=dn;
        dels=dstotal/(nnew[i-1]+1);
      }
      //add in new points interpolated along the spline:
      for (long m=1; m <= nnew[i-1]; m++) {
	gsl_interp_eval_e(xinterp, s, x, s[k]+m*dels, xaccel, xnew+j);
	gsl_interp_eval_e(yinterp, s, y, s[k]+m*dels, yaccel, ynew+j);
        //printf("%g %g\n", xnew[j], ynew[j]);
	j++;
      }
      xnew[j]=x[i];
      ynew[j]=y[i];
      k=i;
      j++;
    } else if (nnew[i-1] == 0) {
      //no new points to add:
      xnew[j]=x[i];
      ynew[j]=y[i];
      k=i;
      j++;
    }
  }

  assert(j==totaln);

  delete [] x;
  delete [] y;
  x=xnew;
  y=ynew;

  delete [] domain;
  domain=new int[totaln];
  for (long i=0; i<totaln; i++) domain[i]=refd;

  n=totaln-wrap_flag;

  //clean up:
  gsl_interp_free(xinterp);
  gsl_interp_free(yinterp);
  gsl_interp_accel_free(xaccel);
  gsl_interp_accel_free(yaccel);

  return n;

}

//to match the different sizes of types encountered on different
//architectures:
typedef float read_real;

template <class real, class vreal>
size_t ctraj_boundary_list<real, vreal>::read(FILE *fs) {
  int16_t yy, mon, dd, hh, min;
  float sec;
  //time_class start;
  size_t nread;
  int32_t readn;
  time_class date;
  char tstring[TFIELD_WIDTH];

  vreal x0[2], x1[2];

  read_real *lon, *lat;

  if (x != NULL) delete [] x;
  if (y != NULL) delete [] y;
  if (domain != NULL) delete [] domain;

  if (dflag) {
    nread=fread(&yy, sizeof(yy), 1, fs);
    nread+=fread(&mon, sizeof(mon), 1, fs);
    nread+=fread(&dd, sizeof(dd), 1, fs);
    nread+=fread(&hh, sizeof(hh), 1, fs);
    nread+=fread(&min, sizeof(min), 1, fs);
    nread+=fread(&sec, sizeof(sec), 1, fs);
    date.init(yy, mon, dd, hh, min, sec);
    //bit silly, convert it to a string, then convert it back again...
    date.write_string(tstring);
  } else {
    nread=fread(tstring, sizeof(char), TFIELD_WIDTH+1, fs);
  }

  tind=vfield->get_tind(tstring);

  nread+=fread(&readn, sizeof(readn), 1, fs);
  n=readn;

  lon=new read_real[n+wrap_flag];
  lat=new read_real[n+wrap_flag];
  domain=new int[n+wrap_flag];

  nread+=fread(lon, sizeof(read_real), n, fs);
  nread+=fread(lat, sizeof(read_real), n, fs);
/*
  for (long i=0; i<n; i++) {
    printf("%f %f\n", lon[i], lat[i]);
  }
*/

  x=new real[n+wrap_flag];
  y=new real[n+wrap_flag];

  printf("%s\n", tstring);

  //convert from lon-lat:
  for (long i=0; i<n; i++) {
    domain[i]=0;
    x0[0]=lon[i];
    x0[1]=lat[i];
    domain[i]=vfield->absolute(-1, x0);
    x[i]=x0[0];
    y[i]=x0[1];
    //printf("(%g, %g); (%d %g, %g)\n", lon[i], lat[i], domain[i], x[i], y[i]);
    //tcoord2_2lonlat(x[i], y[i], -1, domain[i], lond, latd);
  }

/*
  for (long i=0; i<n; i++) {
    printf("%f %f %d\n", x[i], y[i], hemi[i]);
  }
*/

  delete [] lon;
  delete [] lat;

  return nread;

}

template <class real, class vreal>
size_t ctraj_boundary_list<real, vreal>::write(FILE *fs) {
  read_real lon[n], lat[n];
  int16_t yy, mon, dd, hh, min;
  float sec;
  size_t nwrit;
  int32_t readn;
  time_class date;
  char tstring[TFIELD_WIDTH];
  vreal x0[2], x1[2];

  vfield->get_t(tind, tstring);

  if (dflag) {
    date.read_string(tstring);
    date.get_fields(yy, mon, dd, hh, min, sec);

    nwrit=fwrite(&yy, sizeof(yy), 1, fs);
    nwrit+=fwrite(&mon, sizeof(mon), 1, fs);
    nwrit+=fwrite(&dd, sizeof(dd), 1, fs);
    nwrit+=fwrite(&hh, sizeof(hh), 1, fs);
    nwrit+=fwrite(&min, sizeof(min), 1, fs);
    nwrit+=fwrite(&sec, sizeof(sec), 1, fs);
  } else {
    nwrit=fwrite(tstring, sizeof(char), TFIELD_WIDTH+1, fs);
  }

  readn=n;
  nwrit+=fwrite(&readn, sizeof(readn), 1, fs);

  //convert to lon-lat:
  for (long i=0; i<n; i++) {
    x0[0]=x[i];
    x0[1]=y[i];
    vfield->absolute(domain[i], x0);
    lon[i]=x0[0];
    lat[i]=x0[1];
    //tcoord2_2lonlat(x[i], y[i], 1, hemi[i], lon1, lat1);
    //printf("%g %g %d %g %g\n", x[i], y[i], hemi[i], lon[i], lat[i]);
  }

  nwrit+=fwrite(lon, sizeof(read_real), n, fs);
  nwrit+=fwrite(lat, sizeof(read_real), n, fs);

  return nwrit;

}

template <class real, class vreal>
long ctraj_boundary_list<real, vreal>::print(FILE *fs) {
  for (long i=0; i<n; i++) {
    printf("%10.6g %10.6g %d\n", x[i], y[i], domain[i]);
  }
  return n;
}

template class ctraj_boundary_list<double, float>;

} //end namespace ctraj

