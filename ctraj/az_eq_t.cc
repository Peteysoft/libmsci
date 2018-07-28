#include <math.h>
#include "ctraj_defaults.h"

#include "az_eq_t.h"

namespace ctraj {

template<class real>
az_eq_t<real>::az_eq_t(real r) {
  rearth=r;
  cearth=2*r*M_PI;
  kmperdeg=cearth/360;
}

template<class real>
real az_eq_t<real>::ds2(real *x1, real *x2) {
  real xave, yave;	//average of x and y
  real r2;
  real xave2, yave2;
  real sangle;		//sine of the zenith angle

  real xc, yc;		//metric coefficients
  real xdiff, ydiff;

  xave=(x1[0]+x2[0])/2;
  yave=(x1[1]+x2[1])/2;
  xave2=xave*xave;
  yave2=yave*yave;
  r2=xave2+yave2;
  sangle=sin(sqrt(r2)/rearth);
  //printf("r=%g; sin(angle)=%gi; <x>= %g, <y>= %g\n", sqrt(r2), sangle, xave, yave);
  xc=(rearth*rearth*yave2*sangle*sangle/r2+xave2)/r2;
  yc=(rearth*rearth*xave2*sangle*sangle/r2+yave2)/r2;

  //printf("r=%g; xc=%g; yc=%g; ds=%g\n", sqrt(r2), xc, yc, sqrt(xc*xdiff*xdiff+yc*ydiff*ydiff));

  xdiff=x2[0]-x1[0];
  ydiff=x2[1]-x1[1];

  return xc*xdiff*xdiff+yc*ydiff*ydiff;

}

template<class real>
void az_eq_t<real>::mcoef2(real *x, real *c) {
  real r2, sangle;
  real x2, y2;

  x2=x[0]*x[0];
  y2=x[1]*x[1];
  r2=x2+y2;
  sangle=sin(sqrt(r2)/rearth);
  c[0]=(rearth*rearth*y2*sangle*sangle/r2+x2)/r2;
  c[1]=(rearth*rearth*x2*sangle*sangle/r2+y2)/r2;

}

template<class real>
real az_eq_t<real>::gcurv(real *x) {
  return 1/rearth/rearth;
}

template<class real>
int az_eq_t<real>::ndim() {
  return 2;
}

template<class real>
real az_eq_t<real>::getr() {
  return rearth;
}

template<class real>
void az_eq_t<real>::from(int hemi, real *x, real *lonlat) {
  real r;

  r=sqrt(x[0]*x[0]+x[1]*x[1]);
  lonlat[1]=hemi*(90-r/kmperdeg);
  lonlat[0]=atan(x[1]/x[0])*RAD2DEG;

  if (x[0] < 0) lonlat[0]=180+lonlat[0];
  if (lonlat[0] < 0) lonlat[0]=360+lonlat[0];

}

template<class real>
void az_eq_t<real>::to(real *lonlat, int hemi, real *x) {
  real r;

  //do the coordinate transformation:
  r=(90-hemi*lonlat[1])*cearth/360;
  x[0]=r*cos(lonlat[0]*DEG2RAD);
  x[1]=r*sin(lonlat[0]*DEG2RAD);

}

template<class real>
int az_eq_t<real>::to(real *lonlat, real *x) {
  int hemi;
  //if it's right on the equator, stick it in the N. hemisphere:
  if (lonlat[1]<0) hemi=-1; else hemi=1;
  to(lonlat, hemi, x);
  return hemi;
}

template <class real>
int az_eq_t<real>::fix(int hemi, real *x) {
  float r, rnew;

  r=sqrt(x[0]*x[0]+x[1]*x[1]);
  if (r>cearth/4) {
    hemi=-hemi;
    rnew=cearth/2-r;

    x[0]=x[0]*rnew/r;
    x[1]=x[1]*rnew/r;
  }

  return hemi;
}

template <class real>
void az_eq_t<real>::swap(real *x) {
  real r, rnew;

  r=sqrt(x[0]*x[0]+x[1]*x[1]);
  rnew=cearth/2-r;

  x[0]=x[0]*rnew/r;
  x[1]=x[1]*rnew/r;
}


template <class real>
int az_eq_t<real>::fixH(int hemi, real *x, real h[4]) {
  real x2, y2;
  real r2, r, rnew;
  real dxdx, dxdy, dydx, dydy;

  x2=x[0]*x[0];
  y2=x[1]*x[1];
  r2=x2+y2;
  r=sqrt(r);
  if (r>cearth/4) {
    hemi=-hemi;
    rnew=cearth/2-r;

    x[0]=x[0]*rnew/r;
    x[1]=x[1]*rnew/r;

    dxdx=cearth/2/r*(1-x2/r2)-1.;
    dydy=cearth/2/r*(1-y2/r2)-1.;
    dxdy=cearth*y2/r2/r/2;
    dydx=cearth*x2/r2/r/2;

    h[0]=h[0]*dxdx+h[2]*dxdy;
    h[1]=h[1]*dxdx+h[3]*dxdy;
    h[2]=h[0]*dydx+h[2]*dydy;
    h[3]=h[1]*dydx+h[3]*dydy;
  }
  return hemi;
}

template <class real>
void az_eq_t<real>::vtran(int hemi, real *x, real *v, real *dxdt) {
  real r, val;
 
  r=sqrt(x[0]*x[0]+x[1]*x[1]);
  val=rearth*sin(r/rearth);

  dxdt[0]=-hemi*v[1]*x[0]/r-v[0]*x[1]/val;
  dxdt[1]=-hemi*v[1]*x[1]/r+v[0]*x[0]/val;

}

template class az_eq_t<float>;
template class az_eq_t<double>;

} //end namespace ctraj

