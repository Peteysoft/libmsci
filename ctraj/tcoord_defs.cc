#include <math.h>
#include "ctraj_defaults.h"

namespace ctraj {

  template<class real>
  real tcoord_ds2(real x1, real y1, real x2, real y2) {
    real xave, yave;	//average of x and y
    real r2;
    real xave2, yave2;
    real sangle;		//sine of the zenith angle

    real xc, yc;		//metric coefficients
    real xdiff, ydiff;

    xave=(x1+x2)/2;
    yave=(y1+y2)/2;
    xave2=xave*xave;
    yave2=yave*yave;
    r2=xave2+yave2;
    sangle=sin(sqrt(r2)/REARTH);
    //printf("r=%g; sin(angle)=%gi; <x>= %g, <y>= %g\n", sqrt(r2), sangle, xave, yave);
    xc=(REARTH*REARTH*yave2*sangle*sangle/r2+xave2)/r2;
    yc=(REARTH*REARTH*xave2*sangle*sangle/r2+yave2)/r2;

    //printf("r=%g; xc=%g; yc=%g; ds=%g\n", sqrt(r2), xc, yc, sqrt(xc*xdiff*xdiff+yc*ydiff*ydiff));

    xdiff=x2-x1;
    ydiff=y2-y1;

    return xc*xdiff*xdiff+yc*ydiff*ydiff;

  }

  template<class real>
  void tcoord2_2lonlat(real &x, real &y, short dir, short hemi, real &lon, real &lat) {
    real r;

    if (dir > 0) {
      r=sqrt(x*x+y*y);
      lat=hemi*(90-r/KMPERDEG);
      lon=atan(y/x)*RAD2DEG;

      if (x < 0) lon=180+lon;
      if (lon < 0) lon=360+lon;
    } else if (dir < 0) {

      if (hemi == 0) hemi=(short) (lat/fabs(lat));

      //do the coordinate transformation:
      r=(90-hemi*lat)*KMPERDEG;
      x=r*cos(lon*DEG2RAD);
      y=r*sin(lon*DEG2RAD);
    }
  }

  template <class real1, class real2>
  void tcoord_fixH(real1 &x, real1 &y, short &hemi, real2 h[4]) {
    real2 x2, y2;
    real2 r2, r, rnew;
    real2 dxdx, dxdy, dydx, dydy;

    x2=x*x;
    y2=y*y;
    r2=x2+y2;
    r=sqrt(r);
    if (r>10000.) {
      hemi=-hemi;
      rnew=20000.-r;

      x=x*rnew/r;
      y=y*rnew/r;

      dxdx=20000./r*(1-x2/r2)-1.;
      dydy=20000./r*(1-y2/r2)-1.;
      dxdy=20000.*y2/r2/r;
      dydx=20000.*x2/r2/r;

      h[0]=h[0]*dxdx+h[2]*dxdy;
      h[1]=h[1]*dxdx+h[3]*dxdy;
      h[2]=h[0]*dydx+h[2]*dydy;
      h[3]=h[1]*dydx+h[3]*dydy;
    }
  }

  template float tcoord_ds2(float x1, float y1, float x2, float y2);
  template double tcoord_ds2(double x1, double y1, double x2, double y2);

  template void tcoord2_2lonlat(float &x, float &y, short dir, short hemi, float &lon, float &lat);
  template void tcoord2_2lonlat(double &x, double &y, short dir, short hemi, double &lon, double &lat);

  template void tcoord_fixH(float &x, float &y, short &hemi, float h[4]);
  template void tcoord_fixH(float &x, float &y, short &hemi, double h[4]);
  template void tcoord_fixH(double &x, double &y, short &hemi, float h[4]);
  template void tcoord_fixH(double &x, double &y, short &hemi, double h[4]);

} //end namespace ctraj

