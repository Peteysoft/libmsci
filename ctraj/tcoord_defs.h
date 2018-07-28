#ifndef ctraj_TCOORD_DEFS__H
#define ctraj_TCOORD_DEFS__H

#include <math.h>

//this module should be deprecated soon...

namespace ctraj {

  //returns the squared metric, but without using the f2vector type:
  template <class real>
  real tcoord_ds2(real x1, real y1, real x2, real y2);

  template <class real>
  void tcoord2_2lonlat(real &x, real &y, short dir, short hemi, real &lon, real &lat);

  template <class real>
  inline void tcoord_fix(real &x, real &y, short &hemi) {
    float r, rnew;

    r=sqrt(x*x+y*y);
    if (r>10000) {
      hemi=-hemi;
      rnew=20000-r;

      x=x*rnew/r;
      y=y*rnew/r;
    }
  }

  template <class real>
  inline void tcoord_N2S(real &x, real &y) {
    real r, rnew;

    r=sqrt(x*x+y*y);
    rnew=20000-r;

    x=x*rnew/r;
    y=y*rnew/r;
  }

  template <class real1, class real2>
  void tcoord_fixH(real1 &x, real1 &y, short &hemi, real2 h[4]);

} //end namespace ctraj

#endif

