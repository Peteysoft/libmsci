#ifndef ctraj_AZ_EQ_T__H
#define ctraj_AZ_EQ_T__H

#include <math.h>

#include "metric_base.h"

namespace ctraj {

  //azimuthal equidistant coord. transforms contained in an object class
  //(variable radius parameter...)
  template <class real>
  class az_eq_t:public metric_base<real> {
    protected:
      real rearth;
      real cearth;
      real kmperdeg;
    public:
      az_eq_t(real r);

      virtual real ds2(real * x1, real *x2);
      virtual void mcoef2(real *x, real *c);
      virtual real gcurv(real *x);
      virtual int ndim();

      //transform from lon-lat to az. eq. coords
      //**note: hemi is a parameter, not a return value...
      void to(real *lonlat, int hemi, real *x);
      //here hemi is a return value:
      int to(real *lonlat, real *x);

      //transform from az. eq. coords to lon-lat:
      void from(int hemi, real *x, real *lonlat);

      int fix(int hemi, real *x);

      void swap(real *x);

      int fixH(int hemi, real *x, real h[4]);

      void vtran(int hemi, real *x, real *v, real *dxdt);

      real getr();
  };

} //end namespace ctraj

#endif
