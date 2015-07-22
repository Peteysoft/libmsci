#ifndef CTRAJ_BOUNDARY_LIST
#define CTRAJ_BOUNDARY_LIST

#include <stdio.h>

#include "metric_base.h"
#include "ctraj_vfield_base.h"

namespace ctraj {

  template <class real, class vreal>
  class ctraj_boundary_list {
    private:
      ctraj_vfield_base<vreal> *vfield;
      metric_base<vreal> *metric;

      int dflag;			//output time values as dates instead of strings

    protected:
      real *x;
      real *y;
      int32_t *domain;
      long n;
      double tind;

      real thresh_arc;
      real min_spac;
      real max_spac;

      int wrap_flag;

    public:
      ctraj_boundary_list();
      ctraj_boundary_list(ctraj_vfield_base<vreal> *v, metric_base<vreal> *m, real tarc, real mins, real maxs, int df=0);
      virtual ~ctraj_boundary_list();

      virtual long init_circle(real x0, real y0, real r);

      virtual void wrap_on();
      void wrap_off();

      void sett(double t);
      double gett();

      virtual long fix();
      virtual double advance(double tstep, int32_t nrk);

      virtual size_t read(FILE *fs);
      virtual size_t write(FILE *fs);
      virtual long print(FILE *fs);

  };

} //end namespace ctraj

#endif

