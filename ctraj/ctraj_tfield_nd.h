#ifndef CTRAJ_TRACER_2D__H
#define CTRAJ_TRACER_2D__H 1

#include <stdint.h>

#include "error_codes.h"

#include "simple_temp.h"
#include "dependent_temp.h"
#include "dependent_intc.h"

#include "ctraj_defaults.h"
#include "ctraj_tfield_base.h"

namespace ctraj {
  using namespace libpetey::datasets;

  template <class real>
  class ctraj_tfield_nd:public ctraj_tfield_base<real> {
    protected:
      int32_t n_dim;

      simple<real> **grid;
      dependent_intc *tracer;

      sub_1d_type nmap;
      dependent<sub_1d_type> *map;
      sub_1d_type *inverse_map;

      void get_loc_raw(sub_1d_type ind, real *loc);
    public:
      ctraj_tfield_nd(int32_t nd=0);
      virtual ~ctraj_tfield_nd();

      //sets up the object based on command line arguments:
      virtual int setup(int argc, char **argv);

      virtual int32_t nel();		//return number of elements
      virtual int32_t ndim();
      virtual int32_t nwt();

      //gets the location of the <ind>th point in the field:
      virtual int32_t 			//return v-field domain
		get_loc(int32_t ind, 		//index of point
		real *loc);			//relative coordinates

      //returns interpolation coefficients for a given location:
      virtual int32_t 			//number of weighting coefficients
		interpolate(int32_t domain, 	//domain of point
		real *loc, 			//relative coordinates
		int32_t *ind, 			//point indices of weighting coefficients
		double *wt);			//weighting coefficients

      virtual real *to(real *input, real *parm);
      virtual real *from(real *input, real *parm);

      //get raw grids (if appropriate):
      virtual int get_raw_grids(int32_t *grids);
      //get range:
      virtual int get_range(real *low, real *high);

      void help(FILE *fs);

  };

} //end namespace ctraj

#endif

