#ifndef CTRAJ_TRACER_STANDARD__H
#define CTRAJ_TRACER_STANDARD__H 1

#include "simple_temp.h"
#include "dependent_swap.h"
#include "composite_dataset.h"

#include "dependent_intc.h"

#include "ctraj_tfield_base.h"

#include "az_eq_t.h"

namespace ctraj {
  using namespace libpetey;
  using namespace datasets;

  template <class real>
  class ctraj_tfield_standard:public ctraj_tfield_base<real> {
    protected:
      simple<real> *xgrid;               //tracer x grid
      simple<real> *ygrid;               //tracer y grid
      az_eq_t<real> *metric;

      //dummy tracer field for calculating interpolation coefficients:
      dependent_intc *tracer;
      //mapping from tracer field to vector to be held in sparse matrix:
      sub_1d_type nmap;
      dependent<sub_1d_type> *map_map;

      sub_1d_type *inverse_map;

    public:
      ctraj_tfield_standard();
      virtual ~ctraj_tfield_standard();

      //initialize object:
      int init1(int32_t np, real sld2=SIDELENGTH_Q);	//points/side, sidelength/2
      int init2(int32_t n, real sld2=SIDELENGTH_Q);	//total # points, sidelength/2

      //set the metric:
      void set_metric(az_eq_t<real> *m);

      //sets up the object based on command line arguments:
      virtual int setup(int argc, char **argv);

      virtual int32_t nel();		//return number of elements
      virtual int32_t ndim();		//return number of dimensions
      virtual int32_t nwt();		//return number of interpolation weights

      //gets the location of the <ind>th point in the field:
      virtual int32_t 			//return domain
		get_loc(int32_t ind, 		//index of point
		real *loc);			//relative coordinates

      //returns interpolation coefficients for a given location:
      virtual int32_t 			//number of weighting coefficients
		interpolate(int32_t domain, 	//domain of point
		real *loc, 			//relative coordinates
		int32_t *ind, 			//point indices of weighting coefficients
		double *wt,			//weighting coefficients
		real dt=1);

      virtual real *to(real *input, real *parm);
      virtual real *from(real *input, real *parm);

      //get raw grids (if appropriate):
      virtual int get_raw_grids(int32_t *grids);
      //get range:
      virtual int get_range(real *low, real *high);

      virtual void help(FILE *fs);
  };

} //end namespace ctraj

#endif

