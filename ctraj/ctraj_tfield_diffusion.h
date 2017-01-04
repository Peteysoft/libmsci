#ifndef CTRAJ_TFIELD_DIFFUSION__H
#define CTRAJ_TFIELD_DIFFUSION__H 1

#include "simple_temp.h"
#include "dependent_swap.h"
#include "composite_dataset.h"

#include "dependent_intc.h"

#include "ctraj_tfield_standard.h"

namespace ctraj {
  using namespace libpetey;
  using namespace datasets;

  template <class real>
  class ctraj_tfield_diffusion:public ctraj_tfield_standard<real> {
    protected:
      real dcoeff;		//diffusion coefficient
      real cutoff;		//cut-off size for weight

    public:
      ctraj_tfield_diffusion();
      virtual ~ctraj_tfield_diffusion();

      //initialize object:
      int init(int32_t np, 			//points pers side
		      real sld2=SIDELENGTH_Q, 	//sidlength/2
		      real d=DIFFUSION, 	//diffusion coefficient
		      real co=GAUSS_CUTOFF);	//cut-off value for Gaussian kernel

      //sets up the object based on command line arguments:
      virtual int setup(int argc, char **argv);

      //how many weight do you need, maximum?
      virtual int32_t nwt();

      //returns interpolation coefficients for a given location:
      virtual int32_t 			//number of weighting coefficients
		interpolate(int32_t domain, 	//domain of point
		real *loc, 			//relative coordinates
		int32_t *ind, 			//point indices of weighting coefficients
		double *wt,			//weighting coefficients
		real dt=1);

      virtual void help(FILE *fs);
  };

} //end namespace ctraj

#endif

