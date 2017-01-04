#ifndef CTRAJ_TRACER_BASE__H
#define CTRAJ_TRACER_BASE__H 1

#include <stdio.h>

namespace ctraj {

template <class real>
class ctraj_tfield_base {
  public:
    virtual ~ctraj_tfield_base();
    //sets up the object based on command line arguments:
    virtual int setup(int argc, char **argv)=0;

    virtual int32_t nel()=0;		//return number of elements
    virtual int32_t ndim()=0;		//return number of dimensions
    virtual int32_t nwt()=0;		//return maximum number of weights
					//needed for interpolation

    //gets the location of the <ind>th point in the field:
    virtual int32_t 			//return v-field domain
	get_loc(int32_t ind, 		//index of point
	real *loc)=0;			//relative coordinates

    //returns interpolation coefficients for a given location:
    virtual int32_t 			//number of weighting coefficients
	interpolate(int32_t domain, 	//domain of point
	real *loc, 			//relative coordinates
	int32_t *ind, 			//point indices of weighting coefficients
	double *wt,			//weighting coefficients
	real dt=1)=0;			//time-step (for tunable diffusion)

    //convert from standard input format (lon-lat, with corners, etc.)
    //to the vector used by the codes:
    virtual real *to(real *input, real *parm)=0;
    virtual real *from(real *input, real *parm)=0;

    //get raw grids (if appropriate):
    virtual int get_raw_grids(int32_t *grids)=0;
    //get range:
    virtual int get_range(real *low, real *high)=0;

    //print help screen:
    virtual void help(FILE *fs)=0;

};

template <class real>
ctraj_tfield_base<real>::~ctraj_tfield_base() {};

} //end namespace ctraj

#endif

