#ifndef ctraj_TRACER_ANAL_H
#define ctraj_TRACER_ANAL_H

#include <stdint.h>

#include "sparse.h"
#include "meas_data.h"

#include "ctraj_vfield_base.h"
#include "ctraj_tfield_base.h"

namespace ctraj {
  using namespace libpetey;
  using namespace libsparse;

  //this function is completely redundant:
  float ** tracer_multiply(sparse_matrix *a, int32_t N, float *v0);

  //maps a tracer, re-normalizing at each step:
  float ** tracer_multiply_renorm(sparse_matrix *a, 	//tracer mapping
		int32_t N, 			//number of time grids
		float *v0);			//initial tracer field

  //interpolate a tracer
  //assumes standard side-length of 20000 km
  void tracer_interp(float **qall, 		//tracer fields
		time_class *t, 			//time grids
		int32_t nt, 			//number of time grids
		int32_t np, 			//number of horizontal grids
		meas_data *data, 		//list of locations/interpolates
		long nsamp);			//number of interpolates

  //performs the pc-proxy interpolation:
  float ** pc_proxy(sparse_matrix *matall, 	//tracer mapping
		time_class *t, 			//time grids
		int32_t N, 			//"lead" time
		int32_t nall, 			//total number of grids
		meas_data *samp, 		//list of samples
		long nsamp, 			//number of samples
		int32_t nev, 			//number of eigenvalues
		int32_t ncv, 			//number of Arnoldi vectors
		int cflag=0,			//constant term?
		int mflag=0);			//output multiple fields?

  //returns the coefficients:
  //(doesn't work--one of these days I should fix it)
  float * proxy_tracer(float **tall, 
		int32_t np, 
		time_class *t, 
		int32_t nall, 
		meas_data *samp, 
		long nsamp, 
		int32_t order);

  //new paradigm:
  template <class real>
  void tracer_map(ctraj_tfield_base<real> * tracer, 
		ctraj_vfield_base<real> *vfield,
		sparse<int32_t, real> *map,
		real t, real dt, int32_t nt);

  //calculate approximate area of each grid point:
  template <typename real>
  void grid_area(int32_t n, 		//number of points
		  real *area);		//area at each point

  //calculate equivalent latitude:
  template <typename real>
  class eq_lat {
    real *area;			//area of each grid point
    int32_t n;			//number of points
    real total_area;		//total area

    public:
      eq_lat(int32_t np);
      real operator () (real *q, real *el);
      ~eq_lat();
  };

} //end namespace ctraj

#endif
