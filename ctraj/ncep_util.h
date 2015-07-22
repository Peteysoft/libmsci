#ifndef ctraj_NCEP_UTIL__H
#define ctraj_NCEP_UTIL__H

#include <stdlib.h>

#include <netcdfcpp.h>

#include "time_class.h"
#include "simple_temp.h"
#include "dependent_temp.h"

#include "ctraj_defaults.h"

namespace ctraj {

  using namespace libpetey;
  using namespace datasets;

  int get_ncep_dim(NcFile *nci,			//netcdf file handle 
		char *name, 
		simple<float> *&dim); 		//longitude grid

  //read in the grids from a 4-D ncep file:
  int get_ncep_grid(NcFile *nci,			//netcdf file handle 
		simple<float> *&lon, 		//longitude grid
		simple<float> *&lat, 		//latitude grid
		simple<float> *&lev,		//vertical grid
		simple<time_class> *&tgrid);	//time grid

  //read in a 3-D field for one time index:
  int get_ncep_t(NcFile *nci,			//netcdf file handle
		const char *var,
		long tind,			//time index
		dependent<float> *q);		//field

  int get_ncep_surf_t(NcFile *nci,		//netcdf file handle
		const char *var,
		long tind,			//time index
		dependent<float> *q);		//field

  //read in a 3-D field for one time index:
  int get_ncep_theta_level(NcFile *nci,			//netcdf file handle
		const char *var,
		long tind,			//time index
		dependent<double> *c,		//interpolation coeffs
		dependent<float> *q);		//field

  //read in a 2-D field for time and vertical (z) direction:
  int get_ncep_tz(NcFile *nci,			//netcdf file handle
		const char *var,
		long tind,			//time index
		double zindex,			//z index
		dependent<float> *q);		//field

  //read in a 2-field for time and theta level:
  int get_ncep_theta_interp(NcFile *nc_T,		//temperatures
		float theta_level,		//theta level
		long tind,			//time index
		dependent<double> *c, 		//output field
		float kappa=KAPPA);

  //read in a 2-field for time and theta level:
  int get_dthdp(NcFile *nc_T,		//temperatures
		float theta_level,		//theta level
		long tind,			//time index
		dependent<double> *c, 		//output field
		dependent<float> *dthdp, 	//output field
		float kappa=KAPPA);

  //take interpolation coefficients and interpolate within a field:
  int theta_interp(dependent<float> *q1,
		dependent<double> *c,
		dependent<float> *q2);

  //get time grids for multiple files:
  simple<time_class> * get_ncep_tgrid(const char *base, time_class t1, time_class t2,
		long &it1, long &it2);

  //number of time grids in a file:
  inline long ncep_nt(NcFile *nc) {
    NcVar *ncv;

    ncv=nc->get_var("time");
    return ncv->num_vals();
  }

} //end namespace ctraj

#endif
