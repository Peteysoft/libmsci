#ifndef CTRAJ_3D_FIELDS__H
#define CTRAJ_3D_FIELDS__H

#include <stdio.h>
#include <stdint.h>

#include "error_codes.h"
#include "time_class.h"

#include "ctraj_defaults.h"

using namespace libpetey;

namespace ctraj {

  //allocates a 3-D field for holding weather data
  //formated as an array of individually-allocated 2-D matrices
  template <typename scalar>
  scalar ***allocate_3D_field(int nlev, 	//number of levels (matrices)
		  int ny, 			//number of latitude grids (m)
		  int nx);			//number of longitude grids (n)

  //deletes the 3-D field allocated above:
  template <typename scalar>
  void delete_3D_field(scalar ***field, int nlev);

  //calculate potential temperature:
  template <typename real>
  real *** calc_pot_temp(real ***t,			//temperature field
		  real *plev,				//pressure levels
		  int nlev,
		  int nlat,
		  int nlon,
		  real pref=P0);			//reference pressure

  //finds interpolation coefficients for a set of vertical levels
  //returns 3-D field containing coefficients
  template <typename real>
  double *** calc_zlev_coef(real ***z, 			//geopot. height or similar var.
		int nlev, 				//number of vertical grids
		int nlat, 				//number of lat. grids
		int nlon, 				//number of lon. grids
		real *zlev, 				//vertical levels
		int nz); 				//# pot. temp. grids

  //interpolates to vertical levels
  //takes as input returned int. coefficients from above
  //returns interpolated field
  template <typename real>
  real *** interpolate_zlev(real ***field, 	//field in pressure levels
		  double ***coef, 		//interpolation coef.
		  int nlev, 			//number of pot. temp. grids
		  int nlat, 			//number of lat. grids
		  int nlon);			//number of lon. grids

} //end namespace ctraj

#endif
