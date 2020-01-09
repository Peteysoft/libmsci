#ifndef PP_UTIL__H
#define PP_UTIL__H

#include <stdio.h>
#include <stdint.h>

#include "error_codes.h"
#include "time_class.h"

#include "ctraj_defaults.h"

using namespace libpetey;

namespace ctraj {

  //kind of a dumb way to do it; should probably put everything in a structure:
  const int PP_HEADLEN=64;

  //locations of various things in the header:
  const int PP_HEADLOC_CODE=22;
  const int PP_HEADLOC_NLAT=17;
  const int PP_HEADLOC_NLON=18;
  const int PP_HEADLOC_N=14;
  const int PP_HEADLOC_LEV=51;
  const int PP_HEADLOC_LAT0=58;
  const int PP_HEADLOC_LON0=60;
  const int PP_HEADLOC_DLAT=59;
  const int PP_HEADLOC_DLON=61;

  //byte-swaps an array of 4-byte integers:
  void swap_endian (int32_t *data, int n);

  //reads all the headers and data from a file in "pp" format:
  int pp_read_all(char *fname, 			//name of file
		  int32_t **headers_all,	//fields headers
		  float ***fields, 		//fields
		  int nmax);			//maximum number to read in
					//(size of pre-allocated return values)

  //reads all the headers and data from a file in "pp" format:
  int pp_read_field(FILE *fs, 			//name of file
		  int32_t **headers_all,	//header for each field
		  float ***fields, 		//fields
		  int nmax,			//maximum number to read in
		  int field_code,		//code of field to read in
		  float *plev=NULL,		//pressure levels to read in
		  int nlev=0,			//number of levels (0=all)
		  int toendflag=0);		//read to end of file

  //codes for different fields:
  const int PP_U_CODE=56;		//zonal wind
  const int PP_V_CODE=57;		//meridional wind
  const int PP_Z_CODE=1;		//geopotential height
  const int PP_T_CODE=16;		//temperature
  const int PP_W_CODE=40;		//vertical wind

  //extracts fields from pp file
  //assumes following fields are contained in this order: u, v, z, t, w
  //*** note no new fields are allocated, but piggy-back off original data
  //throws error if format is wrong
  int pp_extract_uvwtz(float ***data, 		//data in as-read format
		  int32_t **header, 		//headers
		  int n,			//number of fields
		float ***&u, 			//3-D zonal wind field
		float ***&v, 			//3-D meridional wind field
		float ***&w, 			//  " vertical wind field
		float ***&t, 			//  " temerature
		float ***&z);			//  " geopotential height

  //extracts vertical levels from headers
  //throws error if vertical levels are not all the same
  float * pp_extract_levels(int32_t **header, 	//headers
		  int n, 			//number of fields
		  int &nlev);			//returned number of levels

  //interpolates zonal and meridional winds to the same horizontal grids as
  //temperature and geopotential height fields
  int pp_interpolate_uv(float ***u, 		//zonal wind
		  float ***v, 			//meridional wind
		  int nlev, 			//number of vertical levels
		  int nlat, 			//# lat. grids for u & v
		  int nlon, 			//number of longitude grids
		  float ***unew, 		//returned wind fields
		  float ***vnew);		//**must be pre-allocated**

  //finds interpolation coefficients for a set of potential temperature levels
  //returns 3-D field containing coefficients
  float *** pp_interpolate_pt_levels(float *plev, 	//pressure levels
		  float ***t,				//temp. field
		int nlev, 				//number of z grids
		int nlat, 				//number of lat. grids
		int nlon, 				//number of lon. grids
		float *ptlev, 				//pot. temp. levels
		int npt, 				//# pot. temp. grids
		float pref=1000);			//reference pressure

  //interpolates to potential temperature levels
  //takes as input returned int. coefficients from above
  //returns interpolated field
  float *** pp_zinterpolate(float ***field, 	//field in pressure levels
		  float ***coef, 		//interpolation coef.
		  int nlev, 			//number of pot. temp. grids
		  int nlat, 			//number of lat. grids
		  int nlon);			//number of lon. grids

} //end namespace ctraj

#endif
