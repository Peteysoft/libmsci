//lon-lat ascii to azimuthal-equidistant binary...

#include "simple_temp.h"
#include "dependent_temp.h"

//converts tracer field as
//dataset in lon-lat coords 
//to tracer field as vector in az. eq. coords to tracer field as
float * ll2ae(dependent<float> *q1, 	//input field in lon-lat coords (must wrap)
		simple<float> *lon, 	//longitude grid
		simple<float> *lat, 	//latitude grid
		long ngrid);		//number of grid points per side in az. eq. coords

//converts tracer field as vector in az. eq. coords to tracer field as
//dataset in lon-lat coords
long ae2ll(float *qvec, 		//input field as vector in az. eq. coords
		long ngrid,		//number of grids per side in az. eq. coords
		simple<float> *lon, 	//longitude grid
		simple<float> *lat, 	//latitude grid
		dependent<float> *q1);  //output field in lon-lat coords

