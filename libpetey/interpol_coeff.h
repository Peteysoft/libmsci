#ifndef _INTERPOL_COEFF_H_
#define _INTERPOL_COEFF_H_

//performs multi-linear interpolation for any dimension
int interpol_coeff(int32_t rank, //rank or dimension
		int32_t *dim,			//number of grids each dimension 
		double *indices,		//interpolation indices for each grid
		int64_t *subscripts,		//returned 1-d subscripts
		double *coeffs);		//interpolation coeffs


#endif

