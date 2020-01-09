#include <stdio.h>
#include <stdint.h>

#include <gsl/gsl_fft_complex.h>

#include "peteys_tmpl_lib.h"
#include "error_codes.h"
#include "full_util.h"
#include "time_class.h"

#include "ctraj_3d_fields.h"

#include "ctraj_defaults.h"

using namespace libpetey;

namespace ctraj {

  template <typename scalar>
  scalar ***allocate_3D_field(int nlev, int ny, int nx) {
    scalar ***field;
    field=new scalar **[nlev];
    for (int i=0; i<nlev; i++) field[i]=allocate_matrix<scalar, int32_t>(ny, nx);
    return field;
  }

  template <typename scalar>
  void delete_3D_field(scalar ***field, int nlev) {
    for (int i=0; i<nlev; i++) delete_matrix(field[i]);
    delete [] field;
  }

  //finds interpolation coefficients for a set of potential temperature levels
  template <typename real>
  real *** calc_pot_temp(real ***t, real *plev,
		int nlev, int nlat, int nlon, 
		real pref) {
    real ***pt;

    pt=allocate_3D_field<real>(nlev, nlat, nlon);

    for (int k=0; k<nlev; k++) {
      for (int i=0; i<nlat; i++) {
        for (int j=0; j<nlon; j++) {
          pt[k][i][j]=t[k][i][j]*pow(pref/plev[k], KAPPA);
        }
      }
    }

    return pt;
  }

  //finds interpolation coefficients for a set of vertical levels
  template <typename real>
  double *** calc_zlev_coef(real ***z,
		int nlev, int nlat, int nlon, 
		real *zlev, int nz) {
    double ***coef;			//returned coefficients
    real col[nlev];			//single vertical column
    long sind[nlev];			//sort indices
    long loc;				//interpolation index
    long lastind=-1;

    coef=allocate_3D_field<double>(nz, nlat, nlon);

    for (int i=0; i<nlat; i++) {
      for (int j=0; j<nlon; j++) {
        for (int k=0; k<nlev; k++) {
          col[k]=z[k][i][j];
	  //printf("%g ", pt[k]);
        }
        //printf("\n");
        heapsort(col, sind, nlev);
        map_vector_inplace(col, sind, nlev);
        for (int k=0; k<nz; k++) {
          loc=bin_search(col, nlev, zlev[k], lastind);
          if (sind[loc]>=nlev-1 || loc < 0) {
            fprintf(stderr, "Vertical level (%g) falls out of sky\n", zlev[k]);
            throw PARAMETER_OUT_OF_RANGE;
          }
          coef[k][i][j]=sind[(int) loc]+(zlev[k]-col[loc])/(col[loc+1]-col[loc]);
          //printf("%g ", coef[k][i][j]);
        }
      }
      //printf("\n");
    }

    return coef;
  }

  template <typename real>
  real *** interpolate_zlev(real ***field, double ***coef, int nlev, int nlat, int nlon) {
    real *** fnew;
    fnew=new real**[nlev];
    for (int i=0; i<nlev; i++) {
      fnew[i]=allocate_matrix<real, int32_t>(nlat, nlon);
      for (int j=0; j<nlat; j++) {
        for (int k=0; k<nlon; k++) {
          int pind=(int) coef[i][j][k];
          double frac=coef[i][j][k]-pind;
          fnew[i][j][k]=field[pind][j][k]*(1-frac)
		+field[pind+1][j][k]*frac;
          //printf("%g ", field[pind][j][k]);
        }
        //printf("\n");
      }
    }
    return fnew;
  }

  template float ***allocate_3D_field<float>(int nlev, int ny, int nx);
  template double ***allocate_3D_field<double>(int nlev, int ny, int nx);

  template void delete_3D_field<float>(float ***field, int nlev);
  template void delete_3D_field<double>(double ***field, int nlev);

  template float *** calc_pot_temp<float>(float ***, float *, int, int, int, float); 
  template double *** calc_pot_temp<double>(double ***, double *, int, int, int, double); 

  template double *** calc_zlev_coef<float>(float ***, int, int, int, float *, int);
  template double *** calc_zlev_coef<double>(double ***, int, int, int, double *, int);

  template float *** interpolate_zlev<float>(float ***, double ***, int, int, int);
  template double *** interpolate_zlev<double>(double ***, double ***, int, int, int);
} //end namespace ctraj
