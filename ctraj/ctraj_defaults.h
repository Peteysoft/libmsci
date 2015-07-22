#ifndef CTRAJ_DEFAULTS_H
#define CTRAJ_DEFAULTS_H 1

#include <stdio.h>
#include <stdint.h>

#include "ctraj_vfield_base.h"

//default values, coefficients and a bunch of other stuff:

namespace ctraj {

  //default gridding:
  const int32_t NGRID=120;
  const float SIDELENGTH=12000.;

  //for the tracer field we don't need any overlap:
  const int32_t NGRID_Q=100;		//a smaller grid size is probably a good idea as well...
  const float SIDELENGTH_Q=10000.;

  //default conversion to spherical-polar coords:
  const int32_t NLON=360;
  const int32_t NLAT=181;

  //ECMWF gridding:
  const int32_t ECMWF_NLON=240;
  const int32_t ECMWF_NLAT=121;

  //defaults for contour advection:
  const float MAXARC=1.;
  const float DSMAX=10.;
  const float DSMIN=1.;

  //deprecated:
  //const int64_t C3_BLOCK_SIZE=67108864;		//64 megs

  //defaults for pc-proxy:
  const int32_t NARNOLDI=20;
  const int32_t NEIG=5;

  //number of plot contours:
  const int32_t NZ=21;

  //time steps:
  const float TSTEP_COARSE=1.;
  const int32_t TSTEP_NFINE=6;

  //write interval:
  const int WRITE_INT=1;

  //date field width:
  const int TFIELD_WIDTH=23;

  //maximum size of date string conversions:
  const int TSTRING_LEN=30;

  //do we allow whitespace between the option and it's parameter?
  const int OPT_WHITESPACE=1;

  //for computing uncertainty exponent:
  const float EPS_MIN=1.;
  const float EPS_MAX=1000.;
  const int32_t NEPS=20;
  const int32_t NUNC_MIN=100;
  const int32_t UNC_MAXN=10000;

  //"magic" number for contour advection files:
  const int32_t MAGIC=-55;

  //coefficients:
  //unit conversions:
  const float KMPERDEG=111.11;
  const float RAD2DEG=57.295779511;
  const float REARTH=6366.2;
  const float DEG2RAD=0.017453292521;
  const float HOURSPERDAY=24;
  const float MPERKM=1000;
  const float SECPERDAY=86400;

  const float KAPPA=0.2857;
  const float P0=1000;			//reference pressure for potential temperature
  const float ABS0=273.2;
  const float TOFFS=-2;

  //where to switch from curvilinear metric to straight Cartesian:
  const float LAT_THRESH=87.5;

  template <class real>
  int ctraj_deriv(double t, real *x, real *v, void *param);

  template <class real>
  int ctraj_deriv_L(real t, real *x, real *v, void *param);

  void ctraj_optargs(FILE *fs, const char *optargs, int flag=0);

  ctraj_vfield_base<float> *ctraj_loader(int &argc, char **argv, double &I0, int32_t &n);

}  //end namespace ctraj

#endif

