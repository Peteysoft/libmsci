#include <stdio.h>
#include <stdint.h>

#include <gsl/gsl_fft_complex.h>

#include "peteys_tmpl_lib.h"
#include "error_codes.h"
#include "time_class.h"
#include "parse_command_opts.h"

#include "ctraj_defaults.h"
#include "pp_util.h"
#include "ctraj_3d_fields.h"

using namespace ctraj;
using namespace libpetey;

#define MAXNREC 1000

int main(int argc, char **argv) {
  char *fname;
  float zlev;		//get 2-D wind fields for this vertical level
  float **data[MAXNREC];
  int32_t *header[MAXNREC];
  int nrec;
  float *zgrid;
  int nlev;
  int nlon, nlat;
  //fields:
  float ***u0;		//zonal wind
  float ***v0;		//meridional wind
  float ***w0;		//vertical wind
  float ***t;		//temperature
  float ***z;		//geopotential height
  float ***u;		//on same grid as t and z
  float ***v;
  float ***uout;
  float ***vout;
  int flags[20];
  void *optargs[20];

  //parse command options:
  argc=parse_command_opts(argc, argv, "TQ?", "%%%",
                  optargs, flags, OPT_WHITESPACE);
  if (argc < 0) exit(21);

  fname=argv[1];

  //print out help screen:  
  if (flags[2] || (argc < 3 && flags[1]!=1) || (flags[1] && argc < 2)) {
    FILE *docfs;
    int err;
    if (flags[2]) {
      docfs=stdout;
      err=0;
    } else {
      docfs=stderr;
      err=INSUFFICIENT_COMMAND_ARGS;
    }
    fprintf(docfs, "\n");
    fprintf(docfs, "Syntax:  read_pp [-Q] [-T] file level\n");
    fprintf(docfs, "\n");
    fprintf(docfs, "where:\n");
    fprintf(docfs, "  file     =  file containing UKMO data in 'pp' format\n");
    fprintf(docfs, "  level    =  level to extract\n");
    fprintf(docfs, "\n");
    fprintf(docfs, "options:\n");
    ctraj_optargs(docfs, "QT?");
    fprintf(docfs, "\n");
    return err;
  }

  //read in all the data:
  nrec=pp_read_all(fname, header, data, MAXNREC);

  //just print out a list of fields contained in the file:
  if (flags[1]) {
    printf(" no.       date    code     level  nlat nlon\n");
    for (int i=0; i<nrec; i++) {
      //printf("Record %d:\n", i);
      //field code, vertical level, nlat, nlon, lat0, dlat, lon0, dlon:
        printf("%4d %d/%d/%d-%d %4d %10.3f %4d %4d\n", i+1, 
			header[i][0], header[i][1], header[i][2], header[i][3],
			header[i][PP_HEADLOC_CODE], 
		    ((float *)(header[i]))[PP_HEADLOC_LEV], 
		    header[i][PP_HEADLOC_NLAT], header[i][PP_HEADLOC_NLON]);
        //printf("%g %g %g %g\n",
	//	    ((float *)(header[i]))[PP_HEADLOC_LAT0], 
	//	    ((float *)(header[i]))[PP_HEADLOC_DLAT], 
	//	    ((float *)(header[i]))[PP_HEADLOC_LON0], 
	//	    ((float *)(header[i]))[PP_HEADLOC_DLON]);
    }
    //print out first 38 header records:
    //for (int j=0; j<38; j++) printf("%d %d\n", j, header[i][j]);

    //print out 50-64:
    //for (int j=49; j<HEADLEN; j++) printf("%d %g\n", j, ((float *) header[i])[j]);
    //printf("\n");
    //printf("\n");
    goto finish;
  }
  
  zlev=atof(argv[2]);
  zgrid=pp_extract_levels(header, nrec, nlev);
  //for (int i=0; i<nlev; i++) printf("%g\n", zgrid[i]);

  //we expect the data in specific format pulled out by this function:
  pp_extract_uvwtz(data, header, nrec, u0, v0, w0, t, z);

  //interpolate wind fields to same grid as temperature (one more lat. grid):
  nlon=header[0][PP_HEADLOC_NLON];
  nlat=header[0][PP_HEADLOC_NLAT];
  u=allocate_3D_field<float>(nlev, nlat+1, nlon);
  v=allocate_3D_field<float>(nlev, nlat+1, nlon);

  pp_interpolate_uv(u0, v0, nlev, nlat, nlon, u, v);

  if (flags[0]) {
    //interpolated to potential-temperature level:
    float ***theta;		//pot. temp.
    double ***coef;		//interpolation ceofficients
    //coef=pp_interpolate_pt_levels(zgrid, t, nlev, nlat+1, nlon, &zlev, 1);
    theta=calc_pot_temp(t, zgrid, nlev, nlat+1, nlon);
    coef=calc_zlev_coef(theta, nlev, nlat+1, nlon, &zlev, 1);
    uout=interpolate_zlev(u, coef, 1, nlat+1, nlon);
    vout=interpolate_zlev(v, coef, 1, nlat+1, nlon);
    delete_3D_field(coef, 1);
    delete_3D_field(theta, nlev);
  } else {
    //interpolate to pressure level:
    double coef;
    int ind;
    double frac;
    //in-line interpolation:
    reverse(zgrid, nlev);
    //need single interpolation coef.:
    coef=interpolate(zgrid, nlev, zlev);
    if (coef>nlev-1 || coef<0) {
      fprintf(stderr, "Pressure level (%g) out-of-range\n", zlev);
      exit(PARAMETER_OUT_OF_RANGE);
    }
    ind=(int) coef;
    frac=coef-ind;
    ind=nlev-ind-1;
    uout=allocate_3D_field<float>(1, nlat+1, nlon);
    vout=allocate_3D_field<float>(1, nlat+1, nlon);
    //cycle through grids:
    for (int i=0; i<nlat+1; i++) {
      for (int j=0; j<nlon; j++) {
        uout[0][i][j]=u[ind][i][j]*(1-frac);
        vout[0][i][j]=u[ind][i][j]*(1-frac);
	//returning bottom grid is possible (maybe...)
	if (ind<nlev-1) {
          uout[0][i][j] += u[ind-1][i][j]*frac;
          vout[0][i][j] += u[ind-1][i][j]*frac;
	}
      }
    }
  }

  //print_matrix(stdout, u[0], nlat+1, nlon);
  //print_matrix(stdout, uout[0], nlat+1, nlon);

  //exit(0);

  //print out 2-D wind fields:
  for (int i=nlat; i>=0; i--) {
    for (int j=0; j<nlon; j++) {
      printf("%g ", uout[0][i][j]);
    }
    printf("\n");
  }
  printf("\n");
  for (int i=nlat; i>=0; i--) {
    for (int j=0; j<nlon; j++) {
      printf("%g ", vout[0][i][j]);
    }
    printf("\n");
  }
  printf("\n");

  //interpolated to temp. grid:
  delete_3D_field(u, nlev);
  delete_3D_field(v, nlev);
  //interpolated to single vertical grid:
  delete_3D_field(uout, 1);
  delete_3D_field(vout, 1);

  finish:
    //delete all data read from file:
    for (int i=0; i<nrec; i++) {
      delete [] header[i];
      delete [] data[i][0];
      delete [] data[i];
    }

  return 0;

}
