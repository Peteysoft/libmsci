#include <stdio.h>
#include <stdint.h>

#include <gsl/gsl_fft_complex.h>

#include "peteys_tmpl_lib.h"
#include "error_codes.h"
#include "full_util.h"
#include "time_class.h"
#include "parse_command_opts.h"

#include "ctraj_defaults.h"

using namespace ctraj;
using namespace libpetey;

//kind of a dumb way to do it; should probably put everything in a structure:
#define HEADLEN 64
//locations of various things in the header:
#define PP_HEADLOC_CODE 22
#define PP_HEADLOC_NLAT 17
#define PP_HEADLOC_NLON 18
#define PP_HEADLOC_N 14
#define PP_HEADLOC_LEV 51
#define PP_HEADLOC_LAT0 58
#define PP_HEADLOC_LON0 60
#define PP_HEADLOC_DLAT 59
#define PP_HEADLOC_DLON 61

void swap_endian (int32_t *data, int n) {
  char *fourbyte;
  char swp;
  for (int i=0; i<n; i++) {
    fourbyte=(char *) (data+i);
    swp=fourbyte[0];
    fourbyte[0]=fourbyte[3];
    fourbyte[3]=swp;
    swp=fourbyte[1];
    fourbyte[1]=fourbyte[2];
    fourbyte[2]=swp;
  }
}

template <typename scalar>
scalar ***allocate_3D_field(int nlev, int ny, int nx) {
  scalar ***field;
  field=new float **[nlev];
  for (int i=0; i<nlev; i++) field[i]=allocate_matrix<scalar, int32_t>(ny, nx);
  return field;
}

template <typename scalar>
void delete_3D_field(scalar ***field, int nlev) {
  for (int i=0; i<nlev; i++) delete [] field[i];
  delete [] field;
}


int pp_read_all(char *fname, int32_t **headers_all, float ***fields, int nmax) {

  FILE *fs;
  int32_t f77recsep;
  int32_t *header;
  float **data;
  int readcount;
  int nrec;

  fs=fopen(fname, "r");
  if (fs==NULL) {
    fprintf(stderr, "Error: unable to open %s for input\n", fname);
    throw UNABLE_TO_OPEN_FILE_FOR_READING;
  }

  for (nrec=0; nrec<nmax; nrec++) {
    readcount=fread(&f77recsep, sizeof(f77recsep), 1, fs);
    //printf("%d %d\n", nrec, readcount);
    if (readcount==0) {
      break;
    }
    swap_endian(&f77recsep, 1);

    //should be 256:
    //printf("%d\n", f77recsep);
    if (f77recsep != HEADLEN*4) {
      fprintf(stderr, "Error: error in record separator--wrong file type\n");
      fprintf(stderr, " actual: %d; expected: %d\n", f77recsep, HEADLEN*4);
      throw FILE_READ_ERROR;
    }

    header=new int32_t[HEADLEN];
    fread(header, sizeof(int32_t), HEADLEN, fs);
    swap_endian(header, HEADLEN);
    headers_all[nrec]=header;

    fread(&f77recsep, sizeof(f77recsep), 1, fs);
    fread(&f77recsep, sizeof(f77recsep), 1, fs);
    swap_endian(&f77recsep, 1);

    //should be 6912*4:
    //printf("%d\n", f77recsep);

    if (header[PP_HEADLOC_N]!=header[PP_HEADLOC_NLAT]*header[PP_HEADLOC_NLON]) {
      fprintf(stderr, "Error in record header: %d*%d != %d\n", header[PP_HEADLOC_NLAT], header[PP_HEADLOC_NLON], header[PP_HEADLOC_N]);
      throw FILE_READ_ERROR;
    }

    data=allocate_matrix<float, int32_t>(header[PP_HEADLOC_NLAT], header[PP_HEADLOC_NLON]);

    fread(data[0], sizeof(float), header[PP_HEADLOC_N], fs);
    swap_endian((int32_t *) data[0], header[PP_HEADLOC_N]);
    fields[nrec]=data;

    fread(&f77recsep, sizeof(f77recsep), 1, fs);

  }

  fclose(fs);

  return nrec;

}

#define PP_U_CODE 56
#define PP_V_CODE 57
#define PP_Z_CODE 1
#define PP_T_CODE 16
#define PP_W_CODE 40

int pp_extract_uvwtz(float ***data, int32_t **header, int n,
		float ***&u, float ***&v, float ***&w, float ***&t, float ***&z) {
  int loc=0;

  if (header[loc][PP_HEADLOC_CODE] !=PP_U_CODE) goto fail;
  
  u=data+loc;
  for (loc=0; loc<n; loc++) {
    if (header[loc][PP_HEADLOC_CODE] != PP_U_CODE) break;
  }
    
  if (header[loc][PP_HEADLOC_CODE] !=PP_V_CODE || loc>=n) goto fail;
  v=data+loc;
  for ( ; loc<n; loc++) {
    if (header[loc][PP_HEADLOC_CODE] != PP_V_CODE) break;
  }
    
  if (header[loc][PP_HEADLOC_CODE] !=PP_Z_CODE || loc>=n) goto fail;
  z=data+loc;
  for ( ; loc<n; loc++) {
    if (header[loc][PP_HEADLOC_CODE] != PP_Z_CODE) break;
  }
    
  if (header[loc][PP_HEADLOC_CODE] !=PP_T_CODE || loc>=n) goto fail;
  t=data+loc;
  for ( ; loc<n; loc++) {
    if (header[loc][PP_HEADLOC_CODE] != PP_T_CODE) break;
  }
    
  if (header[loc][PP_HEADLOC_CODE] !=PP_W_CODE || loc>=n) goto fail;
  w=data+loc;
    
  return 0;

  fail:
    fprintf(stderr, "Error: formatting of field data; %d vs. %d\n", PP_V_CODE, header[loc][PP_HEADLOC_CODE]);
    throw FILE_READ_ERROR;
}

float * pp_extract_levels(int32_t **header, int n, int &nlev) {
  int32_t code;
  float *lev;

  nlev=0;
  code=header[0][PP_HEADLOC_CODE];
  for (int i=1; i<n; i++) {
    if (code!=header[i][PP_HEADLOC_CODE]) {
      nlev=i;
      break;
    }
  }

  lev=new float[nlev];
  for (int i=0; i<nlev; i++) {
    lev[i]=((float *) header[i])[PP_HEADLOC_LEV];
  }

  for (int i=0; i<n; i++) {
    if (((float *) header[i])[PP_HEADLOC_LEV] != lev[i%nlev]) goto fail;
  }

  return lev;

  fail:
    fprintf(stderr, "Error: formatting of vertical levels\n");
    throw FILE_READ_ERROR;
}

int pp_interpolate_uv(float ***u, float ***v, int nlev, int nlat, int nlon, float ***unew, float ***vnew) {
  int nfft=1 << int(log(nlon)/log(2));
  double *Spolecirc;
  double *Npolecirc;
  float Nx=0, Ny=0;		//North pole wind
  float Sx=0, Sy=0;		//South pole wind
  float sinth, costh;

  //Spolecirc=new double[nfft*2];
  //Npolecirc=new double[nfft*2];

  for (int i=0; i<nlev; i++) {
    for (int j=1; j<nlat; j++) {
      for (int k=0; k<nlon; k++) {
        unew[i][j][k]=(u[i][j-1][k]+u[i][j][k])/2;
        vnew[i][j][k]=(v[i][j-1][k]+v[i][j][k])/2;
      }
    }
  }

  //do some shit for the poles:
  for (int i=0; i<nlev; i++) {
    for (int j=0; j<nlon; j++) {
      costh=cos(2*M_PI*j/nlon);
      sinth=sin(2*M_PI*j/nlon);
      //average to a single point at the pole by rotating to a common
      //coordinate system:
      Sx+=u[i][0][j]*costh-v[i][0][j]*sinth;
      Sy+=v[i][0][j]*sinth+u[i][0][j]*costh;
      Nx+=u[i][nlat-1][j]*costh-v[i][nlat-1][j]*sinth;
      Ny+=v[i][nlat-1][j]*sinth+u[i][nlat-1][j]*costh;
    }
    Sx/=nlon;
    Sy/=nlon;
    Nx/=nlon;
    Ny/=nlon;
    //printf("%g %g %g %g\n", Sx, Sy, Nx, Ny);
    for (int j=0; j<nlon; j++) {
      costh=cos(2*M_PI*j/nlon);
      sinth=sin(2*M_PI*j/nlon);
      //printf("%g %g\n", costh, sinth);
      //zonal and meridional winds are just this point in rotated
      //coordinates depeding on the longitude:
      unew[i][0][j]=Sx*costh+Sy*sinth;
      vnew[i][0][j]=Sx*sinth-Sy*costh;
      unew[i][nlat][j]=Nx*costh+Ny*sinth;
      vnew[i][nlat][j]=Nx*sinth-Ny*costh;
      //printf("(%g, %g) ", unew[i][0][j], vnew[i][0][j]);
    }
    //printf("\n");
  }
}

//finds interpolation coefficients for a set of potential temperature levels
float *** pp_interpolate_pt_levels(float *plev, float ***t, 
		int nlev, int nlat, int nlon, 
		float *ptlev, int npt, 
		float pref=1000) {
  float ***coef;
  float pt[nlev];
  long *sind;
  long loc;
  long lastind=-1;

  coef=allocate_3D_field<float>(npt, nlat, nlon);

  for (int i=0; i<nlat; i++) {
    for (int j=0; j<nlon; j++) {
      for (int k=0; k<nlev; k++) {
        pt[k]=t[k][i][j]*pow(pref/plev[k], KAPPA);
	//printf("%g ", pt[k]);
      }
      //printf("\n");
      sind=heapsort(pt, nlev);
      map_vector_inplace(pt, sind, nlev);
      for (int k=0; k<npt; k++) {
        loc=bin_search(pt, nlev, ptlev[k], lastind);
	if (sind[loc]>=nlev-1 || loc < 0) {
          fprintf(stderr, "Potential temperature level (%g) falls out of sky\n", ptlev[k]);
	  throw PARAMETER_OUT_OF_RANGE;
	}
        coef[k][i][j]=sind[(int) loc]+(ptlev[k]-pt[loc])/(pt[loc+1]-pt[loc]);
	//printf("%g ", coef[k][i][j]);
      }
    }
    //printf("\n");
  }

  delete [] sind;

  return coef;
}

float *** pp_zinterpolate(float ***field, float ***coef, int nlev, int nlat, int nlon) {
  float *** fnew;
  fnew=new float**[nlev];
  for (int i=0; i<nlev; i++) {
    fnew[i]=allocate_matrix<float, int32_t>(nlat, nlon);
    for (int j=0; j<nlat; j++) {
      for (int k=0; k<nlon; k++) {
        int pind=(int) coef[i][j][k];
	float frac=coef[i][j][k]-pind;
        fnew[i][j][k]=field[pind][j][k]*(1-frac)
		+field[pind+1][j][k]*frac;
	//printf("%g ", field[pind][j][k]);
      }
      //printf("\n");
    }
  }
  return fnew;
}

#define MAXNREC 1000

int main(int argc, char **argv) {
  char *fname;
  float theta;
  float **data[MAXNREC];
  int32_t *header[MAXNREC];
  int nrec;
  float *lev;
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

  nrec=pp_read_all(fname, header, data, MAXNREC);

  if (flags[1]) {
    printf(" no.       date    code     level  nlat nlon\n");
    for (int i=0; i<nrec; i++) {
      //printf("Record %d:\n", i);
      //field code, vertical level, nlat, nlon, lat0, dlat, lon0, dlon:
        printf("%4d %d/%d/%d-%d %4d %10.3f %4d %4d\n", i, 
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
  
  theta=atof(argv[2]);
  lev=pp_extract_levels(header, nrec, nlev);
  //for (int i=0; i<nlev; i++) printf("%g\n", lev[i]);

  pp_extract_uvwtz(data, header, nrec, u0, v0, w0, t, z);

  nlon=header[0][PP_HEADLOC_NLON];
  nlat=header[0][PP_HEADLOC_NLAT];
  u=allocate_3D_field<float>(nlev, nlat+1, nlon);
  v=allocate_3D_field<float>(nlev, nlat+1, nlon);

  pp_interpolate_uv(u0, v0, nlev, nlat, nlon, u, v);

  if (flags[0]) {
    float ***coef;		//interpolation ceofficients
    coef=pp_interpolate_pt_levels(lev, t, nlev, nlat+1, nlon, &theta, 1);
    uout=pp_zinterpolate(u, coef, 1, nlat+1, nlon);
    vout=pp_zinterpolate(v, coef, 1, nlat+1, nlon);
    delete_3D_field(coef, 1);
  } else {
    double coef;
    int ind;
    double frac;
    reverse(lev, nlev);
    coef=interpolate(lev, nlev, theta);
    if (coef>nlev-1 || coef<0) {
      fprintf(stderr, "Pressure level (%g) out-of-range\n", theta);
      exit(PARAMETER_OUT_OF_RANGE);
    }
    ind=(int) coef;
    frac=coef-ind;
    ind=nlev-ind-1;
    uout=allocate_3D_field<float>(1, nlat+1, nlon);
    vout=allocate_3D_field<float>(1, nlat+1, nlon);
    for (int i=0; i<nlat+1; i++) {
      for (int j=0; j<nlon; j++) {
        uout[0][i][j]=u[ind][i][j]*(1-frac);
        vout[0][i][j]=u[ind][i][j]*(1-frac);
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

  delete_3D_field(u, nlev);
  delete_3D_field(v, nlev);
  delete_3D_field(uout, 1);
  delete_3D_field(vout, 1);

  finish:
    for (int i=0; i<nrec; i++) {
      delete [] header[i];
      delete_matrix(data[i]);
    }

  return 0;

}
