#include <stdio.h>
#include <stdint.h>

#include <gsl/gsl_fft_complex.h>

#include "peteys_tmpl_lib.h"
#include "error_codes.h"
#include "full_util.h"
#include "time_class.h"

#include "pp_util.h"

#include "ctraj_defaults.h"

using namespace libpetey;

namespace ctraj {

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
    field=new scalar **[nlev];
    for (int i=0; i<nlev; i++) field[i]=allocate_matrix<scalar, int32_t>(ny, nx);
    return field;
  }

  template <typename scalar>
  void delete_3D_field(scalar ***field, int nlev) {
    for (int i=0; i<nlev; i++) delete_matrix(field[i]);
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
      if (f77recsep != PP_HEADLEN*4) {
        fprintf(stderr, "Error: error in record separator--wrong file type\n");
        fprintf(stderr, " actual: %d; expected: %d\n", f77recsep, PP_HEADLEN*4);
        throw FILE_READ_ERROR;
      }

      header=new int32_t[PP_HEADLEN];
      fread(header, sizeof(int32_t), PP_HEADLEN, fs);
      swap_endian(header, PP_HEADLEN);
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

  int pp_read_field(FILE *fs, int32_t **headers, float ***fields, int nmax, int field_code, float *plev, int nlev, int toendflag) {

    int32_t f77recsep;
    int32_t *header;
    int nhead_check;		//check number of elements read in
    float **data;
    int readcount;
    int nrec;
    int n, n1;		//number of elements read in should match
    int nlon, nlat;
    float level;
    int lastcode=-1;
    int code=-1;

    do {
      readcount=fread(&f77recsep, sizeof(f77recsep), 1, fs);
      //printf("%d %d\n", nrec, readcount);
      if (readcount==0) {
        break;
      }
      swap_endian(&f77recsep, 1);

      //should be 256:
      //printf("%d\n", f77recsep);
      if (f77recsep != PP_HEADLEN*sizeof(int32_t)) {
        fprintf(stderr, "Error: error in record separator--wrong file type\n");
        fprintf(stderr, " actual: %d; expected: %d\n", f77recsep, PP_HEADLEN*4);
        throw FILE_READ_ERROR;
      }

      header=new int32_t[PP_HEADLEN];
      nhead_check=fread(header, sizeof(int32_t), PP_HEADLEN, fs);
      if (nhead_check!=PP_HEADLEN) {
        fprintf(stderr, "Not enough header records read in: %d vs. %d\n", nhead_check, PP_HEADLEN);
        if (nhead_check==0) break; else throw FILE_READ_ERROR;
      }

      swap_endian(header, PP_HEADLEN);
      headers[nrec]=header;

      fread(&f77recsep, sizeof(f77recsep), 1, fs);
      fread(&f77recsep, sizeof(f77recsep), 1, fs);
      swap_endian(&f77recsep, 1);

      //should be 6912*4:
      //printf("%d\n", f77recsep);
      n=header[PP_HEADLOC_N];
      nlat=header[PP_HEADLOC_NLAT];
      nlon=header[PP_HEADLOC_NLON];
      lastcode=code;
      code=header[PP_HEADLOC_CODE];
      level=*(float *)(header+PP_HEADLOC_LEV);

      if (n!=nlon*nlat) {
        fprintf(stderr, "Error in record header: %d*%d != %d\n", nlat, nlon, n);
        throw FILE_READ_ERROR;
      }
      if (n*sizeof(int32_t)!=f77recsep) {
        fprintf(stderr, "Error in record header: %d != %d\n", header[PP_HEADLOC_NLAT], f77recsep);
        throw FILE_READ_ERROR;
      }

      data=NULL;
      if (code == field_code) {
        if (plev!=NULL) {
          for (int i=0; i<nlev; i++) {
            if (plev[i]==level) {
              data=allocate_matrix<float, int32_t>(nlat, nlon);
	      break;
	    }
	  }
	} else {
          data=allocate_matrix<float, int32_t>(nlat, nlon);
	}
      }

      if (data!=NULL) {
        n1=fread(data[0], sizeof(float), n, fs);
        if (n1!=n) {
          fprintf(stderr, "Not enough header records read in: %d vs. %d\n", nhead_check, PP_HEADLEN);
          throw FILE_READ_ERROR;
        }
        swap_endian((int32_t *) data[0], n);
        fields[nrec]=data;
	nrec++;
        fread(&f77recsep, sizeof(f77recsep), 1, fs);
      } else {
        fseek(fs, n*sizeof(int32_t)+1, SEEK_CUR);
      }
    } while (lastcode!=field_code || lastcode==code || toendflag);

    return nrec;

  }

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
		float pref) {
    float ***coef;
    float pt[nlev];
    long sind[nlev];
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
        heapsort(pt, sind, nlev);
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

  template float ***allocate_3D_field<float>(int nlev, int ny, int nx);
  template double ***allocate_3D_field<double>(int nlev, int ny, int nx);

  template void delete_3D_field<float>(float ***field, int nlev);
  template void delete_3D_field<double>(double ***field, int nlev);

} //end namespace ctraj
