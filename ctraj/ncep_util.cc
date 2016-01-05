#include <math.h>
#include <string.h>

#include <netcdfcpp.h>

#include "error_codes.h"
#include "peteys_tmpl_lib.h"
#include "time_class.h"
#include "simple_temp.h"
#include "dependent_temp.h"

#include "ncep_util.h"

using namespace libpetey;
using namespace datasets;

namespace ctraj {

int get_ncep_dim(NcFile *nci, char *name, simple<float> *&dim) {

  NcVar *ncv;		//grid variables

  float *data;		//for reading the grids
  long ndata;

  //read in the grids for the NCEP data:
  //ncep longitude grid:
  ncv=nci->get_var(name);
  ndata=ncv->num_vals();
  data=new float[ndata];
  ncv->get(data, ndata);
  dim=new simple<float>(data, ndata);
  //for (long i=0;i<nlon;i++) printf("%f ", data[i]);
  //printf("\n");
  delete [] data;

}

//read in the grids from a 4-D ncep file:
int get_ncep_grid(NcFile *nci,			//netcdf file handle 
		simple<float> *&lon, 		//longitude grid
		simple<float> *&lat, 		//latitude grid
		simple<float> *&lev,		//vertical grid
		simple<time_class> *&tgrid){		//time grid
  
  NcVar *ncv;		//grid variables

  float *data;		//for reading the grids
  long ndata;

  //time grids:
  double *traw;
  long nt;
  double tconv;

  time_class *time;

  //read in the grids for the NCEP data:
  //ncep longitude grid:
  ncv=nci->get_var("lon");
  ndata=ncv->num_vals();
  data=new float[ndata];
  ncv->get(data, ndata);
  lon=new simple<float>(data, ndata);
  //for (long i=0;i<nlon;i++) printf("%f ", data[i]);
  //printf("\n");
  delete [] data;

  //add 360 degree longitude grid which wraps 0 degree grid:
  lon->add_el(360);
  
  //ncep latitude grid:
  ncv=nci->get_var("lat");		//some rather unfortunate variable names here...
  ndata=ncv->num_vals();
  data=new float[ndata];
  ncv->get(data, ndata);
  //lat=new simple<float>(data, nlat);
  lat=new simple<float>(data, ndata, 0);
  //printf("\n");
  delete [] data;

  //time grid:
  ncv=nci->get_var("time");
  nt=ncv->num_vals();
  traw=new double[nt];
  ncv->get(traw, nt);

  //if (nt > 366) tconv=HOURSPERDAY; else tconv=1;
  if (nt > 366) tconv=HOURSPERDAY;
  
  //convert time to our format:
  fprintf(stderr, "Converting time grids...\n");
  time=new time_class[nt];
  for (long i=0; i<nt; i++) {
    //convert raw times to units of days:
    traw[i]/=tconv;
    time[i].init(1, 1, 1, 0, 0, 0);
    time[i].add((traw[i])+TOFFS);
  }
  tgrid=new simple<time_class>(time, nt);

  //ncep pressure grid:
  ncv=nci->get_var("level");
  ndata=ncv->num_vals();
  data=new float[ndata];
  ncv->get(data, ndata);
  lev=new simple<float>(data, ndata, 0);
  delete [] data;

  delete [] traw;
  delete [] time;

}

//read in a 3-D field for one time index:
//(new format: single-precision floating point--real simple)
//- doesn't wrap the longitude 
int read_ncep_3d(NcVar *ncv,			//netcdf file handle
		long tind,			//time index
		float *q){			//field

  ind_type dim[3];		//horizontal dimensions
  long n2d;			//size of horizontal slice
  NcDim *ncd;

  for (int i=0; i<3; i++) {
    ncd=ncv->get_dim(i);
    dim[i]=ncd->size();
  }

  n2d=(dim[0])*dim[1];

  ncv->set_cur(tind, 0, 0, 0);
  ncv->get(q, 1, dim[2], dim[1], dim[0]);

}

//read in a 2-D field for one time index at a given vertical level:
//(new format: single-precision floating point--real simple)
int read_ncep_2d(NcVar *ncv,			//netcdf variable handle
		long tind,			//time index
		long zind,			//vertical index
		float *q){			//field

  ind_type dim[2];		//horizontal dimensions
  long n2d;			//size of horizontal slice
  NcDim *ncd;

  for (int i=0; i<2; i++) {
    ncd=ncv->get_dim(i);
    dim[i]=ncd->size();
  }

  n2d=dim[0]*dim[1];

  if (zind==-1) {
    ncv->set_cur(tind, 0, 0);
    ncv->get(q, 1, dim[1], dim[0]);
  } else {
    ncv->set_cur(tind, zind, 0, 0);
    ncv->get(q, 1, 1, dim[1], dim[0]);
  }

}

//read in a 3-D field for one time index:
//(old format: data are stored as short integers)
int read_ncep_3d_old(NcVar *ncv,		//netcdf variable handle
		long tind,			//time index
		float *q){		//field
  NcAtt *att;			//attribute handle...

  double q0, qm;		//conversion from integer to float...

  short int ** ncdata;
  float qval;
  ind_type dim[3];		//horizontal dimensions
  long n2d;			//size of horizontal slice
  NcDim *ncd;

  for (int i=0; i<3; i++) {
    ncd=ncv->get_dim(i);
    dim[i]=ncd->size();
  }

  n2d=dim[0]*dim[1];
  ncdata=new short*[dim[2]];
  ncdata[0]=new short[n2d*dim[2]];
  for (long i=1; i<dim[2]; i++) ncdata[i]=ncdata[0]+i*n2d;

  att=ncv->get_att("scale_factor");
  qm=att->as_double(0);
  att=ncv->get_att("add_offset");
  q0=att->as_double(0);

  ncv->set_cur(tind, 0, 0, 0);
  ncv->get(ncdata[0], 1, dim[2], dim[1], dim[0]);

  for (long k=0; k<n2d; k++) {
    long i=k%(dim[0]);
    long j=k/(dim[0]);

    for (long iz=0; iz<dim[2]; iz++) {
      q[i*n2d+iz]=ncdata[iz][k]*qm+q0;
    }
  }

  delete [] ncdata[0];
  delete [] ncdata;
}

//read in a 2-D field for one time index:
//(old format: data are stored as short integers)
int read_ncep_2d_old(NcVar *ncv,		//netcdf variable handle
		long tind,			//time index
		long zind,			//vertical index
		float *q){			//field
  NcAtt *att;			//attribute handle
  double q0, qm;		//conversion from integer to float...

  short int * ncdata;
  float qval;
  ind_type dim[2];		//horizontal dimensions
  long n2d;			//size of horizontal slice
  NcDim *ncd;

  for (int i=0; i<3; i++) {
    ncd=ncv->get_dim(i);
    dim[i]=ncd->size();
  }

  n2d=(dim[0])*dim[1];
  ncdata=new short[n2d];

  att=ncv->get_att("scale_factor");
  qm=att->as_double(0);
  att=ncv->get_att("add_offset");
  q0=att->as_double(0);

  if (zind==-1) {
    ncv->set_cur(tind, 0, 0);
    ncv->get(ncdata, 1, dim[1], dim[0]);
  } else {
    ncv->set_cur(tind, zind, 0, 0);
    ncv->get(ncdata, 1, 1, dim[1], dim[0]);
  }

  for (long k=0; k<n2d; k++) {
    long i=k%(dim[0]);
    long j=k/(dim[0]);
    q[k]=ncdata[k]*qm+q0;
  }

  delete [] ncdata;
}

//read in a 2-D field at a single vertical level on a single time grid:
int get_ncep_tz(NcFile *nci,			//netcdf file handle
		const char *var, 		//netcdf variabe name
		long tind,			//time index
		double lev,			//z index
		dependent<float> *q){		//field

  NcVar *ncv;			//netcdf variable handle

  long llev;
  double frac;

  float ** ncdata;
  ind_type dim[2];		//horizontal dimensions
  long n2d;			//size of horizontal slice

  double qval;
  float qval1;

  llev=(long) lev;
  frac=lev-llev;

  q->get_dim(dim);
  n2d=(dim[0]-1)*dim[1];
  ncdata=new float*[2];
  ncdata[0]=new float[n2d*2];
  ncdata[1]=ncdata[0]+n2d;

  ncv=nci->get_var(var);

  if (ncv->type() == NC_FLOAT) {
    read_ncep_2d(ncv, tind, llev, ncdata[0]);
    read_ncep_2d(ncv, tind, llev+1, ncdata[1]);
  } else if (ncv->type() == NC_SHORT) {
    read_ncep_2d_old(ncv, tind, llev, ncdata[0]);
    read_ncep_2d_old(ncv, tind, llev+1, ncdata[1]);
  } else {
    fprintf(stderr, "get_ncep_tz: error in file format; exiting ...\n");
    throw FILE_READ_ERROR;
  }

  //printf("qval:", qval);
  for (long k=0; k<n2d; k++) {
    long i=k%(dim[0]-1);
    long j=k/(dim[0]-1);
    j=dim[1]-j-1;		//latitude grids are reversed...

    qval=ncdata[0][k]*(1-frac)+ncdata[1][k]*frac;
    //it's all linear, so the order of operations shouldn't matter:
    //qval=qm*qval+q0;
    q->cel((float) qval, i, j);
    //printf(" %f", qval);
  }
  //printf("\n");

  //wrap the longitude grids:
  for (long j=0; j<dim[1]; j++) {
    q->get(qval1, 0, j);
    q->cel(qval1, dim[0]-1, j);
  }

  delete [] ncdata[0];
  delete [] ncdata;
}

//calculate interpolation coefficients for theta levels:
//theta_level is for a reference pressure of 1.
int get_ncep_theta_interp(NcFile *nc_T,		//temperatures
		float theta_level,		//theta level
		long tind,			//time index
		dependent<double> *c, 		//interpolation coeffs.
		float kappa){			//reference pressure

  NcVar *ncv;
  float **nctdata;

  ind_type dim[2], n2d;		//horizontal dimensions

  float *p;			//pressure levels
  ind_type np;
  float *theta;			//potential temperature profile

  ind_type zloc;
  double frac;

  //get pressure levels:
  ncv=nc_T->get_var("level");
  np=ncv->num_vals();
  p=new float[np];
  theta=new float[np];
  ncv->get(p, np);

  //get horizontal dimensions:
  c->get_dim(dim);
  n2d=(dim[0]-1)*dim[1];

  //allocate memory for raw ncep data:
  nctdata=new float *[np];
  nctdata[0]=new float[np*n2d];
  for (long i=1; i<np; i++) nctdata[i]=nctdata[0]+i*n2d;

  //get 3-D temperature field:
  fprintf(stderr, "%ld %ld %ld\n", dim[0], dim[1], np);
  if (ncv->type() == NC_FLOAT) {
    read_ncep_3d(ncv, tind, nctdata[0]);
  } else if (ncv->type() == NC_SHORT) {
    read_ncep_3d_old(ncv, tind, nctdata[0]);
  } else {
    fprintf(stderr, "get_ncep_theta_interp: error in file format; exiting ...\n");
    throw FILE_READ_ERROR;
  }

  //run through each horizontal grid, calculate theta profile,
  //calculate interpolation coefficient...
  for (ind_type k=0; k<n2d; k++) {
    ind_type i=k%(dim[0]-1);
    ind_type j=k/(dim[0]-1);

    theta[0]=nctdata[0][k]*pow(1./p[0], kappa);
    zloc=-1;
    //printf("Theta levels:\n%f\n", theta[0]);
    for (ind_type zind=1; zind<np; zind++) {
      theta[zind]=nctdata[zind][k]*pow(1./p[zind], kappa);
      //printf("%f\n", theta[zind]);
      if (theta[zind-1] < theta_level && theta_level < theta[zind]) {
        zloc=zind-1;
        break;
      }
    }
    if (zloc == -1) {
      fprintf(stderr, "Theta level, %f, out of bounds [%f, %f]\n", 
		      theta_level, theta[0], theta[np-1]);
      throw PARAMETER_OUT_OF_RANGE;
    }

    frac=(theta_level-theta[zloc])/(theta[zloc+1]-theta[zloc]);

    c->cel(zloc+frac, i, j);

  }

  delete [] nctdata[0];
  delete [] nctdata;

  delete [] p;
  delete [] theta;

}

int get_dthdp(NcFile *nc_T,		//temperatures
		float theta_level,		//theta level
		long tind,			//time index
		dependent<double> *c, 		//interpolation coeffs.
		dependent<float> *dthdp,	//d(theta)/dP
		float kappa){			//reference pressure


  NcVar *ncv;
  NcAtt *att;

  short **nctdata;

  double T0, Tm;

  ind_type dim[2], n2d;		//horizontal dimensions

  float *p;			//pressure levels
  ind_type np;
  float *theta;			//potential temperature profile

  ind_type zloc;
  double frac;

  //get pressure levels:
  ncv=nc_T->get_var("level");
  np=ncv->num_vals();
  p=new float[np];
  theta=new float[np];
  ncv->get(p, np);

  //get horizontal dimensions:
  c->get_dim(dim);
  n2d=(dim[0]-1)*dim[1];

  //allocate memory for raw ncep data:
  nctdata=new short *[np];
  nctdata[0]=new short[np*n2d];
  for (long i=1; i<np; i++) nctdata[i]=nctdata[0]+i*n2d;

  //get 3-D temperature field:
  ncv=nc_T->get_var("air");
  att=ncv->get_att("scale_factor");
  Tm=att->as_double(0);
  att=ncv->get_att("add_offset");
  T0=att->as_double(0);
  ncv->set_cur(tind, 0, 0, 0);
  fprintf(stderr, "%ld %ld %ld\n", dim[0], dim[1], np);
  ncv->get(nctdata[0], 1, np, dim[1], dim[0]-1);

  //run through each horizontal grid, calculate theta profile,
  //calculate interpolation coefficient...
  for (ind_type k=0; k<n2d; k++) {
    ind_type i=k%(dim[0]-1);
    ind_type j=k/(dim[0]-1);

    theta[0]=(nctdata[0][k]*Tm+T0)*pow(1/p[0], kappa);
    zloc=-1;
    //printf("Theta levels:\n%f\n", theta[0]);
    for (ind_type zind=0; zind<np; zind++) {
      theta[zind]=(nctdata[zind][k]*Tm+T0)*pow(1/p[zind], kappa);
      //printf("%f\n", theta[zind]);
      if (theta[zind-1] < theta_level && theta_level < theta[zind]) {
        zloc=zind-1;
	break;
      }
    }
    if (zloc == -1 || zloc == np-1) {
      fprintf(stderr, "Theta level, %f, out of bounds [%f, %f]\n", 
		      theta_level, theta[0], theta[np-1]);
      exit(-2);
    }

    frac=(theta_level-theta[zloc])/(theta[zloc+1]-theta[zloc]);
    c->cel(zloc+frac, i, j);

    //calculate d(theta)/dp: (use a second-order method...)
    float x1, x2, x3;
    float y1, y2, y3;
    float x1_2, x2_2, x3_2;
    float denom;
    float a, b;
    float dydx;

    if (zloc == np-2) {
      x1=theta[zloc-1];
      x2=theta[zloc];
      x3=theta[zloc+1];
      y1=p[zloc-1];
      y2=p[zloc];
      y3=p[zloc+1];
    } else if (zloc == 0) {
      x1=theta[0];
      x2=theta[1];
      x3=theta[2];
      y1=p[0];
      y2=p[1];
      y3=p[2];
    } else if (frac < 0.5) {
      x1=theta[zloc-1];
      x2=theta[zloc];
      x3=theta[zloc+1];
      y1=p[zloc-1];
      y2=p[zloc];
      y3=p[zloc+1];
    } else {
      x1=theta[zloc];
      x2=theta[zloc+1];
      x3=theta[zloc+2];
      y1=p[zloc];
      y2=p[zloc+1];
      y3=p[zloc+2];
    }

    x1_2=x1*x1;
    x2_2=x2*x2;
    x3_2=x3*x3;

    denom=x1*(x3_2-x2_2)-x2*x3_2+x2_2*x3+x1_2*(x2-x3);
    a=(x1*(y3-y2)-x2*y3+x3*y2+(x2-x3)*y1)/denom;
    b=(-x1_2*(x3-y2)-x2_2*y3+x3_2*y2+(x2_2-x3_2)*y1)/denom;

    dydx=2*a*theta_level+b;

    //stick it into the dataset:
    dthdp->cel(1/dydx, i, j);

  }

  delete [] nctdata[0];
  delete [] nctdata;

  delete [] p;
  delete [] theta;

}

//I suspect this will be really slow...
int theta_interp(dependent<float> *q1,
		dependent<double> *c,
		dependent<float> *q2) {

  ind_type dim[3];
  float val;
  double cc;

  q1->get_dim(dim);

  for (ind_type i=0; i<dim[0]; i++) {
    for (ind_type j=0; j<dim[1]; j++) {
      c->get(cc, i, j);
      q1->interpol(val, i, j, cc);
      q2->cel(val, i, j);
    }
  }

}

//this makes a bit more sense:
int get_ncep_theta_level(NcFile *nci,
			const char *var,
			long tind,
			dependent<double> *c,
			dependent<float> *q) {

  NcVar *ncv;			//netcdf variable handle

  float ** ncdata;		//"raw" data
  ind_type dim[2];		//horizontal dimensions
  long n2d;			//size of horizontal slice
  long nz;			//vertical dimension

  float val1;
  double val;

  double cval;			//interpolation coefficients
  long ind;
  double frac;

  //get horizontal dimensions:
  q->get_dim(dim);
  n2d=(dim[0]-1)*dim[1];
  //get vertical dimension:
  ncv=nci->get_var("level");
  nz=ncv->num_vals();

  ncdata=new float*[nz];
  ncdata[0]=new float[n2d*nz];
  for (long i=1; i<nz; i++) ncdata[i]=ncdata[0]+i*n2d;

  ncv=nci->get_var(var);

  if (ncv->type() == NC_FLOAT) {
    read_ncep_3d(ncv, tind, ncdata[0]);
  } else if (ncv->type() == NC_SHORT) {
    read_ncep_3d_old(ncv, tind, ncdata[0]);
  } else {
    fprintf(stderr, "get_ncep_theta_interp: error in file format; exiting ...\n");
    throw FILE_READ_ERROR;
  }

  for (long k=0; k<n2d; k++) {
    long i=k%(dim[0]-1);
    long j=k/(dim[0]-1);
    j=dim[1]-j-1;		//latitude grids are reversed

    c->get(cval, i, j);
    ind=(long) cval;
    //printf("%f\n", cval);
    frac=cval-ind;
    val=ncdata[ind][k]*(1-frac)+ncdata[ind+1][k]*frac;
    q->cel(val, i, j);

  }

  //wrap the longitude grids:
  for (long j=0; j<dim[1]; j++) {
    q->get(val1, 0, j);
    q->cel(val1, dim[0]-1, j);
  }

  delete [] ncdata[0];
  delete [] ncdata;

}

int interpolate_field(dependent<float> *q0,
		dependent<double> **c,
		long ndim,
		dependent<float> *q1) {

  double *int_ind;

  ind_type dim[2];

  float val;

  int_ind=new double[ndim];

  //get the dimensions
  //Assume:
  //1.  coefficients matrices are two-dimensional, as is output dataset
  //2.  all coefficients matrices have same dimensions...)
  c[0]->get_dim(dim);

  for (ind_type i=0; i<dim[0]; i++) {
    for (ind_type j=0; j<dim[1]; j++) {
      for (ind_type k=0; k<ndim; k++) {
        c[k]->get(int_ind[k], i, j);
      }
      q0->interpol(val, int_ind);
      q1->cel(val, i, j);
    }
  }

  delete [] int_ind;

}

int transform_vfield(simple<float> *xgrid, 
		simple<float> *ygrid,
		dependent<float> *u,
		dependent<float> *v,
		int hemi,	
		float vscale) {

  long nx, ny;
  float xval, yval;
  float uval, vval;
  float utval, vtval;
  float r;
  float val;

  nx=xgrid->nel();
  ny=ygrid->nel();

  for (ind_type i=0; i<nx; i++) {
    xgrid->get(xval, i);
    for (ind_type j=0; j<ny; j++) {
      ygrid->get(yval, j);
      u->get(uval, i, j);
      v->get(vval, i, j);
      r=sqrt(xval*xval+yval*yval);
      val=REARTH*sin(r/REARTH);
      uval*=vscale;
      vval*=vscale;
      utval=-hemi*vval*xval/r-uval*yval/val;
      vtval=-hemi*vval*yval/r+uval*xval/val;
      u->cel(utval, i, j);
      v->cel(utval, i, j);
    }
  }
}

simple<time_class> * get_ncep_tgrid(const char *base,	//base filename (incl. path)
			time_class t1, 		//initial time
			time_class t2, 		//final time
			long &it1, 		//index into first year
			long &it2) {		//index into last year
  time_class tdiff;
  time_class sixhour(0, 0, 0, 6, 0, 0);
  long nt;
  int32_t ntyr, ntyrmax;
  int32_t yr1, yr2;
  time_class *t;
  simple<time_class> *tgrid;

  NcVar *ncv;		//grid variables
  NcFile *nc;

  char *fname;
  char fend[10];

  double *traw, tconv;

  yr1=t1.year();
  yr2=t2.year();

  tdiff=t2-t1;
  //maximum number of time grids:
  ntyrmax=(int32_t) (yr2-yr1+1)*366*4;
  //printf("%ld\n", ntyrmax);
  t=new time_class[ntyrmax];

  nt=0;

  //first we get all the time grids and put them in an array:
  fname=new char[strlen(base)+10];
  for (int32_t yr=yr1; yr<=yr2; yr++) {
    sprintf(fname, "%s.%4.4d.nc", base, yr);

    fprintf(stderr, "Opening file: %s\n", fname);
    nc=new NcFile(fname);

    ncv=nc->get_var("time");
    ntyr=ncv->num_vals();
    traw=new double[ntyr];
    ncv->get(traw, ntyr);

    //if (ntyr > 366) tconv=HOURSPERDAY; else tconv=1;
    tconv=HOURSPERDAY;
  
    //convert time to our format:
    fprintf(stderr, "Converting time grids...\n");
    for (int32_t i=0; i<ntyr; i++) {
      //convert raw times to units of days:
      traw[i]/=tconv;
      t[i+nt].init(1, 1, 1, 0, 0, 0);
      t[i+nt].add((traw[i])+TOFFS);
    }

    nt=nt+ntyr;

    delete nc;
    delete [] traw;
  }

  it1=bin_search(t, nt, t1, -1);
  it2=(long) ceil(interpolate(t, nt, t2, -1));

  tgrid=new simple<time_class>(t+it1, it2-it1+1);
  it2=it2-nt+ntyr;
  delete [] t;
  delete [] fname;
  return tgrid;

}

} //end namespace ctraj

