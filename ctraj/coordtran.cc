
#include <math.h>
#include <assert.h>

#include "tcoord_defs.h"
#include "coordtran.h"
#include "ctraj_defaults.h"

using namespace libpetey::datasets;
namespace ctraj {

simple<float> * create_grid(ind_type ngrid,
		float slen) {

  float grid[ngrid];
  for (ind_type i=0; i<ngrid; i++) {
    grid[i]=i*slen*2/(ngrid-1)-slen;
  }

  return new simple<float>(grid, ngrid);
}    

//given a pair of grids, calculate interpolation coefficients
//to go from lon-lat to azimuthal equidistant...
void intcoeff(simple<float> *lon,
		simple<float> *lat,
		simple<float> *xgrid,
		simple<float> *ygrid,
		dependent<double> *c1, 
		dependent<double> *c2,
		short hemi) {

  double r;
  double c1val;
  double c2val;
  float lonp;
  float latp;
  float xval;
  float yval;

  ind_type nx, ny;

  //calculate the interpolation coefficients:
  //c1=new dependent<double>(xgrid, ygrid);
  //c2=new dependent<double>(xgrid, ygrid);
 
  nx=xgrid->nel();
  ny=ygrid->nel();
  
  //printf("%d %d\n", nx, ny);
  for (ind_type i=0; i<nx; i++) {
    for (ind_type j=0; j<ny; j++) {
      xgrid->get(xval, i);
      ygrid->get(yval, j);
      //printf("%f %f\n", xval, yval);
      //transform to lon-lat coords:

      //printf("|(%f, %f)|=%f\n", xval, yval, r);
      //printf("(%f, %f)\n", lonp, latp);
      
      tcoord2_2lonlat(xval, yval, 1, hemi, lonp, latp);
      //printf("%f %f\n", lonp, latp);

      //get the interpolation coefficients:
      c1val=lon->interp(lonp);
      c2val=lat->interp(latp);
      //printf("(%f, %f)\n", c1val, c2val);
      c1->cel(c1val, i, j);
      c2->cel(c2val, i, j);

    }
  }
}

//since I haven't figured out a way to do this analytically...
sub_1d_type calc_nmap(ind_type ngrid) {

  simple<float> *xgrid;
  simple<float> *ygrid;
  sub_1d_type nmap;
  float r=SIDELENGTH_Q;
  dependent<sub_1d_type> *map;

  xgrid=new simple<float>(-r, r, ngrid);
  ygrid=new simple<float>(-r, r, ngrid);

  map=new dependent<sub_1d_type>(xgrid, ygrid);

  nocorner_map(xgrid, ygrid, map, nmap);

  delete map;
  delete xgrid;
  delete ygrid;

  return nmap;
}

//given a pair of grids, calculate interpolation coefficients
//to go from lon-lat to azimuthal equidistant...
void nocorner_map(simple<float> *xgrid,
		simple<float> *ygrid,
		dependent<sub_1d_type> *map, sub_1d_type &k, 
		az_eq_t<float> *m) {

  double r;
  float xval;
  float yval;

  ind_type nx, ny;

  float dr;
  float dum;
  float dx, dy;

  float maxr;

  //calculate the interpolation coefficients:
  //c1=new dependent<double>(xgrid, ygrid);
  //c2=new dependent<double>(xgrid, ygrid);
 
  nx=xgrid->nel();
  ny=ygrid->nel();

  xgrid->get(dx, nx-1);
  xgrid->get(dum, (ind_type)0);
  dx=(dx-dum)/(nx-1);
  ygrid->get(dy, ny-1);
  ygrid->get(dum, (ind_type)0);
  dy=(dy-dum)/(ny-1);

  dr=sqrt(dx*dx+dy*dy);

  if (m!=NULL) maxr=m->getr()*M_PI/2;  else maxr=REARTH*M_PI/2;
  maxr+=dr;
  //printf("maxr=%g + %g\n", maxr, dr);

  //printf("%d %d\n", nx, ny);
  k=0;
  for (ind_type j=0; j<ny; j++) {
    for (ind_type i=0; i<nx; i++) {
      xgrid->get(xval, i);
      ygrid->get(yval, j);
      //printf("%f %f\n", xval, yval);
      //transform to lon-lat coords:
      r=sqrt(xval*xval+yval*yval);

      //create the mapping to remove the corners:
      if (r > maxr) {
        map->cel(-1, i, j);
      } else {
        map->cel(k, i, j);
	k++;
      }
    }
  }
}

//given a pair of grids, calculate interpolation coefficients
//to go from azimuthal equidistant to lon-lat...
void intcoeff2(simple<float> *lon,
		simple<float> *lat,
		simple<float> *xgrid,
		simple<float> *ygrid,
		dependent<double> *c1, 
		dependent<double> *c2,
		short hemi) {

  double r;
  double c1val;
  double c2val;
  float lonp;
  float latp;
  float xval;
  float yval;

  ind_type nlon, nlat;

  //calculate the interpolation coefficients:
  //c1=new dependent<double>(lon, lat);
  //c2=new dependent<double>(lon, lat);
 
  nlon=lon->nel();
  nlat=lat->nel();
  
  //printf("%d %d\n", nx, ny);
  for (ind_type i=0; i<nlon; i++) {
    for (ind_type j=0; j<nlat; j++) {
      lon->get(lonp, i);
      lat->get(latp, j);
      //printf("%f %f\n", xval, yval);
      //transform to lon-lat coords:
      tcoord2_2lonlat(xval, yval, -1, hemi, lonp, latp);

      //printf("|(%f, %f)|=%f\n", xval, yval, r);
      //printf("(%f, %f)\n", lonp, latp);

      //get the interpolation coefficients:
      c1val=xgrid->interp(xval);
      c2val=ygrid->interp(yval);
      //printf("(%f, %f)\n", c1val, c2val);
      c1->cel(c1val, i, j);
      c2->cel(c2val, i, j);

    }
  }
}

//interpolate a field:
void intq(dependent<float> *q0,
		dependent<double> *c1,
		dependent<double> *c2,
		dependent<float> *q1) {

  ind_type dim[2];
  double c1val, c2val;
  float val;

  assert(q0->get_rank() == 2);
  assert(q1->get_rank() == 2);

  q0->get_dim(dim);

  for (ind_type i=0; i<dim[0]; i++) {
    for (ind_type j=0; j<dim[1]; j++) {
      //do the interpolation:
      c1->get(c1val, i, j);
      c2->get(c2val, i, j);
      q1->interpol(val, c1val, c2val);
    }
  }

}

//interpolate a velocity field (includes metric transformations):

} //end namespace ctraj


