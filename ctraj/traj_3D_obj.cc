#include "assert.h"

//#include "nr.h"

#include "peteys_tmpl_lib.h"
#include "composite_dataset.h"

#include "tcoord_defs.h"
#include "traj_3D_obj.h"


//global datasets containing the velocity field:
simple<float> *global_xgrid;		//Gridding in the x direction
simple<float> *global_ygrid;		//Gridding in the y direction
simple<float> *global_zgrid;		//Gridding in the z direction

dependent_swap<float> **global_u;	//Velocity in the x direction.
dependent_swap<float> **global_v;	//Velocity in the y direction.
dependent_swap<float> **global_w;	//Velocity in the z direction.

dependent_swap<float> *global_sfc;	//Surface boundary...

short int global_zdir;			//z grid ascending or descending?

void traj_3D_derivs(double tind, float xvec[], float dxdt[], long n)
// name: 	traj_derivs
//
// usage: 	derivs(tind, xvec, dxdt, n);
//
// purpose: 	Returns a set of derivatives for use in a Runge-Kutta or
// 		similar integration scheme.  Returns the interpolated values
// 		of a two-dimensional velocity field at the supplied time and
// 		location.
//
// parameters:
// 		tind: (in) The non-dimensional time as a floating point
// 		variable.  Must be normalized to the grid spacing of the
// 		velocity field.
//
// 		xvec: (in) The position.  Must be a two dimensional floating
// 		point array.
//
// 		dxdt: (out) A two-dimensional floating point array containing
// 		the velocity field at the particular time and location.
//
// 		n: (in) Ignored.  Should be two (2).
//
{
  interpol_index xint, yint, zint;
  ind_type zl, z0l;
  interpol_index frac;
  float z0;
  interpol_index z0int;
  float zval;
  float u1, u2, v1, v2, w1, w2;

  xint=global_xgrid->interp(xvec[0]);
  yint=global_ygrid->interp(xvec[1]);
  zint=global_zgrid->interp(xvec[2]);

  global_sfc->interpol(z0, xint, yint, (interpol_index) tind);

  z0int=global_zgrid->interp(z0);

  zl=(ind_type) zint;
  z0l=(ind_type) z0int;
  if (zl == z0l) {
    if (global_zdir > 0) {
      u1=0;
      v1=0;
      w1=0;
      global_u[zl+1]->interpol(u2, xint, yint, 0., (interpol_index) tind);
      global_v[zl+1]->interpol(v2, xint, yint, 0., (interpol_index) tind);
      global_w[zl+1]->interpol(w2, xint, yint, 0., (interpol_index) tind);
      global_zgrid->get(zval, zl);
      frac=(xvec[2]-z0)/(zval-z0);
    } else {
      global_u[zl]->interpol(u1, xint, yint, 0., (interpol_index) tind);
      global_v[zl]->interpol(v1, xint, yint, 0., (interpol_index) tind);
      global_w[zl]->interpol(w1, xint, yint, 0., (interpol_index) tind);
      u2=0;
      v2=0;
      w2=0;
      global_zgrid->get(zval, zl);
      frac=(xvec[2]-z0)/(zval-z0);
    }
  } else {
    global_u[zl]->interpol(u1, xint, yint, 0., (interpol_index) tind);
    global_v[zl]->interpol(v1, xint, yint, 0., (interpol_index) tind);
    global_w[zl]->interpol(w1, xint, yint, 0., (interpol_index) tind);
    global_u[zl+1]->interpol(u2, xint, yint, 0., (interpol_index) tind);
    global_v[zl+1]->interpol(v2, xint, yint, 0., (interpol_index) tind);
    global_w[zl+1]->interpol(w2, xint, yint, 0., (interpol_index) tind);
    frac=zint-(interpol_index) zl;
  }

  dxdt[0]=(1-frac)*u1+frac*u2;
  dxdt[1]=(1-frac)*v1+frac*v2;
  dxdt[2]=(1-frac)*w1+frac*w2;

}

void traj_3D_obj::set_page_size(sub_1d_type page_size) {
  for (int i=0; i<zgrid->nel(); i++) {
    u[i]->set_page_size(page_size, 1);
    v[i]->set_page_size(page_size, 1);
  }
}

void traj_3D_obj::setresize(long nt) {
  //printf("Integrating: res_size=%d, nres=%d\n", res_size, nres);
  if (nt > res_size) {
    for (long i=0; i<=res_size; i++) delete result[i];
    delete result;
    result=new float * [nt+1];
    for (long i=0; i<=nt; i++) result[i]=new float[3];
    res_size=nt;
    delete [] tt;
    tt=new time_class[nt+1];
  }
  nres=nt;

}

void traj_3D_obj::integrate(float p0[3], float pf[3], double tstart, double tstep, long nt) {
  //move the velocity fields to their global counterparts:
  //(not thread safe...)
  global_xgrid=xgrid;
  global_ygrid=ygrid;
  global_zgrid=zgrid;
  global_u=u;
  global_v=v;
  global_w=w;
  global_sfc=sfc;

  //if there are more time steps than can be held in result, re-size it:
  //printf("Integrating: res_size=%d, nres=%d\n", res_size, nres);
  setresize(nt);

  //do the integration:
  rk_dumb(tstart, p0, 3L, tstep, nt, result, &traj_3D_derivs);

  pf[0]=result[nt][0];
  pf[1]=result[nt][1];
  pf[2]=result[nt][2];

  t0=(interpol_index) tstart;
  dt=(interpol_index) tstep;

}

int traj_3D_obj::get_result(time_class date, float p[3]) {
  interpol_index tind;
  interpol_index frac;
  long ind;
  int err;
  
  tind=get_tind(date);
  tind=(tind-t0)/dt;

  ind=(long) tind;
  err=0;
  if (ind < 0) {
    err=-1;
    ind=0;
  }
  if (ind >= res_size) {
    err=1;
    ind=res_size-2;
  }
  frac=tind-(interpol_index) ind;
  p[0]=frac*result[ind+1][0]+(1-frac)*result[ind][0];
  p[1]=frac*result[ind+1][1]+(1-frac)*result[ind][1];
  p[2]=frac*result[ind+1][2]+(1-frac)*result[ind][2];

  return err;
}

time_class traj_3D_obj::get_time(interpol_index tind) {
  time_class tlow, tdiff, thigh, result;
  ind_type lind;
  interpol_index frac;
  char tstring[32];

  lind=(ind_type) tind;
  frac=tind-(interpol_index) lind;

  tgrid->get(tlow, lind);
  tgrid->get(thigh, lind+1);
  tdiff=thigh-tlow;

  tdiff=tdiff*frac;
  result=tlow+tdiff;

  return result;

}

traj_3D_obj::traj_3D_obj(char **filename, simple<float> *z) {
  long loc, dum;
  long nz;

  simple<float> *xdum;
  simple<float> *ydum;
  simple<time_class> *tdum;

  zgrid=z;
  nz=zgrid->nel();
//  composite_dataset vdata;

  vfield_swap=new FILE *[nz+1];
  vdata=new composite_dataset[nz+1];

  //assume x grid, y grid and time grid are all the same...
  vfield_swap[0]=fopen(filename[0], "r");
  
  //get the x grid:
  loc=vdata[0].search_var("xgrid", dum);
  xgrid=(simple<float> *) vdata[0].get_var(loc);
  //get the y grid:
  loc=vdata[0].search_var("ygrid", dum);
  ygrid=(simple<float> *) vdata[0].get_var(loc);
  //get the time grid:
  loc=vdata[0].search_var("time", dum);
  tgrid=(simple<time_class> *) vdata[0].get_var(loc);

  //get the velocity in the x direction:
  loc=vdata[0].search_var("u", dum);
  u[0]=(dependent_swap<float> *) vdata[0].get_var(loc);
  //get the velocity in the y direction:
  loc=vdata[0].search_var("v", dum);
  v[0]=(dependent_swap<float> *) vdata[0].get_var(loc);
  //get the velocity in the z direction:
  loc=vdata[0].search_var("w", dum);
  w[0]=(dependent_swap<float> *) vdata[0].get_var(loc);

  for (long i=1; i<nz; i++) {
    //read in the velocity field:
    vfield_swap[i]=fopen(filename[i], "r");
    vdata[i].read(vfield_swap[i]);

    //get the x grid:
    loc=vdata[i].search_var("xgrid", dum);
    xdum=(simple<float> *) vdata[i].get_var(loc);
    assert(xdum==xgrid);
    //get the y grid:
    loc=vdata[i].search_var("ygrid", dum);
    ydum=(simple<float> *) vdata[i].get_var(loc);
    assert(ydum==ygrid);
    //get the time grid:
    loc=vdata[i].search_var("time", dum);
    tdum=(simple<time_class> *) vdata[i].get_var(loc);
    assert(tdum==tgrid);

    //get the velocity in the x direction:
    loc=vdata[i].search_var("u", dum);
    u[i]=(dependent_swap<float> *) vdata[i].get_var(loc);
    //get the velocity in the y direction:
    loc=vdata[i].search_var("v", dum);
    v[i]=(dependent_swap<float> *) vdata[i].get_var(loc);
    //get the velocity in the z direction:
    loc=vdata[i].search_var("w", dum);
    w[i]=(dependent_swap<float> *) vdata[i].get_var(loc);

    //tgrid->print();
  }

  //read in the velocity field:
  vfield_swap[nz]=fopen(filename[nz], "r");
  vdata[nz].read(vfield_swap[nz]);

  //get the x grid:
  loc=vdata[nz].search_var("xgrid", dum);
  xdum=(simple<float> *) vdata[nz].get_var(loc);
  assert(xdum==xgrid);
  //get the y grid:
  loc=vdata[nz].search_var("ygrid", dum);
  ydum=(simple<float> *) vdata[nz].get_var(loc);
  assert(ydum==ygrid);
  //get the time grid:
  loc=vdata[nz].search_var("time", dum);
  tdum=(simple<time_class> *) vdata[nz].get_var(loc);
  assert(tdum==tgrid);

  //get the surface coordinate:
  loc=vdata[nz].search_var("sfc", dum);
  sfc=(dependent_swap<float> *) vdata[nz].get_var(loc);

  nres=0;
  res_size=0;
  result=new float * [1];
  result[0]=new float[2];

  tt=new time_class[1];

}

long traj_3D_obj::print_result(FILE *fs) {
  time_class tcur;
  char tstr[MAX_TSTR_WIDTH];

  for (long i=0; i<=nres; i++) {
    tcur=get_time(t0+dt*i);
    tcur.write_string(tstr);
    fprintf(fs, "%s %10.6g %10.6g %10.6g\n", tstr, result[i][0], result[i][1], result[i][2]);
  }

  return nres;
}

long traj_3D_obj::print_result(FILE *fs, short hemi) {
  time_class tcur;
  char tstr[MAX_TSTR_WIDTH];
  float lon, lat;

  for (long i=1; i<=nres; i++) {
    tcur=get_time(t0+dt*i);
    tcur.write_string(tstr);
    tcoord2_2lonlat(result[i][0], result[i][1], 1, hemi, lon, lat);
    fprintf(fs, "%s %10.6g %10.6g %10.6\n", tstr, lon, lat, result[i][2]);
  }

  return nres;

}

traj_3D_obj::~traj_3D_obj() {

  for (long i=0; i<=nres; i++) delete result[i];
  delete result;

//  delete u;
//  delete v;
  
  for (long i=0; i<=zgrid->nel(); i++) {
    fclose(vfield_swap[i]);
  }

  delete [] vdata;
  delete [] vfield_swap;

//  delete xgrid;
//  delete ygrid;
//  delete tgrid;

}
