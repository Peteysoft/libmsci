//#include "nr.h"

#include "peteys_tmpl_lib.h"
#include "composite_dataset.h"

#include "tcoord_defs.h"
#include "traj_int_obj.h"

/*
//global datasets containing the velocity field:
simple<float> *global_xgrid;		//Gridding in the x direction
simple<float> *global_ygrid;		//Gridding in the y direction

dependent_swap<float> *global_u;	//Velocity in the x direction.
dependent_swap<float> *global_v;	//Velocity in the y direction.
*/

int traj_derivs(double tind, float xvec[], float dxdt[], void *param)
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
  interpol_index xint, yint;

  simple<float> *xgrid;		//Gridding in the x direction
  simple<float> *ygrid;		//Gridding in the y direction

  dependent_swap<float> *u;	//Velocity in the x direction.
  dependent_swap<float> *v;	//Velocity in the y direction.

  dataset **param2=(dataset **) param;

  xgrid=(simple<float> *) param2[0];
  ygrid=(simple<float> *) param2[1];
  u=(dependent_swap<float> *) param2[2];
  v=(dependent_swap<float> *) param2[3];

  xint=xgrid->interp(xvec[0]);
  yint=ygrid->interp(xvec[1]);
  u->interpol(dxdt[0], xint, yint, 0., (interpol_index) tind);
  v->interpol(dxdt[1], xint, yint, 0., (interpol_index) tind);

}

void traj_int_obj::set_page_size(sub_1d_type page_size) {
  u->set_page_size(page_size, 1);
  v->set_page_size(page_size, 1);
}

void traj_int_obj::setresize(long nt) {
  //printf("Integrating: res_size=%d, nres=%d\n", res_size, nres);
  if (nt > res_size) {
    for (long i=0; i<=res_size; i++) delete result[i];
    delete [] result;
    result=new float * [nt+1];
    for (long i=0; i<=nt; i++) result[i]=new float[2];
    res_size=nt;
    delete [] tt;
    tt=new time_class[nt+1];
  }
  nres=nt;

}

void traj_int_obj::integrate(float p0[2], float pf[2], double tstart, double tstep, long nt) {
  /*
  //move the velocity fields to their global counterparts:
  global_xgrid=xgrid;
  global_ygrid=ygrid;
  global_u=u;
  global_v=v;
  */

  //if there are more time steps than can be held in result, re-size it:
  //printf("Integrating: res_size=%d, nres=%d\n", res_size, nres);
  setresize(nt);

  //do the integration:
  rk_dumb(tstart, p0, 2L, tstep, nt, result, (void *) param, &traj_derivs);

  pf[0]=result[nt][0];
  pf[1]=result[nt][1];

  t0=(interpol_index) tstart;
  dt=(interpol_index) tstep;

}

void traj_int_obj::integrate(float p0[2], float **result1, double tstart, double tstep, long nt) {
  //move the velocity fields to their global counterparts:
/*
  global_xgrid=xgrid;
  global_ygrid=ygrid;
  global_u=u;
  global_v=v;
*/

  //do the integration:
  rk_dumb(tstart, p0, 2L, tstep, nt, result1, (void *) param, &traj_derivs);

  t0=(interpol_index) tstart;
  dt=(interpol_index) tstep;

}

int traj_int_obj::get_result(time_class date, float p[2]) {
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

  return err;
}

float ** traj_int_obj::int_tarray(float p0[2], time_class tstart, time_class *tarr, long nt) {

  double tind[nt];
  double tind0;
  long ind0;

  /*
  //move the velocity fields to their global counterparts:
  global_xgrid=xgrid;
  global_ygrid=ygrid;
  global_u=u;
  global_v=v;
  */

  setresize(nt-1);

  for (long i=0; i<nt; i++) tind[i]=get_tind(tarr[i]);
  tind0=get_tind(tstart);
  ind0=bin_search(tind, nt, tind0, -1);

  if (ind0 < 0) ind0=0;
  if (ind0 >= nt) ind0=nt;

  //integrate forward:
  if (tind0 < tind[nt-1]) {
    rk_dumb(tind0, p0, 2, tind[ind0+1]-tind0, 1, result+ind0, (void *) param, &traj_derivs);
    for (long i=ind0+1; i<nt-1; i++) {
      rk_dumb(tind[i], result[i], 2, tind[i+1]-tind[i], 1, result+i, (void *) param, &traj_derivs);
    }
  }

  //integrate backwards:
  if (tind0 > tind[0]) {
    float **dum;
    dum=new float *[2];
    dum[0]=result[ind0];
    dum[1]=result[ind0-1];
    rk_dumb(tind0, p0, 2, tind[ind0]-tind0, 1, dum, (void *) param, &traj_derivs);
    for (long i=ind0; i>0; i--) {
      dum[0]=result[i];
      dum[1]=result[i-1];
      rk_dumb(tind[i], result[i], 2, tind[i-1]-tind[i], 1, dum, (void *) param, &traj_derivs);
    }
    delete [] dum;
  }


  return result;

}

/*
float ** traj_int_obj::int_tarray_approx(float p0[2], time_class tstart, time_class *tarr, long nt, double toffs) {

  float t[3];
  float tind0;
  float x[3], y[3];		//we interpolate between these three points
  float **intres;
  float tind;

  //double t[3];
  //double tind0;
  //double x[3], y[3];		//we interpolate between these three points
  //double tind;

  //polynomial coefficients:
  float cofx[3];
  float cofy[3];

  //double cofx[3];
  //double cofy[3];

  //move the velocity fields to their global counterparts:
  global_xgrid=xgrid;
  global_ygrid=ygrid;
  global_u=u;
  global_v=v;

  setresize(nt-1);

  //select three points, integrate between them and do a 2nd order 
  //interpolation between them:
  tind0=(float) get_tind(tstart);
  t[1]=0.;
  t[0]=t[1]-toffs;
  t[2]=t[1]+toffs;
  x[1]=0.;
  y[1]=0.;
  intres=rk_dumb((double) tind0, p0, 2, (double) toffs, 1, &traj_derivs);
  x[2]=intres[1][0]-p0[0];
  y[2]=intres[1][1]-p0[1];
  delete [] intres[0];
  delete [] intres[1];
  delete [] intres;
  intres=rk_dumb((double) tind0, p0, 2, -(double) toffs, 1, &traj_derivs);
  x[0]=intres[1][0]-p0[0];
  y[0]=intres[1][1]-p0[1];
  delete [] intres[0];
  delete [] intres[1];
  delete [] intres;

  polcoe(t, x, 2, cofx);
  polcoe(t, y, 2, cofy);

//  printf("x0= %f; y0= %f\n", p0[0], p0[1]);
//  printf("x_0= %f; y_0= %f\n", x[0]+p0[0], y[0]+p0[1]);

  for (long i=0; i<nt; i++) {
    tind=get_tind(tarr[i])-tind0;
//    printf("tind=%f\n", tind);
    result[i][0]=cofx[0]+tind*(cofx[1]+tind*cofx[2])+p0[0];
    result[i][1]=cofy[0]+tind*(cofy[1]+tind*cofy[2])+p0[1];
  }

  return result;

}
*/

time_class traj_int_obj::get_time(interpol_index tind) {
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

traj_int_obj::traj_int_obj(char *filename) {
  long loc, dum;
//  composite_dataset vdata;

  //read in the velocity field:
  vfield_swap=fopen(filename, "r");
  vdata.read(vfield_swap);

  //get the x grid:
  loc=vdata.search_var("xgrid", dum);
  xgrid=(simple<float> *) vdata.get_var(loc);
  param[0]=xgrid;
  //get the y grid:
  loc=vdata.search_var("ygrid", dum);
  ygrid=(simple<float> *) vdata.get_var(loc);
  param[1]=ygrid;
  //get the time grid:
  loc=vdata.search_var("time", dum);
  tgrid=(simple<time_class> *) vdata.get_var(loc);
  //get the velocity in the x direction:
  loc=vdata.search_var("u", dum);
  u=(dependent_swap<float> *) vdata.get_var(loc);
  param[2]=u;
  //get the velocity in the y direction:
  loc=vdata.search_var("v", dum);
  v=(dependent_swap<float> *) vdata.get_var(loc);
  param[3]=v;

  //tgrid->print();

  nres=0;
  res_size=0;
  result=new float * [1];
  result[0]=new float[2];

  tt=new time_class[1];

}

ind_type traj_int_obj::nt() {
  return tgrid->nel();
}

void traj_int_obj::get_Jmatrix(float x, float y, interpol_index tind, float *j) {
  //interpolation indices for x and y:
  interpol_index xind, yind;
  //winds at all four corners:
  float u1, u2, u11, u12, u21, u22;
  float v1, v2, v11, v12, v21, v22;
  //x and y coords at all four corners:
  float x1, x2;
  float y1, y2;
  //average values of x and y:
  //change in x and y:
  float dx, dy;
  //integer indices:
  ind_type lx, ly, lt;
  interpol_index frac;

  //metric coefficients:
  float cx, cy;
  //for calculating the metric coefficients:
  float r2, r;
  float xp2, yp2;
  float sinRor, sin2Ror;

  //for the time interval:
  time_class t1, t2;
  double delt;

  xind=xgrid->interp(x);
  yind=ygrid->interp(y);

  lx=(ind_type) xind;
  ly=(ind_type) yind;
  lt=(ind_type) tind;

  xgrid->get(x1, lx);
  xgrid->get(x2, lx+1);
  ygrid->get(y1, ly);
  ygrid->get(y2, ly+1);

  //calculate metric coefficients:
  xp2=x*x;
  yp2=y*y;
  r2=xp2+yp2;
  r=sqrt(r2);
  sinRor=sin(REARTH/r);
  sin2Ror=sinRor*sinRor;
  cx=1./r2*(REARTH*REARTH/r2*sin2Ror*yp2+xp2);
  cx=sqrt(cx);
  cy=1./r2*(REARTH*REARTH/r2*sin2Ror*xp2+yp2);
  cy=sqrt(cy);

  dx=(x2-x1)*cx;
  dy=(y2-y1)*cy;

  //printf("lx=%d, ly=%d, lt=%d\n", lx, ly, lt);

  //get the winds at all four corners:
  frac=tind-(interpol_index) lt;
  u->get(u1, lx, ly, 0, lt);
  u->get(u2, lx, ly, 0, lt+1);
  //printf("u1=%f, u2=%f\n", u1, u2);
  u11=u1*(1-frac)+u2*frac;
  u->get(u1, lx, ly+1, 0, lt);
  u->get(u2, lx, ly+1, 0, lt+1);
  //printf("u1=%f, u2=%f\n", u1, u2);
  u12=u1*(1-frac)+u2*frac;
  u->get(u1, lx+1, ly, 0, lt);
  u->get(u2, lx+1, ly, 0, lt+1);
  u21=u1*(1-frac)+u2*frac;
  u->get(u1, lx+1, ly+1, 0, lt);
  u->get(u2, lx+1, ly+1, 0, lt+1);
  u22=u1*(1-frac)+u2*frac;
  //printf("u11=%f, u12=%f, u21=%f, u22=%f\n", u11, u12, u12, u22);

  v->get(v1, lx, ly, 0, lt);
  v->get(v2, lx, ly, 0, lt+1);
  v11=v1*(1-frac)+v2*frac;
  v->get(v1, lx, ly+1, 0, lt);
  v->get(v2, lx, ly+1, 0, lt+1);
  v12=v1*(1-frac)+v2*frac;
  v->get(v1, lx+1, ly, 0, lt);
  v->get(v2, lx+1, ly, 0, lt+1);
  v21=v1*(1-frac)+v2*frac;
  v->get(v1, lx+1, ly+1, 0, lt);
  v->get(v2, lx+1, ly+1, 0, lt+1);
  v22=v1*(1-frac)+v2*frac;

  frac=xind-(interpol_index) lx;
  j[1]=((u12-u11)*(1-frac)+(u22-u21)*frac)/dy;
  j[3]=((v12-v11)*(1-frac)+(v22-v21)*frac)/dy;

  frac=yind-(interpol_index) ly;
  j[0]=((u21-u11)*(1-frac)+(u22-u12)*frac)/dx;
  j[2]=((v21-v11)*(1-frac)+(v22-v12)*frac)/dx;

/*
  tgrid->get(t1, lt);
  tgrid->get(t2, lt+1);
  delt=(double) t2 - (double) t1;
  delt*=86400;
*/

  j[0]*=cx;
  j[1]*=cy;
  j[2]*=cx;
  j[3]*=cy;

}

void traj_int_obj::integrate_Hmatrix(float *hmatrix) {
  float jmatrix[4];
  float dh[4];
  interpol_index tind;

  //take an Euler step:
  get_Jmatrix((result[0][0]+result[1][0])/2, 
		(result[0][1]+result[1][1])/2, tind, jmatrix);
  printf("J=%f %f %f %f\n", jmatrix[0], jmatrix[1], jmatrix[2], jmatrix[3]);

  dh[0]=jmatrix[0]*hmatrix[0]+jmatrix[1]*hmatrix[2];
  dh[1]=jmatrix[0]*hmatrix[1]+jmatrix[1]*hmatrix[3];
  dh[2]=jmatrix[2]*hmatrix[0]+jmatrix[3]*hmatrix[2];
  dh[3]=jmatrix[2]*hmatrix[1]+jmatrix[3]*hmatrix[3];

  hmatrix[0]+=dh[0]*dt/2;
  hmatrix[1]+=dh[1]*dt/2;
  hmatrix[2]+=dh[2]*dt/2;
  hmatrix[3]+=dh[3]*dt/2;

  for (long i=0; i<nres; i++) {
    tind=t0+dt*i;
    get_Jmatrix(result[i][0], result[i][1], tind, jmatrix);
    dh[0]=jmatrix[0]*hmatrix[0]+jmatrix[1]*hmatrix[2];
    dh[1]=jmatrix[0]*hmatrix[1]+jmatrix[1]*hmatrix[3];
    dh[2]=jmatrix[2]*hmatrix[0]+jmatrix[3]*hmatrix[2];
    dh[3]=jmatrix[2]*hmatrix[1]+jmatrix[3]*hmatrix[3];

    hmatrix[0]+=dh[0]*dt;
    hmatrix[1]+=dh[1]*dt;
    hmatrix[2]+=dh[2]*dt;
    hmatrix[3]+=dh[3]*dt;

  }

  //take a final Euler step:
  get_Jmatrix(result[nres][0], result[nres][1], tind, jmatrix);
  dh[0]=jmatrix[0]*hmatrix[0]+jmatrix[1]*hmatrix[2];
  dh[1]=jmatrix[0]*hmatrix[1]+jmatrix[1]*hmatrix[3];
  dh[2]=jmatrix[2]*hmatrix[0]+jmatrix[3]*hmatrix[2];
  dh[3]=jmatrix[2]*hmatrix[1]+jmatrix[3]*hmatrix[3];

  hmatrix[0]+=dh[0]*dt/2;
  hmatrix[1]+=dh[1]*dt/2;
  hmatrix[2]+=dh[2]*dt/2;
  hmatrix[3]+=dh[3]*dt/2;

}

void traj_int_obj::integrate_Hcomp_matrix(float *hmatrix) {
  float jmatrix[4];
  float dh[4];
  interpol_index tind;

  //take an Euler step:
  get_Jmatrix((result[0][0]+result[1][0])/2, 
		(result[0][1]+result[1][1])/2, tind, jmatrix);
  dh[0]=hmatrix[0]*jmatrix[0]+hmatrix[1]*jmatrix[2];
  dh[1]=hmatrix[0]*jmatrix[1]+hmatrix[1]*jmatrix[3];
  dh[2]=hmatrix[2]*jmatrix[0]+hmatrix[3]*jmatrix[2];
  dh[3]=hmatrix[2]*jmatrix[1]+hmatrix[3]*jmatrix[3];

  hmatrix[0]-=dh[0]*dt/2;
  hmatrix[1]-=dh[1]*dt/2;
  hmatrix[2]-=dh[2]*dt/2;
  hmatrix[3]-=dh[3]*dt/2;

  for (long i=0; i<nres; i++) {
    tind=t0+dt*i;
    get_Jmatrix(result[i][0], result[i][1], tind, jmatrix);
    dh[0]=hmatrix[0]*jmatrix[0]+hmatrix[1]*jmatrix[2];
    dh[1]=hmatrix[0]*jmatrix[1]+hmatrix[1]*jmatrix[3];
    dh[2]=hmatrix[2]*jmatrix[0]+hmatrix[3]*jmatrix[2];
    dh[3]=hmatrix[2]*jmatrix[1]+hmatrix[3]*jmatrix[3];

    hmatrix[0]-=dh[0]*dt;
    hmatrix[1]-=dh[1]*dt;
    hmatrix[2]-=dh[2]*dt;
    hmatrix[3]-=dh[3]*dt;

  }

  //take a final Euler step:
  get_Jmatrix(result[nres][0], result[nres][1], tind, jmatrix);
  dh[0]=hmatrix[0]*jmatrix[0]+hmatrix[1]*jmatrix[2];
  dh[1]=hmatrix[0]*jmatrix[1]+hmatrix[1]*jmatrix[3];
  dh[2]=hmatrix[2]*jmatrix[0]+hmatrix[3]*jmatrix[2];
  dh[3]=hmatrix[2]*jmatrix[1]+hmatrix[3]*jmatrix[3];

  hmatrix[0]-=dh[0]*dt/2;
  hmatrix[1]-=dh[1]*dt/2;
  hmatrix[2]-=dh[2]*dt/2;
  hmatrix[3]-=dh[3]*dt/2;

}

void traj_int_obj::integrate_Hmatrix(double *hmatrix) {
  float jmatrix[4];
  float dh[4];
  interpol_index tind;

  //take an Euler step:
  get_Jmatrix((result[0][0]+result[1][0])/2, 
		(result[0][1]+result[1][1])/2, tind, jmatrix);
  //printf("J=%f %f %f %f\n", jmatrix[0], jmatrix[1], jmatrix[2], jmatrix[3]);

  dh[0]=jmatrix[0]*hmatrix[0]+jmatrix[1]*hmatrix[2];
  dh[1]=jmatrix[0]*hmatrix[1]+jmatrix[1]*hmatrix[3];
  dh[2]=jmatrix[2]*hmatrix[0]+jmatrix[3]*hmatrix[2];
  dh[3]=jmatrix[2]*hmatrix[1]+jmatrix[3]*hmatrix[3];

  hmatrix[0]+=dh[0]*dt/2;
  hmatrix[1]+=dh[1]*dt/2;
  hmatrix[2]+=dh[2]*dt/2;
  hmatrix[3]+=dh[3]*dt/2;

  for (long i=0; i<nres; i++) {
    tind=t0+dt*i;
    get_Jmatrix(result[i][0], result[i][1], tind, jmatrix);
    dh[0]=jmatrix[0]*hmatrix[0]+jmatrix[1]*hmatrix[2];
    dh[1]=jmatrix[0]*hmatrix[1]+jmatrix[1]*hmatrix[3];
    dh[2]=jmatrix[2]*hmatrix[0]+jmatrix[3]*hmatrix[2];
    dh[3]=jmatrix[2]*hmatrix[1]+jmatrix[3]*hmatrix[3];

    hmatrix[0]+=dh[0]*dt;
    hmatrix[1]+=dh[1]*dt;
    hmatrix[2]+=dh[2]*dt;
    hmatrix[3]+=dh[3]*dt;

  }

  //take a final Euler step:
  get_Jmatrix(result[nres][0], result[nres][1], tind, jmatrix);
  dh[0]=jmatrix[0]*hmatrix[0]+jmatrix[1]*hmatrix[2];
  dh[1]=jmatrix[0]*hmatrix[1]+jmatrix[1]*hmatrix[3];
  dh[2]=jmatrix[2]*hmatrix[0]+jmatrix[3]*hmatrix[2];
  dh[3]=jmatrix[2]*hmatrix[1]+jmatrix[3]*hmatrix[3];

  hmatrix[0]+=dh[0]*dt/2;
  hmatrix[1]+=dh[1]*dt/2;
  hmatrix[2]+=dh[2]*dt/2;
  hmatrix[3]+=dh[3]*dt/2;

}

void traj_int_obj::integrate_Hcomp_matrix(double *hmatrix) {
  float jmatrix[4];
  float dh[4];
  interpol_index tind;

  //take an Euler step:
  get_Jmatrix((result[0][0]+result[1][0])/2, 
		(result[0][1]+result[1][1])/2, tind, jmatrix);
  dh[0]=hmatrix[0]*jmatrix[0]+hmatrix[1]*jmatrix[2];
  dh[1]=hmatrix[0]*jmatrix[1]+hmatrix[1]*jmatrix[3];
  dh[2]=hmatrix[2]*jmatrix[0]+hmatrix[3]*jmatrix[2];
  dh[3]=hmatrix[2]*jmatrix[1]+hmatrix[3]*jmatrix[3];

  hmatrix[0]-=dh[0]*dt/2;
  hmatrix[1]-=dh[1]*dt/2;
  hmatrix[2]-=dh[2]*dt/2;
  hmatrix[3]-=dh[3]*dt/2;

  for (long i=0; i<nres; i++) {
    tind=t0+dt*i;
    get_Jmatrix(result[i][0], result[i][1], tind, jmatrix);
    dh[0]=hmatrix[0]*jmatrix[0]+hmatrix[1]*jmatrix[2];
    dh[1]=hmatrix[0]*jmatrix[1]+hmatrix[1]*jmatrix[3];
    dh[2]=hmatrix[2]*jmatrix[0]+hmatrix[3]*jmatrix[2];
    dh[3]=hmatrix[2]*jmatrix[1]+hmatrix[3]*jmatrix[3];

    hmatrix[0]-=dh[0]*dt;
    hmatrix[1]-=dh[1]*dt;
    hmatrix[2]-=dh[2]*dt;
    hmatrix[3]-=dh[3]*dt;

  }

  //take a final Euler step:
  get_Jmatrix(result[nres][0], result[nres][1], tind, jmatrix);
  dh[0]=hmatrix[0]*jmatrix[0]+hmatrix[1]*jmatrix[2];
  dh[1]=hmatrix[0]*jmatrix[1]+hmatrix[1]*jmatrix[3];
  dh[2]=hmatrix[2]*jmatrix[0]+hmatrix[3]*jmatrix[2];
  dh[3]=hmatrix[2]*jmatrix[1]+hmatrix[3]*jmatrix[3];

  hmatrix[0]-=dh[0]*dt/2;
  hmatrix[1]-=dh[1]*dt/2;
  hmatrix[2]-=dh[2]*dt/2;
  hmatrix[3]-=dh[3]*dt/2;

}


long traj_int_obj::print_result(FILE *fs) {
  time_class tcur;
  char tstr[MAX_TSTR_WIDTH];

  for (long i=0; i<=nres; i++) {
    tcur=get_time(t0+dt*i);
    tcur.write_string(tstr);
    fprintf(fs, "%s %10.6g %10.6g\n", tstr, result[i][0], result[i][1]);
  }

  return nres;
}

long traj_int_obj::print_result(FILE *fs, short hemi) {
  time_class tcur;
  char tstr[MAX_TSTR_WIDTH];
  float lon, lat;

  for (long i=1; i<=nres; i++) {
    tcur=get_time(t0+dt*i);
    tcur.write_string(tstr);
    tcoord2_2lonlat(result[i][0], result[i][1], 1, hemi, lon, lat);
    fprintf(fs, "%s %10.6g %10.6g\n", tstr, lon, lat);
  }

  return nres;

}

traj_int_obj::~traj_int_obj() {

  for (long i=0; i<=nres; i++) delete [] result[i];
  delete result;

//  delete u;
//  delete v;
  
  fclose(vfield_swap);

//  delete xgrid;
//  delete ygrid;
//  delete tgrid;

}
