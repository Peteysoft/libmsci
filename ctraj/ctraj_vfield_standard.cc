#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "parse_command_opts.h"
#include "error_codes.h"

#include "ctraj_defaults.h"
#include "coordtran.h"
#include "ctraj_vfield_standard.h"

using namespace libpetey;
using namespace datasets;

namespace ctraj {

template <class real>
int ctraj_vfield_standard<real>::setup(int argc, char **argv) {
  void *optarg[10];
  int flag[10];
  int ref=1;			//reference domain
  time_class t1, t2;
  int ret;
  int err=1;
  int err2;
  char *fname[2];
  int64_t page_size=-1;		//vfield page size in bytes

  optarg[2]=&page_size;

  ret=parse_command_opts(argc, argv, "-+B", "%%%ld", optarg, flag, OPT_WHITESPACE+2);
  if (ret < 0) {
    fprintf(stderr, "vfield_standard: error parsing command line");
    err=-1;
    ret=ret*err;
  }
  if (ret < 3) {
    fprintf(stderr, "vfield_standard: insufficient command line arguments\n");
    return -ret;
  }
  if (flag[0]) ref=0;
  if (flag[1]) ref=1;

  fname[0]=argv[1];
  fname[1]=argv[2];

  err2=init(fname, page_size, "r", ref);

  if (err2!=0) err=-1;

  for (int i=2; i<ret; i++) {
    argv[i-2]=argv[i];
  }
  argv[ret-1]=fname[1];
  argv[ret-2]=fname[0];

  return (ret-2)*err;
}

template <class real>
void ctraj_vfield_standard<real>::help(FILE *fs) {
  fprintf(fs, "Velocity field: 2-D global on azimuthal equidistant coords\n");
  fprintf(fs, "arg1 Southern hemisphere velocity field\n");
  fprintf(fs, "arg2 Northern hemisphere velocity field\n");
  ctraj_optargs(fs, "-+B", 1);
}

template <class real>
int ctraj_vfield_standard<real>::init(char **fname, int64_t page_size, const char *mode, int32_t ref) {
  simple<real> *z1[2];
  simple<real> *h[2];
  real h1, h2;
  long loc, dum;
  time_class t0, t1;
  double dt0, dt1;

  refd=ref;

  for (int i=0; i<2; i++) {
    fprintf(stderr, "Opening: %s\n", fname[i]);
    fs[i]=fopen(fname[i], mode);
    if (fs[i]==NULL) {
      fprintf(stderr, "vfield_standard: could not open, %s, for reading\n", fname[i]);
      return UNABLE_TO_OPEN_FILE_FOR_READING;
    }
  }

  if (strcmp(mode, "r")==0 || strcmp(mode, "r+")==0) {
    for (int i=0; i<2; i++) {
      all[i]->read(fs[i]);
      loc=all[i]->search_var("xgrid", dum);
      if (loc<0) {
        fprintf(stderr, "vfield_standard: variable, %s, not found in, %s\n", "xgrid", fname[i]);
        return FILE_READ_ERROR;
      }
      x[i]=(simple<real> *) all[i]->get_var(loc);
      loc=all[i]->search_var("ygrid", dum);
      if (loc<0) {
        fprintf(stderr, "vfield_standard: variable, %s, not found in, %s\n", "ygrid", fname[i]);
        return FILE_READ_ERROR;
      }
      y[i]=(simple<real> *) all[i]->get_var(loc);
      loc=all[i]->search_var("zgrid", dum);
      if (loc<0) {
        fprintf(stderr, "vfield_standard: variable, %s, not found in, %s\n", "zgrid", fname[i]);
        return FILE_READ_ERROR;
      }
      z1[i]=(simple<real> *) all[i]->get_var(loc);
      loc=all[i]->search_var("height", dum);
      if (loc<0) {
        h[i]=NULL;
      } else {
        h[i]=(simple<real> *) all[i]->get_var(loc);
      }
      loc=all[i]->search_var("time", dum);
      if (loc<0) {
        fprintf(stderr, "vfield_standard: variable, %s, not found in, %s\n", "time", fname[i]);
        return FILE_READ_ERROR;
      }
      t[i]=(simple<time_class> *) all[i]->get_var(loc);
      loc=all[i]->search_var("u", dum);
      if (loc<0) {
        fprintf(stderr, "vfield_standard: variable, %s, not found in, %s\n", "u", fname[i]);
        return FILE_READ_ERROR;
      }
      U[i]=(dependent_swap<real> *) all[i]->get_var(loc);
      loc=all[i]->search_var("v", dum);
      if (loc<0) {
        fprintf(stderr, "vfield_standard: variable, %s, not found in, %s\n", "v", fname[i]);
        return FILE_READ_ERROR;
      }
      V[i]=(dependent_swap<real> *) all[i]->get_var(loc);
  
      //set page size:
      if (page_size > 0) {
        U[i]->set_page_size(page_size, 1);
        V[i]->set_page_size(page_size, 1);
      }

    }
    if (!(*z1[0]==*z1[1])) {
      fprintf(stderr, "vfield_standard: velocity fields must be on the same vertical level\n");
      return FILE_READ_ERROR;
    }
    z=z1[0];
 
    if (h[0] != NULL && h[1] != NULL) {
      h[0]->get(h1, 0);
      h[1]->get(h2, 0);
      metric=new az_eq_t<real>(REARTH+(h1+h2)/2);
    } else {
      metric=new az_eq_t<real>(REARTH);
    }

    //synchronize the two time grids:
    t[0]->get(t0, 0);
    t[0]->get(t1, 1);
    dt0=(double) (t1-t0);
    t[1]->get(t0, 0);
    t[1]->get(t1, 1);
    dt1=(double) (t1-t0);
    t[0]->get(t0, 0);
    t[1]->get(t1, 0);
    if (refd) {
      tmult=dt0/dt1;
      if (t1>t0) {
        tdiff=t[0]->interp(t1);
      } else {
        tdiff=-t[1]->interp(t0)*tmult;
      }
    } else {
      tmult=dt1/dt0;
      if (t0>t1) {
        tdiff=t[1]->interp(t0);
      } else {
        tdiff=-t[0]->interp(t1)*tmult;
      }
    }
    //printf("tmult=%g; tdiff=%g\n", tmult, tdiff);
    //printf("time grids:\n");
    //t[0]->print();
    //t[1]->print();

    //search in the "reference" dataset for input field grids
    //and interpolation coefficients:
    if (strcmp(mode, "r+")==0) {
      simple<real> *loncheck=NULL;
      simple<real> *latcheck=NULL;

      lon=NULL;
      lat=NULL;
      c1[0]=NULL;
      c2[0]=NULL;
      c1[1]=NULL;
      c2[1]=NULL;

      loc=all[refd]->search_var("lon", dum);
      if (loc>=0) {
        lon=(simple<real> *) all[refd]->get_var(loc);
        loc=all[refd]->search_var("lat", dum);
        lat=(simple<real> *) all[refd]->get_var(loc);
        fprintf(stdin, "ctraj_field_standard: found input grid\n");
        fprintf(stdin, "lon:");
        lon->print();
        fprintf(stdin, "lat:");
        lat->print();
      }

      if (refd==1) ref=0; else ref=1;
      loc=all[ref]->search_var("lon", dum);
      if (loc>=0) {
        loncheck=(simple<real> *) all[ref]->get_var(loc);
        loc=all[ref]->search_var("lat", dum);
        latcheck=(simple<real> *) all[ref]->get_var(loc);
      }

      if (lon!=NULL) {
        loc=all[refd]->search_var("c1", dum);
        if (loc >= 0) {
          c1[refd]=(dependent_swap<interpol_index> *) all[refd]->get_var(loc);
          loc=all[refd]->search_var("c2", dum);
          c2[refd]=(dependent_swap<interpol_index> *) all[refd]->get_var(loc);
        } else {
          c1[refd]=new dependent<interpol_index>(x[refd], y[refd]);
          c2[refd]=new dependent<interpol_index>(x[refd], y[refd]);
          //calculate the interpolation coefficients:
          intcoeff(lon, lat, x[refd], y[refd], c1[refd], c2[refd], 2*refd-1);
        }

        if (refd==1) ref=0; else ref=1;
        if (*lon==*loncheck && *lat==*latcheck) {
          loc=all[ref]->search_var("c1", dum);
          if (loc > 0) {
            c1[ref]=(dependent_swap<interpol_index> *) all[ref]->get_var(loc);
            loc=all[ref]->search_var("c2", dum);
            c2[ref]=(dependent_swap<interpol_index> *) all[ref]->get_var(loc);
          }
        }
        if (c1[ref]==NULL) {
          c1[ref]=new dependent<interpol_index>(x[ref], y[ref]);
          c2[ref]=new dependent<interpol_index>(x[ref], y[ref]);
          //calculate the interpolation coefficients:
          intcoeff(lon, lat, x[ref], y[ref], c1[ref], c2[ref], 2*ref-1);
        }
      }
    }
  } else {
    for (int i=0; i<2; i++) {
      t[i]=NULL;
      x[i]=NULL;
      y[i]=NULL;
      U[i]=NULL;
      V[i]=NULL;
      c1[i]=NULL;
      c2[i]=NULL;
    }
    z=NULL;
    lon=NULL;
    lat=NULL;
    metric=new az_eq_t<real>(REARTH);
  }
  return 0;
  
}
  
template <class real>
ctraj_vfield_standard<real>::ctraj_vfield_standard() {
  all[0]=new composite_dataset();
  all[1]=new composite_dataset();
  fs[0]=NULL;
  fs[1]=NULL;
  metric=NULL;
}

template <class real>
ctraj_vfield_standard<real>::~ctraj_vfield_standard() {
  delete all[0];
  delete all[1];
  if (fs[0]!=NULL) fclose(fs[0]);
  if (fs[1]!=NULL) fclose(fs[1]);
  delete metric;
}

template <class real>
int ctraj_vfield_standard<real>::v(int32_t domain, double tind, real *x1, real *v) {
  double xind, yind;

  if (domain!=refd) tind=tind*tmult+tdiff;

  xind=x[domain]->interp(x1[0]);
  yind=y[domain]->interp(x1[1]);

  U[domain]->interpol(v[0], xind, yind, 0, tind);
  V[domain]->interpol(v[1], xind, yind, 0, tind);

}

template <class real>
int32_t ctraj_vfield_standard<real>::fix(int32_t domain, double tind, real *x) {
  int hemi;
  int32_t domain_new;
  time_class date;

  hemi=2*domain-1;
  //tcoord_fix(x[0], x[1], hemi);
  hemi=metric->fix(hemi, x);
  domain_new=(hemi+1)/2;

  return domain_new;
}

template <class real>
int32_t ctraj_vfield_standard<real>::absolute(int32_t domain, real *x) {
  int hemi;
  real x1[2];

  //tcoord2_2lonlat(x1[0], x1[1], dir, hemi, x2[0], x2[1]);
  if (domain < 0) {
    hemi=metric->to(x, x1);
    x[0]=x1[0];
    x[1]=x1[1];
    domain=(hemi+1)/2;
  } else {
    hemi=2*domain-1;
    metric->from(hemi, x, x1);
    x[0]=x1[0];
    x[1]=x1[1];
    domain=-1;
  }

  return domain;

}

template <class real>
int32_t ctraj_vfield_standard<real>::reference(int32_t domain, real *x) {
  int hemi;

  if (domain<0) {
    //if (domain!=refd) tcoord_N2S(x2[0], x2[1]);
    hemi=2*refd-1;
    hemi=metric->fix(hemi, x);
    //tcoord_fix(x2[0], x2[1], hemi);
    domain=(hemi+1)/2;
  } else {
    if (domain!=refd) metric->swap(x);
    domain=refd;
  }

  return domain;
}

template <class real>
int32_t ctraj_vfield_standard<real>::ndim() {
  return 2;
}

template <class real>
double ctraj_vfield_standard<real>::maxt() {
  time_class t1, t2;
  t[0]->get(t1, t[0]->nel()-1);
  t[1]->get(t2, t[1]->nel()-1);
  if (t1 > t2) t1=t2;
  return t[refd]->interp(t1);
}

template <class real>
double ctraj_vfield_standard<real>::get_tind(time_class date) {
  interpol_index tind;
  return t[refd]->interp(date);
}

template <class real>
time_class ctraj_vfield_standard<real>::get_t(double tind) {
  time_class t0;
  t[refd]->get(t0, tind);
  return t0;
}

template <class real>
double ctraj_vfield_standard<real>::get_tind(char *date) {
  time_class t0;
  interpol_index tind;
  t0.read_string(date);
  return t[refd]->interp(t0);
}

template <class real>
int ctraj_vfield_standard<real>::jmat(int32_t d, double tind, real *loc, real *j) {
  //interpolation indices for x and y:
  interpol_index xind, yind;
  //winds at all four corners:
  real u1, u2, u11, u12, u21, u22;
  real v1, v2, v11, v12, v21, v22;
  //x and y coords at all four corners:
  real x1, x2;
  real y1, y2;
  //average values of x and y:
  //change in x and y:
  real dx, dy;
  //integer indices:
  ind_type lx, ly, lt;
  interpol_index frac;

  //metric coefficients:
  real mcoef2[2];
  real cx, cy;

  //for the time interval:
  time_class t1, t2;
  double delt;

  xind=x[d]->interp(loc[0]);
  yind=y[d]->interp(loc[1]);

  lx=(ind_type) xind;
  ly=(ind_type) yind;
  lt=(ind_type) tind;

  x[d]->get(x1, lx);
  x[d]->get(x2, lx+1);
  y[d]->get(y1, ly);
  y[d]->get(y2, ly+1);

  //calculate metric coefficients:
  metric->mcoef2(loc, mcoef2);
  cx=sqrt(mcoef2[0]);
  cy=sqrt(mcoef2[1]);

  dx=(x2-x1)*cx;
  dy=(y2-y1)*cy;

  //get the winds at all four corners:
  frac=tind-(interpol_index) lt;
  U[d]->get(u1, lx, ly, 0, lt);
  U[d]->get(u2, lx, ly, 0, lt+1);
  //printf("u1=%f, u2=%f\n", u1, u2);
  u11=u1*(1-frac)+u2*frac;
  U[d]->get(u1, lx, ly+1, 0, lt);
  U[d]->get(u2, lx, ly+1, 0, lt+1);
  //printf("u1=%f, u2=%f\n", u1, u2);
  u12=u1*(1-frac)+u2*frac;
  U[d]->get(u1, lx+1, ly, 0, lt);
  U[d]->get(u2, lx+1, ly, 0, lt+1);
  u21=u1*(1-frac)+u2*frac;
  U[d]->get(u1, lx+1, ly+1, 0, lt);
  U[d]->get(u2, lx+1, ly+1, 0, lt+1);
  u22=u1*(1-frac)+u2*frac;
  //printf("u11=%f, u12=%f, u21=%f, u22=%f\n", u11, u12, u12, u22);

  V[d]->get(v1, lx, ly, 0, lt);
  V[d]->get(v2, lx, ly, 0, lt+1);
  v11=v1*(1-frac)+v2*frac;
  V[d]->get(v1, lx, ly+1, 0, lt);
  V[d]->get(v2, lx, ly+1, 0, lt+1);
  v12=v1*(1-frac)+v2*frac;
  V[d]->get(v1, lx+1, ly, 0, lt);
  V[d]->get(v2, lx+1, ly, 0, lt+1);
  v21=v1*(1-frac)+v2*frac;
  V[d]->get(v1, lx+1, ly+1, 0, lt);
  V[d]->get(v2, lx+1, ly+1, 0, lt+1);
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

template <class real>
int ctraj_vfield_standard<real>::get_t(double tind, char *date) {
  time_class t0;
  t[refd]->get(t0, tind);
  t0.write_string(date);

  return date!=NULL;
}

template <class real>
az_eq_t<real> * ctraj_vfield_standard<real>::get_metric() {
  return metric;
}

template <class real>
int ctraj_vfield_standard<real>::set_outgrid(real zlev, int32_t ngrid, real sl2) {
  long loc, dum;
  simple<real> *z1;

  if (x[0]!=NULL) {
    fprintf(stderr, "ctraj_vfield_standard: Output grids already set\n");
    return -1;
  }

  z=new simple<real>(&zlev, 1);
  z1=new simple<real>(&zlev, 1);
  loc=all[0]->add_var("zgrid");
  all[0]->cvar(loc, (dataset *) z);
  loc=all[1]->add_var("zgrid");
  all[1]->cvar(loc, (dataset *) z1);

  for (int i=0; i<2; i++) {
    x[i]=new simple<real>(-sl2, sl2, ngrid);
    y[i]=new simple<real>(-sl2, sl2, ngrid);

    loc=all[i]->add_var("xgrid");
    all[i]->cvar(loc, (dataset *) x[i]);
    loc=all[i]->add_var("ygrid");
    all[i]->cvar(loc, (dataset *) y[i]);

  }

  return 0;
}

template <class real>
int ctraj_vfield_standard<real>::set_tgrid(const simple<time_class> &t1) {
  long loc, dum;

  if (t[0]!=NULL) {
    fprintf(stderr, "ctraj_vfield_standard: time grid already defined\n");
    return -1;
  }

  for (int i=0; i<2; i++) {
    t[i]=new simple<time_class>(t1);
    loc=all[i]->add_var("time");
    all[i]->cvar(loc, (dataset *) t[i]);
  }
  tdiff=0;
  tmult=1;

  return 0;
}

template <class real>
int ctraj_vfield_standard<real>::set_ingrid(const simple<real> &lonin, const simple<real> &latin) {
  long loc, dum;

  if (lon!=NULL) {
    delete lon;
    delete lat;
  }
  if (c1[0]!=NULL) {
    delete c1[0];
    delete c2[0];
  }
  if (c1[1]!=NULL) {
    delete c1[1];
    delete c2[1];
  }

  //bit stupid to copy all this shit from one object to another
  //but then again these "datasets" are pretty stupid anyways
  //at least the way they are currently implemented...
  lon=new simple<real>(lonin);
  lat=new simple<real>(latin);

  for (int i=0; i<2; i++) {
    c1[i]=new dependent<interpol_index>(x[i], y[i]);
    c2[i]=new dependent<interpol_index>(x[i], y[i]);
    //calculate the interpolation coefficients:
    intcoeff(lon, lat, x[i], y[i], c1[i], c2[i], 2*i-1);
  }

}

template <class real>
int ctraj_vfield_standard<real>::write(int keep_input_grids) {
  long loc, dum; 
  int32_t unref;
  simple<real> *lon1, *lat1;
  simple<real> *z1;

  if (keep_input_grids) {
    lon1=new simple<real>(*lon);
    lat1=new simple<real>(*lat);

    loc=all[refd]->search_var("lon", dum);
    if (loc<0) loc=all[refd]->add_var("lon");
    all[refd]->cvar(loc, (dataset *) lon);
    loc=all[refd]->search_var("lat", dum);
    if (loc<0) loc=all[refd]->add_var("lat");
    all[refd]->cvar(loc, (dataset *) lat);

    if (refd==0) unref=1; else unref=0;
    loc=all[unref]->search_var("lon", dum);
    if (loc<0) loc=all[unref]->add_var("lon");
    all[unref]->cvar(loc, (dataset *) lon1);
    loc=all[unref]->search_var("lat", dum);
    if (loc<0) loc=all[unref]->add_var("lat");
    all[unref]->cvar(loc, (dataset *) lat1);

    for (int i=0; i<2; i++) {
      loc=all[i]->search_var("c1", dum);
      if (loc<0) loc=all[i]->add_var("c1");
      all[i]->cvar(loc, (dataset *) c1[i]);
      loc=all[i]->search_var("c2", dum);
      if (loc<0) loc=all[i]->add_var("c2");
      all[i]->cvar(loc, (dataset *) c2[i]);
    }
  }

  for (int i=0; i<2; i++) {
    loc=all[i]->search_var("zgrid", dum);
    z1=(simple<real> *) all[i]->get_var(loc);
    U[i]=new dependent_swap<real>(x[i], y[i], z1, t[i]);
    V[i]=new dependent_swap<real>(x[i], y[i], z1, t[i]);
    loc=all[i]->add_var("u");
    all[i]->cvar(loc, (dataset *) U[i]);
    loc=all[i]->add_var("v");
    all[i]->cvar(loc, (dataset *) V[i]);
  }

  all[0]->write(fs[0]);
  all[1]->write(fs[1]);

}

template <class real>
dependent<real> * ctraj_vfield_standard<real>::empty_field() {
  return new dependent<real>(lon, lat);
}

template <class real>
int ctraj_vfield_standard<real>::add_field(int32_t tind, dependent<real> *uu, dependent<real> *vv) {
  double c1val, c2val;
  real vval[2], vtval[2];
  real xval[2];
  time_class ttmpa, ttmpb;
  int32_t nt;
  double tint;
  real vunit;
  int32_t nx, ny;
  double tindf;
  int32_t tind1;

  for (int dom=0; dom<2; dom++) {
    //calculate time interval which is used to normalize the velocities
    nt=t[dom]->nel();
    t[dom]->get(ttmpa, tind);
    if (tind == 0) {
      t[dom]->get(ttmpb, (ind_type) 1);
      tint=ttmpb.diff(ttmpa);
    } else if (tind == nt-1) {
      t[dom]->get(ttmpb, nt-2);
      tint=ttmpa.diff(ttmpb);
    } else {
      t[dom]->get(ttmpa, tind-1);
      t[dom]->get(ttmpb, tind+1);
      tint=ttmpb.diff(ttmpa)/2;
    }
    vunit=(real) tint*SECPERDAY/MPERKM;//convert from m/s to km/(delta t=6 h
    //printf("vunit=%g\n", vunit);

    nx=x[dom]->nel();
    ny=y[dom]->nel();

    for (ind_type j=0; j<ny; j++) {
      for (ind_type i=0; i<nx; i++) {
        //do the interpolation:
        c1[dom]->get(c1val, i, j);
        c2[dom]->get(c2val, i, j);
        uu->interpol(vval[0], c1val, c2val);
        vv->interpol(vval[1], c1val, c2val);
        //transform the velocity:
        x[dom]->get(xval[0], i);
        y[dom]->get(xval[1], j);
        metric->vtran(2*dom-1, xval, vval, vtval);
        vtval[0]*=vunit;
        vtval[1]*=vunit;

        if (i==refd) {
          tind1=tind;
        } else {
          tindf=tmult*tind+tdiff;
          tind1=tindf;
          if (tind1!=tindf) {
            fprintf(stderr, "ctraj_vfield_standard: N & S hemisphere time grids do not line up\n");
          }
        }

        U[dom]->cel(vtval[0], i, j, 0, tind1);
        V[dom]->cel(vtval[1], i, j, 0, tind1);
        //printf("(%g, %g)\n", vtval[0], vtval[1]);
      }
    }
    //printf("\n");
  }

  return 0;

}

template class ctraj_vfield_standard<float>;

} //end namespace ctraj

