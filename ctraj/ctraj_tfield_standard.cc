#include "parse_command_opts.h"
#include "peteys_tmpl_lib.h"
#include "tcoord_defs.h"

#include "ctraj_defaults.h"
#include "coordtran.h"
#include "ctraj_tfield_standard.h"

using namespace libpetey;
using namespace datasets;

namespace ctraj {

template <class real>
ctraj_tfield_standard<real>::ctraj_tfield_standard() {
  xgrid=NULL;
  ygrid=NULL;
  map_map=NULL;
  tracer=NULL;
  inverse_map=NULL;
}

template <class real>
ctraj_tfield_standard<real>::~ctraj_tfield_standard() {
  if (map_map!=NULL) delete map_map;
  if (tracer!=NULL) delete tracer;
  if (inverse_map!=NULL) delete [] inverse_map;
  if (xgrid!=NULL) delete xgrid;
  if (ygrid!=NULL) delete ygrid;
}

template <class real>
void ctraj_tfield_standard<real>::help(FILE *fs) {
  fprintf(fs, "Tracer model: 2-D global on azimuthal-equidistant coords\n");
  ctraj_optargs(fs, "rn", OPT_WHITESPACE+2);

}

template <class real>
int ctraj_tfield_standard<real>::init1(int32_t np, real sdl2) {
  sub_1d_type ind;

  metric=new az_eq_t<real>(REARTH);

  xgrid=new simple<real>(-sdl2, sdl2, np);
  ygrid=new simple<real>(-sdl2, sdl2, np);
  //x_qgrid->print();

  tracer=new dependent_intc(xgrid, ygrid);

  map_map=new dependent<sub_1d_type>(xgrid, ygrid);

  nocorner_map(xgrid, ygrid, map_map, nmap, metric);

  //create the inverse mapping:
  inverse_map=new sub_1d_type[nmap];
  for (sub_1d_type i=0; i<np*np; i++) {
    map_map->get_1d(ind, i);
    if (ind >= 0) inverse_map[ind]=i;
  }
  return 0;
}

template <class real>
int ctraj_tfield_standard<real>::init2(int32_t n, real sld2) {
  int32_t np;
  np=2*(int32_t) sqrt(n/M_PI/2);

  return init1(np, sld2);

}

template <class real>
int ctraj_tfield_standard<real>::setup(int argc, char **argv) {
  void *optarg[10];
  int flag[10];

  real sidelength=SIDELENGTH_Q;
  int32_t np=NGRID_Q;
  int nleft;
  sub_1d_type ind;

  optarg[0]=&sidelength;
  optarg[1]=&np;

  nleft=parse_command_opts(argc, argv, "rn", "%g%d", optarg, flag, OPT_WHITESPACE+2);

  if (nleft<0) {
    fprintf(stderr, "tracer_standard: command option parse error\n");
  }

  init1(np, sidelength);

  return nleft;

}

template <class real>
void ctraj_tfield_standard<real>::set_metric(az_eq_t<real> *m) {
  if (metric!=NULL) delete metric;
  metric=m;
}


template <class real>
int32_t ctraj_tfield_standard<real>::get_loc(int32_t ind, real *loc) {
  sub_1d_type ind1;
  int32_t xyind[2];
  int32_t domain;

  domain=ind/nmap;
  ind1=inverse_map[ind%nmap];
  tracer->indices(ind1, xyind);
  //printf("[%5d %5d]\n", xyind[0], xyind[1]);

  xgrid->get(loc[0], xyind[0]);
  ygrid->get(loc[1], xyind[1]);

  return domain;

}

template <class real>
int32_t ctraj_tfield_standard<real>::interpolate(int32_t domain, real *loc, int32_t *ind, double *wt, real dt) {
  short hemi;
  real loc1[2]={loc[0], loc[1]};
  double xyind[2];
  sub_1d_type sub[4];
  sub_1d_type ind0;
  int32_t d1;

  hemi=2*domain-1;

  //tcoord_fix(xf, yf, hemi);
  hemi=metric->fix(hemi, loc1);
  d1=(hemi+1)/2;
  xyind[0]=xgrid->interp(loc1[0]);
  xyind[1]=ygrid->interp(loc1[1]);
  tracer->interpol_coeff(xyind, sub, wt);

  for (int i=0; i<4; i++) {
    map_map->get_1d(ind0, sub[i]);
    ind[i]=ind0+nmap*d1;
  }

  return 4;
}

template <class real>
real *ctraj_tfield_standard<real>::to(real *input, real *parm) {
  real *qvec;
  simple<real> *lon;
  simple<real> *lat;
  //interpolation coefficients:
  dependent<interpol_index> *c1;
  dependent<interpol_index> *c2;
  dependent<real> *q;
  int nx, ny;
  interpol_index c1val, c2val;
  sub_1d_type mapval;
  real val;
  int nlon, nlat;

  nx=xgrid->nel();
  ny=ygrid->nel();

  nlon=parm[0];
  nlat=parm[1];

  qvec=new real[nmap*2];

  c1=new dependent<interpol_index>(xgrid, ygrid);
  //c1->print_meta();
  c2=new dependent<interpol_index>(xgrid, ygrid);
  //c2->print_meta();

  lon=new simple<real>(0., 360., nlon+1);
  lat=new simple<real>(-90., 90., nlat);
  q=new dependent<real>(lon, lat);

  for (ind_type j=0; j<nlat; j++) {
    for (ind_type i=0; i<nlon; i++) {
      q->cel(input[i+j*nlon], i, j);
    }
    //"wrap" input grids for interpolation:
    q->get(val, 0, j);
    q->cel(val, nlon, j);
  }

  //Southern Hemisphere:

  //interpolation coefficients:
  intcoeff(lon, lat, xgrid, ygrid, c1, c2, -1);

  //perform the interpolation and write the result to the output file:
  for (ind_type i=0; i<nx; i++) {
    for (ind_type j=0; j<ny; j++) {
      c1->get(c1val, i, j);
      c2->get(c2val, i, j);
      q->interpol(val, c1val, c2val);
      //printf("%d %d\n", i, j);
      //printf("%lg %lg\n", c1val, c2val);
      //printf("%f\n", val);
      //fwrite(&val, 1, sizeof(val), fs);
      map_map->get(mapval, i, j);
      //printf("%d %d %d %d\n", i, j, mapval, nmap);
      if (mapval >= 0) {
        qvec[mapval]=val;
        //printf("%f\n", val);
      }
    }
  }

  //interpolation coefficients:
  intcoeff(lon, lat, xgrid, ygrid, c1, c2, 1);

  //perform the interpolation and write the result to the output file:
  for (ind_type i=0; i<nx; i++) {
    for (ind_type j=0; j<ny; j++) {
      c1->get(c1val, i, j);
      c2->get(c2val, i, j);
      q->interpol(val, c1val, c2val);
      //printf("%d %d\n", i, j);
      //printf("%lg %lg\n", c1val, c2val);
      //printf("%f\n", val);
      //fwrite(&val, 1, sizeof(val), fs);
      map_map->get(mapval, i, j);
      //printf("%d %d %d %d\n", i, j, mapval, nmap);
      if (mapval >= 0) {
        qvec[mapval+nmap]=val;
        //printf("%f\n", val);
      }
    }
  }

  delete q;
  delete c1;
  delete c2;
  delete lon;
  delete lat;

  return qvec;

}

template <class real>
real *ctraj_tfield_standard<real>::from(real *input, real *parm) {
  real *qout;
  simple<real> *lon;
  simple<real> *lat;
  dependent<real> *q1;
  //interpolation coefficients:
  dependent<interpol_index> *c1;
  dependent<interpol_index> *c2;
  int nx, ny;
  int nlon, nlat, nlats;
  int k;
  interpol_index c1val, c2val;

  nx=xgrid->nel();
  nx=xgrid->nel();
  nlon=parm[0];
  nlat=parm[1];
  nlats=nlat/2;

  qout=new real[nlon*nlat];

  lon=new simple<real>(0., 360., nlon+1);
  lat=new simple<real>(-90., 90., nlat);

  c1=new dependent<interpol_index>(lon, lat);
  //c1->print_meta();
  c2=new dependent<interpol_index>(lon, lat);

  q1=new dependent<real>(xgrid, ygrid);

  //Southern Hemisphere:
  for (int i=0; i<nmap; i++) {
    q1->cel_1d(input[i], inverse_map[i]);
  }

  intcoeff2(lon, lat, xgrid, ygrid, c1, c2, -1);

  //perform the interpolation and write the result to the output file:
  k=0;			//lazy way...
  for (ind_type j=0; j<nlats; j++) {
    for (ind_type i=0; i<nlon; i++) {
      c1->get(c1val, i, j);
      c2->get(c2val, i, j);
      q1->interpol(qout[k], c1val, c2val);
      //printf("%d %d\n", i, j);
      //printf("%lg %lg\n", c1val, c2val);
      //printf("%f\n", val);
      k++;
    }
  }

  //Northern Hemisphere:
  for (int i=0; i<nmap; i++) {
    q1->cel_1d(input[i+nmap], inverse_map[i]);
  }

  //create the output grids and interpolation coefficients:
  intcoeff2(lon, lat, xgrid, ygrid, c1, c2, 1);

  //perform the interpolation and write the result to the output file:
  for (ind_type j=nlats; j<nlat; j++) {
    for (ind_type i=0; i<nlon; i++) {
      //printf("%d %d\n", i, j);
      c1->get(c1val, i, j);
      c2->get(c2val, i, j);
      q1->interpol(qout[k], c1val, c2val);
      k++;
    }
  }

  delete q1;
  delete c1;
  delete c2;
  delete lon;
  delete lat;

  return qout;

}

template <class real>
int ctraj_tfield_standard<real>::get_raw_grids(int32_t *grids) {
  grids[0]=xgrid->nel();
  grids[1]=ygrid->nel();
  return 0;
}

template <class real>
int ctraj_tfield_standard<real>::get_range(real *low, real *high) {
  xgrid->get(low[0], 0);
  ygrid->get(low[1], 0);
  xgrid->get(high[0], xgrid->nel()-1);
  ygrid->get(high[1], ygrid->nel()-1);
  return 0;
}


template <class real>
int32_t ctraj_tfield_standard<real>::nel() {
  return nmap*2;
}

template <class real>
int32_t ctraj_tfield_standard<real>::ndim() {
  return 2;
}

template <class real>
int32_t ctraj_tfield_standard<real>::nwt() {
  return 4;
}

template class ctraj_tfield_standard<float>;

} //end namespace ctraj

