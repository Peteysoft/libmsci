#include <stdlib.h>

#include "ctraj_tfield_nd.h"
#include "parse_command_opts.h"
#include "peteys_tmpl_lib.h"
#include "tcoord_defs.h"

using namespace libpetey;
using namespace datasets;

namespace ctraj {

template <class real>
ctraj_tfield_nd<real>::ctraj_tfield_nd(int32_t nd) {
  grid=NULL;
  tracer=NULL;
  n_dim=nd;
  map=NULL;
  inverse_map=NULL;
  nmap=0;
}

template <class real>
ctraj_tfield_nd<real>::~ctraj_tfield_nd() {
  if (tracer!=NULL) delete tracer;
  if (map!=NULL) delete map;
  if (grid!=NULL) {
    for (int32_t i=0; i<n_dim; i++) if (grid[i]!=NULL) delete grid[i];
    delete [] grid;
  }
  if (inverse_map!=NULL) delete [] inverse_map;
}

template <class real>
void ctraj_tfield_nd<real>::help(FILE *fs) {
  fprintf(fs, "Tracer model: n-dimensional Cartesian\n");
  fprintf(fs, "    Specify number of grid points and/or sidelengths as (e.g.):\n");
  fprintf(fs, "    -n n1xn2xn3...\n");
  fprintf(fs, "    -r option also specifies tracer domain radius\n");
  fprintf(fs, "    and overrides -q option\n");
  ctraj_optargs(fs, "rnxyXYq", 1);
}


//-r sidelength/2 + domain radius
//-n n1 x n2 x n3 ...
//-x nx
//-y ny
//-X dx
//-Y dy
//-q d1 x d2 x d3 ...

template <class real>
int ctraj_tfield_nd<real>::setup(int argc, char **argv) {
  void *optarg[10];
  int flag[10];
  char *str;

  float *sidelength=NULL;	//sidelength (all sides)
  float sl;
  int32_t ns;		//number of sidelengths parsed
  int32_t *np=NULL;		//number of gridpoints
  int32_t n;
  int32_t nnp;		//number of gridpoint numbers parsed
  int nleft;

  int32_t nx;		//number of x grids
  int32_t ny;		//number of y grids
  real dx;		//x sidelength
  real dy;		//y sidelength

  //set defaults:
  sl=SIDELENGTH_Q;
  n=NGRID_Q;
  real radius=0;

  optarg[0]=&radius;
  optarg[2]=&nx;
  optarg[3]=&ny;
  optarg[4]=&dx;
  optarg[5]=&dy;

  nleft=parse_command_opts(argc, argv, "rnxyXYq", "%g%s%d%d%g%g%s", optarg, flag, OPT_WHITESPACE+2);

  if (nleft<0) {
    fprintf(stderr, "tfield_nd: command option parse error\n");
  }

  //parse the number of grid points:
  if (flag[1]) {
    str=(char *) optarg[1];
    nnp=1;
    for (int i=0; str[i]!='\0'; i++) 
		if (str[i]=='x' || str[i]=='X') nnp++;
    if (n_dim!=0 && nnp!=1 && nnp!=n_dim) {
      fprintf(stderr, "Error: tfield_nd; error parsing command options\n");
      fprintf(stderr, "      number of dimensions must match\n");
      exit(FATAL_COMMAND_OPTION_PARSE_ERROR);
    }
    if (nnp==1) {
      sscanf(str, "%d", &n);
      np=NULL;
    } else {
      np=new int32_t[nnp];
      for (int32_t i=0; i<nnp; i++) {
        sscanf(str, "%d", np+i);
        for (; str[i]!='\0'; str++) {
          if (str[i]=='x' || str[i]=='X') break;
        }
        str++;
      }
    }
  } else {
    nnp=0;
  }

  //if the -r option has been set, we assume equal sidelengths all around
  if (flag[0]) {
    if (n_dim==0) if (nnp<=1) n_dim=2; else n_dim=nnp;
    sidelength=new float[n_dim];
    for (int32_t i=0; i<n_dim; i++) {
      sidelength[i]=radius*2;
    }
    
  } else if (flag[6]) {
    str=(char *) optarg[6];
    ns=0;
    for (int i=0; str[i]!='\0'; i++) 
		if (str[i]=='x' || str[i]=='X') ns++;
    if (ns!=1 && ((n_dim!=0 && ns!=n_dim) || ns!=nnp))  {
      fprintf(stderr, "Error: tfield_nd; error parsing command options\n");
      fprintf(stderr, "      number of dimensions must match\n");
      exit(FATAL_COMMAND_OPTION_PARSE_ERROR);
    }
    if (ns==1) {
      sscanf(str, "%g", &sl);
      sidelength=NULL;
    } else {
      sidelength=new float[ns];
      for (int32_t i=0; i<ns; i++) {
        sscanf(str, "%g", sidelength+i);
        for (; str[i]!='\0'; str++) {
          if (str[i]=='x' || str[i]=='X') break;
        }
        str++;
      }
    }
  }

//*FLAG* -- logic is confusing: source of bug in lla2aeb and extract_field
  //how many dimensions?
  if (n_dim==0) {
    if (nnp==1 && ns==1) n_dim=2;
    else if (nnp > ns) n_dim=nnp; else n_dim=ns;
  }

  if (np==NULL) {
    np=new int32_t[n_dim];
    for (int32_t i=0; i<n_dim; i++) np[i]=n;
  }
  if (flag[2]) np[0]=nx;
  if (flag[3]) np[1]=ny;

  if (sidelength==NULL) {
    sidelength=new float[n_dim];
    for (int32_t i=0; i<n_dim; i++) sidelength[i]=sl;
  }
  if (flag[4]) sidelength[0]=dx;
  if (flag[5]) sidelength[1]=dy;

  grid=new simple<real>*[n_dim];
  for (int32_t i=0; i<n_dim; i++) {
    grid[i]=new simple<real>(-sidelength[i]/2, sidelength[i]/2, np[i]);
  }

  tracer=new dependent_intc((simple_dataset **) grid, n_dim);

  //we want to exclude the points outside the radius:
  if (flag[0]) {
    real *testloc;
    real r;
    real dr;
    sub_1d_type nel1, ind;

    dr=0;
    for (int32_t i=0; i<n_dim; i++) {
      dr+=sidelength[i]*sidelength[i]/np[i]/np[i];
    }
    dr=sqrt(dr);

    nmap=0;
    nel1=nel();
    map=new dependent<sub_1d_type>((simple_dataset **) grid, n_dim);
    testloc=new real[n_dim];
    for (sub_1d_type i=0; i<nel1; i++) {
      get_loc_raw(i, testloc);
      r=0;
      for (int32_t j=0; j<n_dim; j++) r+=testloc[j]*testloc[j];
      r=sqrt(r);
      if (r <= radius+dr) {
        map->cel_1d(nmap, i);
        nmap++;
      } else {
        map->cel_1d(-1, i);
      }
    }
    inverse_map=new sub_1d_type[nel1];
    for (sub_1d_type i=0; i<nel1; i++) {
      map->get_1d(ind, i);
      if (ind >=0) inverse_map[ind]=i;
    }

    delete [] testloc;
  }

  //clean up:
  delete [] np;
  delete [] sidelength;

  return nleft;

}

template <class real>
void ctraj_tfield_nd<real>::get_loc_raw(sub_1d_type ind, real *loc) {
  ind_type indices[n_dim];
  tracer->indices(ind, indices);
  for (int32_t i=0; i<n_dim; i++) grid[i]->get(loc[i], indices[i]);
}

template <class real>
int32_t ctraj_tfield_nd<real>::get_loc(int32_t ind, real *loc) {
  sub_1d_type ind1;

  if (nmap>0) ind1=inverse_map[ind]; else ind1=ind;
  get_loc_raw(ind1, loc);

  return 0;

}

template <class real>
int32_t ctraj_tfield_nd<real>::interpolate(int32_t domain, real *loc, int32_t *ind, double *wt) {
  double intind[n_dim];
  int32_t d1=nwt();

  for (int32_t i=0; i<n_dim; i++) intind[i]=grid[i]->interp(loc[i]);
  tracer->interpol_coeff(intind, ind, wt);

  if (nmap>0) {
    sub_1d_type ind0;
    int32_t nvert=d1;

    for (int32_t i=0; i<nvert; i++) {
      map->get_1d(ind0, ind[i]);
      ind[i]=ind0;
      if (ind0<0) d1=-1;
    }
  } else {
    real min, max;
    for (int32_t i=0; i<n_dim; i++) {
      grid[i]->get(min, 0);
      grid[i]->get(max, grid[i]->nel()-1);
      if (loc[i] < min || loc[i] > max) {
        d1=-1;
        break;
      }
    }
  }

  return d1;
}

template <class real>
real *ctraj_tfield_nd<real>::to(real *input, real *parm) {
  real *qvec;
  int32_t nin;
  int32_t nout;
  sub_1d_type ind;

  nout=nel();
  qvec=new real[nout];

  if (nmap > 0) {
    nin=1;
    for (int32_t i=0; i<n_dim; i++) nin=grid[i]->nel()*nin;

    for (int i=0; i<nin; i++) {
      map->get_1d(ind, i);
      if (ind>=0) qvec[ind]=input[i];
    }
  } else {
    for (int i=0; i<nout; i++) qvec[i]=input[i];
  }

  return qvec;

}

template <class real>
real *ctraj_tfield_nd<real>::from(real *input, real *parm) {
  real *qout;
  int32_t nin;
  int32_t nout;
  sub_1d_type ind;
  real missing=parm[0];

  nin=nel();

  if (nmap > 0) {
    nout=1;
    for (int32_t i=0; i<n_dim; i++) nout=grid[i]->nel()*nout;
    qout=new real[nout];
    for (int i=0; i<nout; i++) {
      map->get_1d(ind, i);
      if (ind<0) {
        qout[i]=missing;
      } else {
        qout[i]=input[ind];
      }
      //qout[inverse_map[i]]=input[i];
    }
    parm[1]=nout;
  } else {
    qout=new real[nin];
    for (int i=0; i<nin; i++) qout[i]=input[i];
  }

  return qout;

}

template <class real>
int ctraj_tfield_nd<real>::get_raw_grids(int32_t *grids) {
  for (int i=0; i<n_dim; i++) {
    grids[i]=grid[i]->nel();
  }
  return 0;
}

template <class real>
int ctraj_tfield_nd<real>::get_range(real *low, real *high) {
  for (int i=0; i<n_dim; i++) {
    grid[i]->get(low[i], 0);
    grid[i]->get(high[i], grid[i]->nel()-1);
  }
  return 0;
}

template <class real>
int32_t ctraj_tfield_nd<real>::nel() {
  real nel;

  if (nmap > 0) return nmap;
  nel=1;
  for (int32_t i=0; i<n_dim; i++) nel=grid[i]->nel()*nel;
  return nel;
}

template <class real>
int32_t ctraj_tfield_nd<real>::ndim() {
  return n_dim;
}

template <class real>
int32_t ctraj_tfield_nd<real>::nwt() {
  int32_t nvert;
  nvert=1;
  for (int32_t i=0; i<n_dim; i++) nvert*=2;
  return nvert;
}

template class ctraj_tfield_nd<float>;

} //end namespace ctraj

