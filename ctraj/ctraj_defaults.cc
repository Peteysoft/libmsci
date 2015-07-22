#include <math.h>

#include "parse_command_opts.h"

#include "ctraj_defaults.h"
#include "ctraj_vfield_standard.h"
#include "ctraj_vfield_2d.h"
#include "ctraj_vfield_anal.h"

using namespace libpetey;

namespace ctraj {

template <class real>
int ctraj_deriv(double t, real *x, real *v, void *param) {
  void **params;
  int32_t domain;
  int err;
  ctraj_vfield_base<real> *vfield;

  params=(void **) param;
  vfield=(ctraj_vfield_base<real> *) params[0];
  domain=* (int32_t *) params[1];

  err=vfield->v(domain, t, x, v);

  //printf("%d %g (%g, %g): (%g, %g)\n", domain, t, x[0], x[1], v[0], v[1]);

  return err;
}

template int ctraj_deriv<float>(double t, float *x, float *v, void *param);

template <class real>
int ctraj_deriv_L(real t, real *x, real *v, void *param) {
  void **params;
  int32_t domain;
  int err;
  int32_t ndim;
  real *jmat;
  int32_t k;
  ctraj_vfield_base<real> *vfield;

  params=(void **) param;
  vfield=(ctraj_vfield_base<real> *) params[0];
  ndim=vfield->ndim();
  domain=* (int32_t *) params[1];

  err=vfield->v(domain, t, x, v);
  if (err != 0) return err;

  jmat=new real[ndim*ndim];
  err=vfield->jmat(domain, t, x, jmat);

  v=v+ndim;
  x=x+ndim;
  for (int i=0; i<ndim; i++) {
    v[i]=0;
    for (int j=0; j<ndim; j++) {
      k=j+i*ndim;
      v[i]+=jmat[k]*x[j];
    }
  }

  delete [] jmat;

  //printf("%d %g (%g, %g): (%g, %g)\n", domain, t, x[0], x[1], v[0], v[1]);

  return err;
}

template int ctraj_deriv_L<float>(float t, float *x, float *v, void *param);

ctraj_vfield_base<float> *ctraj_loader(int &argc, char **argv, double &I0, int32_t &n) {
  void *optarg[20];
  int flag[20];
  int32_t vtype=0;
  int32_t ndim=2;
  int err;
  ctraj_vfield_base<float> *vfield;

  //I0=0;

  optarg[0]=&I0;
  optarg[1]=&n;
  optarg[2]=&vtype;
  optarg[5]=&ndim;

  argc=parse_command_opts(argc, argv, "0NVifL", "%lg%d%d%s%s%d", optarg, flag, OPT_WHITESPACE+2);
  if (argc < 0) {
    fprintf(stderr, "ctraj_loader: error parsing command line\n");
    argc=-argc;
  }

  switch (vtype) {
    case (0):
      vfield=new ctraj_vfield_standard<float>();
      break;
    case (1):
      vfield=new ctraj_vfield_2d<float>();
      break;
    case (2):
      vfield=new ctraj_vfield_anal<float>(ndim);
      break;
    default:
      vfield=new ctraj_vfield_standard<float>();
      break;
  }

  //note: assumes vfield needs at least one mandatory argument 
  //  & main routine needs a least 2 otherwise main will go to the help screen
  if (argc>3) {
    err=vfield->setup(argc, argv);
    if (err<0) {
      fprintf(stderr, "Error encountered while loading velocity field\n");
      if (argc>0) argc=-argc;
    }
  } else {
    return vfield;
  }

  if (flag[3]) {
    I0=vfield->get_tind((char *) optarg[3]);
    delete [] (char *) optarg[3];
  }
  if (flag[4]) {
    n=ceil(vfield->get_tind((char *) optarg[4])-I0);
    delete [] (char *) optarg[4];
  } else if (flag[1]==0) {
    //n=vfield->maxt()-I0;
  }

  return vfield;
}

}

