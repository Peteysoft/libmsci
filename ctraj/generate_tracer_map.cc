#include <stdint.h>

#include "sparse.h"
#include "ctraj_tracer_base.h"

int tracer_map_deriv(float t, float *x, float *v, void *param) {
  void **params;
  int32_t domain;
  ctraj_vfield_base *vfield;

  params=(void **) param;
  vfield=(ctraj_vfield_base *) params[0];
  domain=(int32_t) params[1];

  return vfield->v(domain, t, x, v);
}
  

int generate_tracer_map(ctraj_vfield_base *vfield,	//velocity field
		ctraj_tracer_base *tracer, 		//tracer interpolates
		double t,				//initial time
		double h,				//time step
		int32_t nt,				//number of time steps
		sparse<int32_t, float> *map,		//matrix mapping
		int32_t maxninterp) {		//maximum number of interpolation
						//coefficients

  void *param[2];
  int32_t n;			//number of points in tracer field
  int32_t ndim, ndim1;		//number of dimensions
  float **x;			//trajectory
  int32_t domain1, domain2;	//velocity/tracer field domain
  int32_t sub[maxninterp];	//index for interpolation neighbour
  float wt[maxninterp];		//interpolation weight
  int32_t ninterp;		//number of interpolation neighbours

  n=tracer->nel();
  ndim=tracer->ndim();
  ndim1=vfield->ndim();
  assert(ndim==ndim1);

  x=new float *[nt];
  x[0]=new float[ndim*nt];
  for (int32_t i=1; i<=nt; i++) x[i]=x[0]+i*ndim;

  param[0]=vfield;

  for (int32_t i=0; i<n; i++) {
    domain1=tracer->get_loc(i, domain2, x[0]);
    param[1]=&domain1;
    rk_dumb((float) t, x[0], 2L, h, nt, x, (void *) param, &tracer_map_derivs);
    ninterp=tracer->interpolate(domain2, x[nt], sub, wt);
    for (int32_t j=0; j<ninterp; j++) {
      map->add_el(wt[j], i, sub[j]);
    }
  }

  delete x[0];
  delete x;
    
}

