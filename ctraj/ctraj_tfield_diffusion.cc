#include <math.h>

#include "parse_command_opts.h"

#include "ctraj_defaults.h"
#include "ctraj_tfield_diffusion.h"

namespace ctraj {
  template <class real>
  ctraj_tfield_diffusion<real>::ctraj_tfield_diffusion() {
    //do I need to do this? I can never remember...
    this->xgrid=NULL;
    this->ygrid=NULL;
    this->map_map=NULL;
    this->tracer=NULL;
    this->inverse_map=NULL;
    dcoeff=DIFFUSION;
    cutoff=GAUSS_CUTOFF;
  }

  template <class real>
  ctraj_tfield_diffusion<real>::~ctraj_tfield_diffusion() {
    //pretty sure I don't need to do anything here...
  }

  template <class real>
  int ctraj_tfield_diffusion<real>::init(int32_t np, real sdl, real d, real co) {
     dcoeff=d;
     cutoff=co;
     return this->ctraj_tfield_standard<real>::init1(np, sdl);
  }

  template <class real>
  int ctraj_tfield_diffusion<real>::setup(int argc, char **argv) {
    void *optarg[10];
    int flag[10];

    real sidelength=SIDELENGTH_Q;
    int32_t np=NGRID_Q;
    real d=DIFFUSION;
    real co=GAUSS_CUTOFF;
    int nleft;
    sub_1d_type ind;

    optarg[0]=&sidelength;
    optarg[1]=&np;
    optarg[2]=&d;
    optarg[3]=&co;

    nleft=parse_command_opts(argc, argv, "rnCG", "%g%d%g%g", optarg, flag, OPT_WHITESPACE+2);

    if (nleft<0) {
      fprintf(stderr, "ctraj_tfield_diffusion: command option parse error\n");
    }

    init(np, sidelength, d, co);

    return nleft;

  }

  template <class real>
  int32_t ctraj_tfield_diffusion<real>::nwt() {
    //maximum number of weights: all the grids
    return this->nmap;
  }

  template <class real>
  int32_t ctraj_tfield_diffusion<real>::interpolate(int32_t domain, real *loc, int32_t *ind, double *wt, real dt) {
    short hemi;
    real loc1[2]={loc[0], loc[1]};
    real loc2[2];
    int32_t d1;
    sub_1d_type ind0;
    sub_1d_type sub;
    real maxd;		//maximum distance
    real k;		//exp(-k*x^2)
    real norm=0;	//normalization coefficient
    int xind1, xind2;
    int yind1, yind2;
    int nx, ny;
    real c[2];
    real ds2;
    int32_t m=0;			//count

    nx=this->xgrid->nel();
    ny=this->ygrid->nel();

    k=1./dt/dcoeff/4;
    maxd=sqrt(-log(cutoff)/k);

    //this hemisphere:
    hemi=2*domain-1;
    hemi=this->metric->fix(hemi, loc1);
    d1=(hemi+1)/2;

    this->metric->mcoef2(loc1, c);
    xind1=ceil(this->xgrid->interp(loc1[0]-sqrt(c[0])*maxd));
    if (xind1<0) xind1=0;
    xind2=this->xgrid->interp(loc1[0]+sqrt(c[0])*maxd);
    if (xind2>=nx) xind2=nx-1;
    yind1=ceil(this->ygrid->interp(loc1[1]-sqrt(c[1])*maxd));
    if (yind1<0) yind1=0;
    yind2=this->ygrid->interp(loc1[1]+sqrt(c[1])*maxd);
    if (yind2>=ny) yind2=ny-1;
    
    for (int i=xind1; i<xind2; i++) {
      this->xgrid->get(loc2[0], i);
      for (int j=yind1; j<yind2; j++) {
        this->ygrid->get(loc2[1], j);
	ds2=this->metric->ds2(loc2, loc1);
	if (ds2<=maxd*maxd) {
          wt[m]=exp(-k*ds2);
	  norm+=wt[m];
	  sub=ny*j+i;
	  this->map_map->get_1d(ind0, sub);
	  ind[m]=ind0+this->nmap*d1;
	  m++;
	}
      }
    }
        
    //the other hemisphere:
    this->metric->swap(loc1);
    hemi=-hemi;
    d1=(hemi+1)/2;

    xind1=ceil(this->xgrid->interp(loc1[0]-sqrt(c[0])*maxd));
    if (xind1<0) xind1=0;
    xind2=this->xgrid->interp(loc1[0]+sqrt(c[0])*maxd);
    if (xind2>=nx) xind2=nx-1;
    yind1=ceil(this->ygrid->interp(loc1[1]-sqrt(c[1])*maxd));
    if (yind1<0) yind1=0;
    yind2=this->ygrid->interp(loc1[1]+sqrt(c[1])*maxd);
    if (yind2>=ny) yind2=ny-1;
    
    for (int i=xind1; i<xind2; i++) {
      this->xgrid->get(loc2[0], i);
      for (int j=yind1; j<yind2; j++) {
        this->ygrid->get(loc2[1], j);
	ds2=this->metric->ds2(loc1, loc2);
	if (ds2<=maxd*maxd) {
          wt[m]=exp(-k*ds2);
	  norm+=wt[m];
	  sub=ny*j+i;
	  this->map_map->get_1d(ind0, sub);
	  ind[m]=ind0+this->nmap*d1;
	  m++;
	}
      }
    }
    for (int i=0; i<m; i++) wt[i]/=norm;
    return m;
  }

  template <class real>
  void ctraj_tfield_diffusion<real>::help(FILE *fs) {
    fprintf(fs, "Tracer model: 2-D global on azimuthal-equidistant coords\n");
    fprintf(fs, "              with tunable diffusion\n");
    ctraj_optargs(fs, "rnCG", OPT_WHITESPACE+2);
  }

  template class ctraj_tfield_diffusion<float>;

} //end namespace ctraj


