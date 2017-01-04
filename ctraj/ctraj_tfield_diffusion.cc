#include <math.h>

#include "ctraj_tfield_diffusion.h"

namespace ctraj {
  template <class real>
  int ctraj_tfield_diffusion<real>::init(int32_t np, real sdl, real d) {
     dcoeff=d;
     return this->ctraj_tfield_standard::init1(np, sdl);
  }

  template <class real>
  int ctraj_tfield_diffusion<real>::setup(int argc, char **argv) {
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
  int32_t ctraj_tfield_diffusion<real>::interpolate(int32_t domain, real *loc, int32_t *ind, double *wt, real dt) {
    short hemi;
    real loc1[2]={loc[0], loc[1]};
    int32_t d1;
    sub_1d_type ind0;
    sub_1d_type sub;
    real maxd;		//maximum distance
    real k;		//exp(-k*x^2)
    real norm=0;	//normalization coefficient
    int xind1, xind2;
    int yind1, yind2;
    int nx, ny;
    int32_t m=0;			//count

    nx=this->xgrid->n_data();
    ny=this->ygrid->n_data();

    k=1./dt/dcoeff/4;
    maxd=sqrt(-alog(cutoff)/k);

    //this hemisphere:
    hemi=2*domain-1;
    hemi=metric->fix(hemi, loc1);
    d1=(hemi+1)/2;

    xind1=ceil(this->xgrid->interpol(loc1[0]-maxd));
    if (xind1<0) xind1=0;
    xind2=this->xgrid->interpol(loc1[0]-maxd);
    if (xind2>=nx) xind2=nx-1;
    yind1=ceil(this->ygrid->interpol(loc1[1]-maxd));
    if (yind1<0) yind1=0;
    yind2=this->ygrid->interpol(loc1[1]-maxd);
    if (yind2>=ny) yind2=ny-1;
    
    for (int i=xind1; i<xind2; i++) {
      real dx2;
      this->xgrid->get(dx2, i);
      dx2-=loc1[0];
      dx2*=dx2;
      for (int j=yind1; j<yind2; j++) {
        real dy2;
        this->ygrid->get(dy2, i);
        dy2-=loc1[1];
        dy2*=dy2;
	if (dx2+dy2<=maxd*maxd) {
          wt[m]=exp(-k*(dx2+dy2));
	  norm+=wt[m];
	  sub=ny*j+i;
	  this->map_map->get_1d(ind0, sub);
	  ind[m]=ind0+nmap*d1;
	  m++;
	}
      }
    }
        
    //the other hemisphere:
    metric->swap(loc1);
    hemi=-hemi;
    d1=(hemi+1)/2;

    xind1=ceil(this->xgrid->interpol(loc1[0]-maxd));
    if (xind1<0) xind1=0;
    xind2=this->xgrid->interpol(loc1[0]-maxd);
    if (xind2>=nx) xind2=nx-1;
    yind1=ceil(this->ygrid->interpol(loc1[1]-maxd));
    if (yind1<0) yind1=0;
    yind2=this->ygrid->interpol(loc1[1]-maxd);
    if (yind2>=ny) yind2=ny-1;
    
    for (int i=xind1; i<xind2; i++) {
      real dx2;
      this->xgrid->get(dx2, i);
      dx2-=loc1[0];
      dx2*=dx2;
      for (int j=yind1; j<yind2; j++) {
        real dy2;
        this->ygrid->get(dy2, i);
        dy2-=loc1[1];
        dy2*=dy2;
	if (dx2+dy2<=maxd*maxd) {
          wt[m]=exp(-k*(dx2+dy2));
	  norm+=wt[m];
	  sub=ny*j+i;
	  this->map_map->get_1d(ind0, sub);
	  ind[m]=ind0+nmap*d1;
	  m++;
	}
      }
    }
    return m
  }

} //end namespace ctraj


