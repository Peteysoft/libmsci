#include <assert.h>
#include <math.h>

#include "peteys_tmpl_lib.h"

#include "nr.h"

#include "tcoord_defs.h"
#include "boundary_list_g.h"

#define PI 3.141592657

extern simple<float> *global_xgrid;
extern simple<float> *global_ygrid;
extern dependent_swap<float> *global_u;
extern dependent_swap<float> *global_v;

extern void derivs(float tind, float x[], float dxdt[], long n);

//default constructor:
boundary_list_g::boundary_list_g() {
  xgrid_N=NULL;
  ygrid_N=NULL;
  xgrid_S=NULL;
  ygrid_S=NULL;
  zgrid=NULL;
  tgrid=NULL;

  u_N=NULL;
  v_N=NULL;
  u_S=NULL;
  v_S=NULL;

  N_fs=NULL;
  S_fs=NULL;

  result=NULL;

  x=NULL;
  y=NULL;
  hemi=NULL;
  n=0;

  tind=0;
  tstep=0;
  nrk=0;
  min_spac=0;
  thresh_arc=0;
  wrap_flag=0;
}

boundary_list_g::boundary_list_g(char *nfile, char *sfile) {
  simple<float> *check;
  simple<time_class> *checkt;
  long loc, dum;

  printf("Opening: %s\n", nfile);
  N_fs=fopen(nfile, "r");

  //get the N. hemisphere velocity field:
  ndata.read(N_fs);
  loc=ndata.search_var("xgrid", dum);
  xgrid_N=(simple<float> *) ndata.get_var(loc);
  loc=ndata.search_var("ygrid", dum);
  ygrid_N=(simple<float> *) ndata.get_var(loc);
  loc=ndata.search_var("zgrid", dum);
  zgrid=(simple<float> *) ndata.get_var(loc);
  loc=ndata.search_var("time", dum);
  tgrid=(simple<time_class> *) ndata.get_var(loc);
  loc=ndata.search_var("u", dum);
  u_N=(dependent_swap<float> *) ndata.get_var(loc);
  loc=ndata.search_var("v", dum);
  v_N=(dependent_swap<float> *) ndata.get_var(loc);

  //get the S. hemisphere velocity field:
  printf("Opening: %s\n", sfile);
  S_fs=fopen(sfile, "r");
  sdata.read(S_fs);
  loc=sdata.search_var("u", dum);
  u_S=(dependent_swap<float> *) sdata.get_var(loc);
  loc=sdata.search_var("v", dum);
  v_S=(dependent_swap<float> *) sdata.get_var(loc);

  //check the gridding of the two files to make sure that they agree:
  loc=sdata.search_var("xgrid", dum);
  xgrid_S=(simple<float> *) sdata.get_var(loc);
  loc=sdata.search_var("ygrid", dum);
  ygrid_S=(simple<float> *) sdata.get_var(loc);
  loc=sdata.search_var("zgrid", dum);
  check=(simple<float> *) sdata.get_var(loc);
  assert(*check == *zgrid);
  loc=sdata.search_var("time", dum);
  checkt=(simple<time_class> *) sdata.get_var(loc);
  assert(*checkt == *tgrid);

  x=NULL;
  y=NULL;
  hemi=NULL;
  n=0;

  result=NULL;
  nrk=0;
  tstep=0;

  wrap_flag=1;

}

void boundary_list_g::set(char *t0, float tcoarse, long nfine, float tarc,
			float minspc) {
  float *temp;

  tind=tgrid->interp(t0);
  printf("Time index: %f\n", tind);

  tstep=tcoarse;
  nrk=nfine;
  thresh_arc=tarc;
  min_spac=minspc;

  if (result != NULL) {
    delete [] result[0];
    delete [] result;
  }

  temp=new float[(nrk+1)*2];
  result=new float *[nrk+1];
  for (long i=0; i<=nrk; i++) result[i]=&temp[2*i];
}

long boundary_list_g::init_circle(float x0, float y0, float r) {
  float ang;

  n=(long) (2*PI/thresh_arc);
  x=new float[n+1*wrap_flag];
  y=new float[n+1*wrap_flag];

  for (long i=0; i<n; i++) {
    ang=i*thresh_arc;
    //no attempt to apply metric corrections:
    x[i]=x0+r*cos(ang);
    y[i]=y0+r*sin(ang);
  }

  return n;
}


void boundary_list_g::wrap_on() {
  float *xnew;
  float *ynew;
  short *hnew;

  if (!wrap_flag) {
    wrap_flag=1;
    if (n != 0) {
      xnew=new float[n+1];
      ynew=new float[n+1];
      hnew=new short[n+1];
      for (long i=1; i<=n; i++) {
        xnew[i]=x[i];
        ynew[i]=y[i];
	hnew[i]=hemi[i];
      }
      delete [] x;
      delete [] y;
      delete [] hemi;
      x=xnew;
      y=ynew;
      hemi=hnew;
    }
  }
}

/*
void boundary_list_g::wrap_off() {
//  if (wrap_flag) {
    wrap_flag=0;
    if (n!=0) {
      for (long i=0; i<n; i++) {
        x[i]=x[i+1];
        y[i]=y[i+1];
      }
    }
  }
}
*/


boundary_list_g::~boundary_list_g() {
//  delete x;
//  delete y;
  delete [] hemi;

  //delete u_N;
  //delete v_N;
  //delete u_S;
  //delete v_S;
  //delete xgrid_N;
  //delete ygrid_N;
  //delete xgrid_S;
  //delete ygrid_S;
  //delete zgrid;
  //delete tgrid;

  fclose(N_fs);
  fclose(S_fs);
}

time_class boundary_list_g::advance() {
  short hold;
  float x0[2];
  float tstepfine=tstep/nrk;
  long tindl;
  double frac;
  time_class t0, t1;

  //pass the grids to the global variables:
  global_xgrid=xgrid_N;
  global_ygrid=ygrid_N;
  global_u=u_N;
  global_v=v_N;
  hold=1;

  for (long i=0; i<n; i++) {
    tcoord_fix(x[i], y[i], hemi[i]);
    if (hemi[i] != hold) {
      hold=hemi[i];
      //pass the velocity fields to the global variables:
      if (hold == 1) {
        global_u=u_N;
        global_v=v_N;
        global_xgrid=xgrid_N;
        global_ygrid=ygrid_N;
      } else {
        global_u=u_S;
        global_v=v_S;
        global_xgrid=xgrid_S;
        global_ygrid=ygrid_S;
      }
    }
    x0[0]=x[i];
    x0[1]=y[i];
    rk_dumb((float) tind, x0, 2L, tstepfine, nrk, (float **) result, &derivs);
    x[i]=result[nrk][0];
    y[i]=result[nrk][1];
  }

  tind+=tstep;

  tindl=(long) tind;
  frac=tind-(double) tindl;

  tgrid->get(t0, tindl);
  tgrid->get(t1, tindl+1);

  t1=t1-t0;
  t1=t1*frac;
  t1=t0+t1;

  return t1;
}

long boundary_list_g::fix() {
  long nn=n+wrap_flag;
  float s[nn];		//the parametric distance along the curve
  float dx0, dy0;
  float dx2[nn];	//2nd derivative of x wrt s
  float dy2[nn];	//2nd derivative of y wrt s
  float dx2av, dy2av;
  float dels;		//the distance between a pair of points
  float dstotal;
  float dangle;		//the radians of arc (^2) between a pair of points
  float curv2;		//the curvature squared
  int nnew[nn-1];	//number of new points to insert between each pair of points
  long totaln, dn;	//final total number of points
  float * xnew;
  float * ynew;
  short * hnew;
  long j, k;
  long offset;		//for removing redundant nodes

  if (wrap_flag) {
    x[n]=x[0];
    y[n]=y[0];
    hemi[n]=hemi[0];
  }

  s[0]=0;
  //we work in the N. hemisphere:
  if (hemi[0] == -1) tcoord_N2S(x[0], y[0]);
  offset=0;
  //printf("nn=%d\n", nn);
  for (long i=1; i<nn; i++) {
    if (offset != 0) {
      //because of any redundant nodes, all the nodes must be shifted:
      j=i-offset;
      x[j]=x[i];
      y[j]=y[i];
      hemi[j]=hemi[i];
    } else {
      j=i;
    }
    if (hemi[j] == -1) tcoord_N2S(x[j], y[j]);
    dels=tcoord_ds2(x[j-1], y[j-1], x[j], y[j]);
    s[j]=s[j-1]+sqrt(dels);
    //printf("%g %g\n", 
    //		sqrt((x[j-1]-x[j])*(x[j-1]-x[j])+(y[j-1]-y[j])*(y[j-1]-y[j])), 
    //		sqrt(dels));
    //check for nodes made redundant by round-off:
    if (s[j]-s[j-1] <= 0) {
      printf("Warning: removing %dth node made redundant by round-off\n", i);
      offset++;
    }
  }
  n-=offset;
  nn-=offset;

  //do a cubic spline interpolation on both variables:
  if (wrap_flag) {
    dx0=(x[n-1]-x[1])/(s[1]+s[n]-s[n-1]);
    dy0=(y[1]-y[n-1])/(s[1]+s[n]-s[n-1]);
//    printf("dx0ds=%g; dy0ds=%g\n", dx0, dy0);
  } else {
    dx0=0;
    dy0=0;
  }
  spline(s-1, x-1, nn, dx0, dx0, dx2-1);
  spline(s-1, y-1, nn, dy0, dy0, dy2-1);

  //calculate the angle traced out by each pair of points:
  //and use to figure out how many new points to add...
  totaln=1;
  dangle=0;
  dstotal=0;
//  printf("i, s, dx2ds2, dy2ds2, dangle, #new\n");
  for (long i=1; i<nn; i++) {
    dels=s[i]-s[i-1];
    dstotal+=dels;
    assert(dels>0);
    //if the point spacing is less than the threshold, either
    //remove next point or add no new:
    if (dstotal < min_spac) {
      if (2*dstotal < min_spac) {
        nnew[i-1]=-1;		//remove point
      } else {
        nnew[i-1]=0;
      }
    } else {
      dstotal=0;
      dx2av=(dx2[i-1]+dx2[i])/2;
      dy2av=(dy2[i-1]+dy2[i])/2;
      dangle+=sqrt(dx2av*dx2av+dy2av*dy2av)*dels;

      //printf("r=%g\n", dels/dangle);

      //if the accumulated angle of arc is greater than the threshold,
      //add new points, otherwise, remove the next point...
      if (dangle*2 < thresh_arc) {
        nnew[i-1]=-1;
      } else {
        nnew[i-1]=(int) (dangle/thresh_arc-0.5);
        //dangle=dangle-(nnew[i-1]+0.5)*thresh_arc;
	dangle=0;
      }
    }
    totaln+=nnew[i-1]+1;
    //printf("%d, %g, %g, %g, %g, %d\n", i, s[i], dx2[i], dy2[i], dangle, nnew[i-1]);
  }

  xnew=new float[totaln];
  ynew=new float[totaln];

  j=1;			//index of new nodes
  xnew[0]=x[0];		//always keeps the first point
  ynew[0]=y[0];

  k=0;			//start interpolating from this node
  for (long i=1; i<nn; i++) {
    if (nnew[i-1] > 0) {
      dstotal=s[i]-s[k];
      dels=dstotal/(nnew[i-1]+1);
      if (2*dels < min_spac) {
        dn=nnew[i-1]-(short) (dstotal/min_spac*2);
        nnew[i-1]-=dn;
	totaln-=dn;
        dels=dstotal/(nnew[i-1]+1);
      }
      //add in new points interpolated along the spline:
      for (long m=1; m <= nnew[i-1]; m++) {
        splint(s-1, x-1, dx2-1, nn, s[k]+m*dels, &xnew[j]);
        splint(s-1, y-1, dy2-1, nn, s[k]+m*dels, &ynew[j]);
        //printf("%g %g\n", xnew[j], ynew[j]);
	j++;
      }
      xnew[j]=x[i];
      ynew[j]=y[i];
      k=i;
      j++;
    } else if (nnew[i-1] == 0) {
      //no new points to add:
      xnew[j]=x[i];
      ynew[j]=y[i];
      k=i;
      j++;
    }
  }

  assert(j==totaln);

  delete [] x;
  delete [] y;
  x=xnew;
  y=ynew;

  delete [] hemi;
  hemi=new short[totaln];
  for (long i=0; i<totaln; i++) hemi[i]=1;

  n=totaln-wrap_flag;

  return n;

}

long boundary_list_g::read(FILE *fs) {
  short yy, mon, dd, hh, min;
  float sec;
  time_class start;
  long nread;

  if (x != NULL) delete [] x;
  if (y != NULL) delete [] y;
  if (hemi != NULL) delete [] hemi;

  nread=fread(&yy, sizeof(yy), 1, fs);
  nread+=fread(&mon, sizeof(mon), 1, fs);
  nread+=fread(&dd, sizeof(dd), 1, fs);
  nread+=fread(&hh, sizeof(hh), 1, fs);
  nread+=fread(&min, sizeof(min), 1, fs);
  nread+=fread(&sec, sizeof(sec), 1, fs);

  start.init(yy, mon, dd, hh, min, sec);

  tind=tgrid->interp(start);

  nread+=fread(&n, sizeof(n), 1, fs);

  x=new float[n+wrap_flag];
  y=new float[n+wrap_flag];
  hemi=new short[n+wrap_flag];

  nread+=fread(x, sizeof(float), n, fs);
  nread+=fread(y, sizeof(float), n, fs);

  //convert from lon-lat:
  for (long i=0; i<n; i++) {
    hemi[i]=0;
    tcoord_from_lonlat(x[i], y[i], hemi[i], x[i], y[i]);
  }

  return nread;

}

long boundary_list_g::write(FILE *fs) {
  float lon[n], lat[n];
  short yy, mon, dd, hh, min;
  float sec;
  time_class tcur, t0, t1;
  long nwrit;
  long tindl;
  double frac;

  tindl=(long) tind;
  frac=tind-(double) tindl;

  tgrid->get(t0, tindl);
  tgrid->get(t1, tindl+1);

  t1=(t1-t0);
  t1=t1*frac;
  tcur=t0+t1;
  tcur.get_fields(yy, mon, dd, hh, min, sec);

  nwrit=fwrite(&yy, sizeof(yy), 1, fs);
  nwrit+=fwrite(&mon, sizeof(mon), 1, fs);
  nwrit+=fwrite(&dd, sizeof(dd), 1, fs);
  nwrit+=fwrite(&hh, sizeof(hh), 1, fs);
  nwrit+=fwrite(&min, sizeof(min), 1, fs);
  nwrit+=fwrite(&sec, sizeof(sec), 1, fs);

  nwrit+=fwrite(&n, sizeof(n), 1, fs);

  //convert to lon-lat:
  for (long i=0; i<n; i++) {
    tcoord_2lonlat(x[i], y[i], hemi[i], lon[i], lat[i]);
    //printf("%g %g %d %g %g\n", x[i], y[i], hemi[i], lon[i], lat[i]);
  }

  nwrit+=fwrite(lon, sizeof(float), n, fs);
  nwrit+=fwrite(lat, sizeof(float), n, fs);

  return nwrit;

}

long boundary_list_g::print(FILE *fs) {
  for (long i=0; i<n; i++) {
    printf("%10.6g %10.6g %d\n", x[i], y[i], hemi[i]);
  }
  return n;
}

