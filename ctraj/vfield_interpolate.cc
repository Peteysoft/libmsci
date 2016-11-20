#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>

#include "error_codes.h"
#include "time_class.h"
#include "parse_command_opts.h"
#include "read_ascii_all.h"

#include "ctraj_defaults.h"

#include "tracer_anal.h"
#include "ctraj_vfield_standard.h"

#define TLEN 30

using namespace libpetey;
using namespace ctraj;

int main(int argc, char **argv) {
  FILE *docfs=stderr;

  FILE *fs;
  char *measurement_file;
  char *outfile;
  char c;
  size_t ncon;

  ctraj_vfield_standard<float> *vfield;

  //for reading and storing the samples:
  long nsamp;
  meas_data *samp;
  meas_data *samp2;
  meas_data *samp3;
  long nsamp2;
  int32_t i0;
  int32_t N;
  int32_t nt; 

  short hemi=0;
  int cflag=0;

  //for calculating correlation coefficient:
  float ave1, ave2;
  float diff1, diff2;
  float var1, var2;
  float cov;

  float thresh;
  float latmin=-90;
  float latmax=90;

  time_class t0;
  time_class tf;
  time_class ttest;
  char tstr[30];

  size_t fsize;
  int32_t nvar;

  void *optarg[20];
  int flag[20];

  int Hflag;		//for generationg histograms
  int tflag;		//we're cheating...

  int64_t page_size=-1;
  int32_t dwid=TFIELD_WIDTH;
  int component=-1;

  i0=0;

  optarg[2]=&i0;
  optarg[3]=&dwid;
  optarg[6]=&N;
  optarg[9]=&thresh;
  optarg[10]=&latmin;
  optarg[11]=&latmax;
  optarg[12]=&page_size;
  argc=parse_command_opts(argc, argv, "-+0tifNPHGIFB?", "%%%d%d%s%s%d%%%g%g%g%ld%", optarg, flag, OPT_WHITESPACE);
  if (argc < 0) exit(411);
  if (flag[0]) hemi=-1; else if (flag[1]) hemi=1;
  if (flag[4]) t0.read_string((char *) optarg[4]);
  if (flag[5]) tf.read_string((char *) optarg[5]);
  cflag=flag[7];
  Hflag=flag[8];
  tflag=flag[9];

  if (argc<4 || flag[13]) {
    int err;
    if (flag[13]) {
      docfs=stdout;
      err=0;
    } else {
      docfs=stderr;
      err=INSUFFICIENT_COMMAND_ARGS;
    }

    fprintf(docfs, "\n");
    fprintf(docfs, "usage: vfield_interpolate [-t dwid] [--] [-+] [-P] [-I latmin] [-F latmax]\n");
    fprintf(docfs, "                          [-0 i0|-i t0] [-N n|-f tf]\n");
    fprintf(docfs, "                          [component] vfileS vfileN measurements\n");
    fprintf(docfs, "\n");
    fprintf(docfs, "where:\n");
    fprintf(docfs, "  component     = u  for zonal wind\n");
    fprintf(docfs, "                  v  for meridional wind\n");
    fprintf(docfs, "                  w  for vertical wind\n");
    fprintf(docfs, "  vfileS        = binary file containing S. hemisphere velocity field\n");
    fprintf(docfs, "  vfileN        = binary file containing N. hemisphere velocity field\n");
    fprintf(docfs, "  measurements  = ASCII file containing measurement locations\n");
    printf("\n");
    fprintf(docfs, "options:\n");
    fprintf(docfs, "  -I   select measurements above this latitude\n");
    fprintf(docfs, "  -F   select measurements below this latitude\n"); 
    ctraj_optargs(docfs, "t-+if0NPHB?", 1);
    printf("\n");
    return err;
  }

  if (argc>4) {
    if (strcmp(argv[1], "u")==0) {
      component=0;
    } else if (strcmp(argv[1], "v")==0) {
      component=1;
    } else if (strcmp(argv[1], "w")==0) {
      component=2;
    }
    argv++;
  }
  measurement_file=argv[3];

  //read in the velocity fields:
  fprintf(docfs, "Reading files: %s %s\n", argv[1], argv[2]);
  vfield=new ctraj_vfield_standard<float>();
  vfield->init(argv+1, page_size, "r");
  nt=vfield->maxt()+1;

  //normally we have to choose between date and index specification
  //(assuming user might try to use both)
  //but here we can do something more interesting:
  //indices select tracer fields, dates select measurements...
  if (i0 < 0) i0=0;
  if (flag[6]==0) N=nt-i0-1;
  if (N >= nt-i0) N=nt-i0-1;

  if (flag[4]==0) t0=vfield->get_t(i0);
  if (flag[5]==0) tf=vfield->get_t(N+i0);
/*
  for (int32_t i=0; i<nt; i++) {
	  t[i].write_string(tstring);
	  printf("%s\n", tstring);
  }
  tf.write_string(tstring);*/
  //printf("%s\n", tstring);
  //printf("i0=%d, N=%d\n", i0, N);

  //read in measurements:
  fprintf(docfs, "Reading in measurements from, %s\n", measurement_file);
  samp=read_meas(measurement_file, &nsamp, dwid);

  //exclude measurements outside the date range:
  samp2=select_meas(t0, tf, samp, nsamp, &nsamp2, hemi);

  if (flag[10] || flag[11]) {
    meas_data *samp4;
    samp4=select_lat_range(samp2, nsamp2, latmin, latmax, &nsamp2);
    delete [] samp2;
    samp2=samp4;
  }

/*  for (long i=0; i<nsamp; i++) {
	  dates[i].write_string(tstring);
	  printf("%s\n", tstring);
  }*/

  //printf("nmap=%d\n", nmap);

  //do the interpolation:
  fprintf(docfs, "Performing interpolation...\n");

  if (cflag || Hflag) {
    samp3=copy_meas(samp2, nsamp2);
  }

  for (int i=0; i<nsamp2; i++) {
    float loc[2];		//location in transformed coords
    float v[3];			//velocity
    int32_t domain;		//hemisphere
    int hemi;
    float r2;			//distance from pole
    double tind;		//time index
    double vunit;		//number of days per time grid
    double q;
    time_class t1, t2;
    loc[0]=samp2[i].lon;
    loc[1]=samp2[i].lat;
    //figure out time stuff:
    tind=vfield->get_tind(samp2[i].t);
    if (tind>=nt-1) {
      t1=vfield->get_t((int) tind-1);
      t2=vfield->get_t((int) tind);
    } else if (tind<=0) {
      t1=vfield->get_t(0);
      t2=vfield->get_t(1);
    } else {
      t1=vfield->get_t((int) tind);
      t2=vfield->get_t((int) tind+1);
    }
    vunit=t2.diff(t1)*SECPERDAY/MPERKM;

    //printf("%g %g\n", loc[0], loc[1]);
    domain=vfield->absolute(-1, loc);
    hemi=domain*2-1;
    //printf("%g %g %g %d\n", loc[0], loc[1], vfield->get_tind(samp2[i].t), domain);
    v[2]=0;
    vfield->v(domain, vfield->get_tind(samp2[i].t), loc, v);
    r2=loc[0]*loc[0]+loc[1]*loc[1];
    switch (component) {
      case (0):
	if (loc[0]!=0 || loc[1]!=0) {
          q=REARTH*sin(sqrt(r2)/REARTH)*(v[1]*loc[0]-v[0]*loc[1])/r2;
	} else {
          q=v[1]*cos(M_PI*samp2[i].lon/180)-v[0]*sin(M_PI*samp2[i].lon/180);
	}
	q/=vunit;
        break;
      case (1):
	if (loc[0]!=0 || loc[1]!=0) {
          q=(v[0]*loc[0]+v[1]*loc[1])/sqrt(r2);
	} else {
          q=v[0]*cos(M_PI*samp2[i].lon/180)+v[1]*sin(M_PI*samp2[i].lon/180);
	}
        q*=-hemi/vunit;
        break;
      case (2):
        samp2[i].q=v[2];
	break;
      default:
	if (loc[0]!=0 || loc[1]!=0) {
          samp2[i].q=REARTH*sin(sqrt(r2)/REARTH)*(v[1]*loc[0]-v[0]*loc[1])/r2;
	} else {
          samp2[i].q=v[1]*cos(M_PI*samp2[i].lon/180)-v[0]*sin(M_PI*samp2[i].lon/180);
	}
	samp2[i].q/=vunit;
	if (loc[0]!=0 || loc[1]!=0) {
          samp2[i].qerr=(v[0]*loc[0]+v[1]*loc[1])/sqrt(r2);
	} else {
          samp2[i].qerr=v[0]*cos(M_PI*samp2[i].lon/180)+v[1]*sin(M_PI*samp2[i].lon/180);
	}
        samp2[i].qerr*=-hemi/vunit;
    }
    if (component>=0) { 
      samp2[i].q=q;
      samp2[i].qerr=0;
    }
  }

  if (cflag) {
    for (int i=0; i<nsamp2; i++) {
      samp2[i].t.write_string(tstr);
      printf("%23s %9.3f %9.3f %14.7lg %14.7lg\n", tstr, 
		      samp2[i].lon, samp2[i].lat, samp3[i].q, samp2[i].q);
    }
    printf("r=%lg\n", correlate_meas(samp2, samp3, nsamp2));
  } else {
    int flag=0;
    if (component<0) flag=1;
    write_meas(samp2, nsamp2, stdout, flag);
  }

  delete [] samp;
  delete [] samp2;

  if (cflag || Hflag) delete [] samp3;
}

