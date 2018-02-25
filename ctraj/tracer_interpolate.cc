#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include "error_codes.h"
#include "time_class.h"
#include "parse_command_opts.h"
#include "read_ascii_all.h"

#include "ctraj_defaults.h"

#include "tracer_anal.h"

#define TLEN 30

using namespace libpetey;
using namespace ctraj;

int main(int argc, char **argv) {
  FILE *docfs=stderr;

  FILE *fs;
  char *tfile;
  char *datefile;
  char *measurement_file;
  char *outfile;
  char c;
  size_t ncon;

  int32_t nt;

  time_class *t;
  char tstring[TLEN];
  char *line;

  //for reading and storing the samples:
  long nsamp;
  meas_data *samp;		//read directly from measurement file
  meas_data *samp2;		//interpolates
  meas_data *samp3;		//restricted to specified range
  long nsamp2;
  int32_t i0;
  int32_t N;

  //the array tracers
  float **qall;

  int32_t dwid=TFIELD_WIDTH;
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

  //for generating the interpolation coefficients:
  int32_t np;
  //sub_1d_type nmap;

  np=NGRID_Q;
  i0=0;

  optarg[2]=&i0;
  optarg[3]=&dwid;
  optarg[6]=&N;
  optarg[9]=&thresh;
  optarg[10]=&latmin;
  optarg[11]=&latmax;
  argc=parse_command_opts(argc, argv, "-+0tifNPHGIF?", "%%%d%d%s%s%d%%%g%g%g%", optarg, flag, OPT_WHITESPACE);
  if (argc < 0) exit(411);
  if (flag[0]) hemi=-1; else if (flag[1]) hemi=1;
  if (flag[4]) t0.read_string((char *) optarg[4]);
  if (flag[5]) tf.read_string((char *) optarg[5]);
  cflag=flag[7];
  Hflag=flag[8];
  tflag=flag[9];

  if (argc<4 || flag[12]) {
    int err;
    if (flag[12]) {
      docfs=stdout;
      err=0;
    } else {
      docfs=stderr;
      err=INSUFFICIENT_COMMAND_ARGS;
    }

    fprintf(docfs, "\n");
    fprintf(docfs, "usage: tracer_interpolate [-t dwid] [--] [-+] [-P] [-I latmin] [-F latmax]\n");
    fprintf(docfs, "                          [-0 i0|-i t0] [-N n|-f tf]\n");
    fprintf(docfs, "                          tfile dates measurements\n");
    fprintf(docfs, "\n");
    fprintf(docfs, "where:\n");
    fprintf(docfs, "  tfile         = binary file containing tracer fields\n");
    fprintf(docfs, "  dates         = ASCII file containing dates corresponding to each field\n");
    fprintf(docfs, "  measurements  = ASCII file containing measurement locations\n");
    printf("\n");
    fprintf(docfs, "options:\n");
    fprintf(docfs, "  -I   select measurements above this latitude\n");
    fprintf(docfs, "  -F   select measurements below this latitude\n"); 
    ctraj_optargs(docfs, "t-+if0NPH?", 1);
    printf("\n");
    return err;
  }

  tfile=argv[1];
  datefile=argv[2];
  measurement_file=argv[3];

  //read in the array of tracer fields:
  //nmap=calc_nmap(np);
  fprintf(docfs, "Reading file: %s\n", tfile);
  fs=fopen(tfile, "r");
  fread(&nvar, sizeof(nvar), 1, fs);
  np=2*(int32_t) sqrt(nvar/M_PI/2);
  fseek(fs, 0, SEEK_END);
  fsize=ftell(fs);
  nt=(fsize-sizeof(nvar))/(nvar*sizeof(float));
  assert((fsize-sizeof(nvar))%(nvar*sizeof(float)) == 0);

  fprintf(docfs, "Found %d records in %s\n", nt, tfile);

  qall=new float *[nt];
  qall[0]=new float [nt*nvar];
  for (int32_t i=1; i<nt; i++) qall[i]=qall[0]+i*nvar;
  fseek(fs, sizeof(nvar), SEEK_SET);
  fread(qall[0], sizeof(float), nt*nvar, fs);
  fclose(fs);
/*
  for (int i=0; i<nt; i++) {
	  for (int j=0; j<nmap*2; j++) printf("%g ", qall[i][j]);
	  printf("\n");
  }
  */

  //read in the dates:
  fprintf(docfs, "Reading dates from %s\n", datefile);
  t=new time_class[nt];

  fs=fopen(datefile, "r");
  //fgets(line, MAXLL, fs);		//"tracer" does not print a header

  for (int32_t i=0; i<nt; i++) {
    int32_t ind;
    line=fget_line(fs);
    sscanf(line, "%d %s", &ind, tstring);
    t[i].read_string(tstring);
    delete [] line;
  }
  fclose(fs);

  //normally we have to choose between date and index specification
  //(assuming user might try to use both)
  //but here we can do something more interesting:
  //indices select tracer fields, dates select measurements...
  if (i0 < 0) i0=0;
  if (flag[6]==0) N=nt-i0-1;
  if (N >= nt-i0) N=nt-i0-1;

  if (flag[4]==0) t0=t[i0];
  if (flag[5]==0) tf=t[N+i0];
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

  tracer_interp(qall+i0, t+i0, N, np, samp2, nsamp2);

  if (Hflag) {
    for (long i=0; i<nsamp2; i++) {
      printf("%14.7g\n", (samp2[i].q-samp3[i].q)/samp3[i].qerr);
    }
  }

  if (cflag) {
    //inefficient ... who cares...
    if (Hflag==0) {
      for (long i=0; i<nsamp2; i++) {
        samp2[i].t.write_string(tstr);
        printf("%23s %9.3f %9.3f %14.7g %14.7g\n", tstr, 
		      samp2[i].lon, samp2[i].lat, samp3[i].q, samp2[i].q);
      }
    }
    //calculate averages:
    ave1=0;
    ave2=0;
    int nbad=0;
    //and bias:
    float bias=0;
    float rms=0;
    for (long i=0; i<nsamp2; i++) {
      //if (tflag && (samp2[i].q <= thresh || samp3[i].q <= thresh)) {
      if (tflag && samp3[i].q <= thresh) {
	nbad++;
      } else {
        ave1+=samp2[i].q;
        ave2+=samp3[i].q;
	bias+=samp2[i].q-samp3[i].q;
      }
    }
    ave1/=nsamp2-nbad;
    ave2/=nsamp2-nbad;
    bias/=nsamp2-nbad;
    //calculate covariance and standard deviations:
    cov=0;
    var1=0;
    var2=0;
    int fac2=0;
    //and rms
    for (long i=0; i<nsamp2; i++) {
      //if (tflag==0 || (samp2[i].q > thresh && samp3[i].q > thresh)) {
      //zero values for ozone sonde data are bad for sure:
      if (tflag==0 || (samp3[i].q > thresh)) {
        float diff=samp2[i].q-samp3[i].q-bias;
        diff1=samp2[i].q-ave1;
        diff2=samp3[i].q-ave2;
        cov+=diff1*diff2;
        var1+=diff1*diff1;
        var2+=diff2*diff2;
	rms+=diff*diff;
	if (samp2[i].q/samp3[i].q < 2 && samp3[i].q/samp2[i].q < 2) fac2++;
      }
    }
    float std1=sqrt(var1/(nsamp2-nbad-1));	//interpolates
    float std2=sqrt(var2/(nsamp2-nbad-1));	//measurements
    rms=sqrt(rms/(nsamp2-nbad-1));
    printf("r    = %g\n", cov/std1/std2/(nsamp2-nbad-1));
    printf("bias = %g\n", bias);
    printf("rms  = %g\n", rms);
    printf("fractional bias = %g\n", bias/std2);
    printf("normalized rms  = %g\n", rms/std2);
    printf("fac2 = %g\n", 1.*fac2/(nsamp2-nbad));
  } else if (Hflag==0) {
    write_meas(samp2, nsamp2, stdout);
  }

  delete [] qall[0];
  delete [] qall;

  delete [] t;
  delete [] samp;
  delete [] samp2;

  if (cflag || Hflag) delete [] samp3;
}

