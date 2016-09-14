#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include "error_codes.h"
#include "parse_command_opts.h"

#include "time_class.h"
#include "peteys_tmpl_lib.h"

#include "coordtran.h"
#include "tracer_anal.h"
#include "ctraj_defaults.h"
#include "av.h"

#define MAXN 10000
#define MAXLL 200
#define TLEN 30

using namespace libpetey;
using namespace libsparse;
using namespace ctraj;

int main(int argc, char **argv) {

  FILE *fs;
  char *tfile;
  char *datefile;
  char *measurement_file;
  char *outfile;
  char c;
  long ncon;

  long fsize;

  float **q;

  int32_t ngrid=NGRID_Q;
  int32_t np;
  int32_t nall;			//total number of tracer fields

  ind_type i0, N;		//start point, number of field to interpolate

  time_class *t;
  char tstring[TLEN];
  char line[MAXLL];

  //for reading and storing the samples:
  long nsamp;
  meas_data *samp;

  float *qint;		//interpolated vector at approximate lead time

  //for calculating lead times:
  float ind;
  float window;
  time_class tf;		//reference date
  time_class t1, t2;		//measurement time window
  int hemi;			//hemisphere
  double l2;			//for truncating matrix array
  int32_t dwid=TFIELD_WIDTH;	//width of date field

  //measurement data within window:
  meas_data *samp_w;
  long nsamp_w;

  int order=2;			//order of regression
  float *coeff;			//regression coefficients
  float qval_p;			//tracer value to various powers

  int err=0;

  void *optarg[20];
  int flag[20];

  i0=0;
  N=-1;

  //parse the command line arguments:
  optarg[0]=&dwid;
  optarg[1]=&order;
  optarg[2]=&i0;
  optarg[3]=&N;
  argc=parse_command_opts(argc, argv, "do0Nif?-+", "%d%d%d%d%s%s%%%", optarg, flag, OPT_WHITESPACE);
  hemi=flag[8]-flag[7];

  if (argc<5 || flag[6]) {
    printf("\n");
    printf("usage: proxy_tracer [-O i0] [-N N] [-o order] [-n ngrid] \n");
    printf("                  tfile dates measurements window outfile\n");
    printf("\n");
    printf("where:\n");
    printf("  matfile      = binary file containing tracer field\n");
    printf("  dates        = ASCII file containing dates corresponding to each field\n");
    printf("  measurements = ASCII file containing measurements and locations\n");
    printf("  window       = measurement window in days\n");
    printf("  outfile      = binary file containing interpolated tracer field\n");
    printf("\n");
    printf("  i0           = index for start point [%d]\n", i0);
    printf("  N            = interpolate this number of fields\n");
    printf("  order        = order of fit [%d]\n", order);
    printf("  n            = grid points per side [%d]\n", ngrid);
    printf("\n");
    return -1;
  }

  tfile=argv[1];
  datefile=argv[2];
  measurement_file=argv[3];
  sscanf(argv[4], "%f", &window);
  outfile=argv[5];

  //read in the array of tracers:
  fprintf(stderr, "Reading file: %s\n", tfile);
  fs=fopen(tfile, "r");
  fread(&np, sizeof(np), 1, fs);
  fseek(fs, 0, SEEK_END);
  fsize=ftell(fs);

  if ((fsize-sizeof(np)) % (sizeof(float)*np) != 0) {
    fprintf(stderr, "Warning: not an even number of tracer fields in file, %s\n", tfile);
  }

  nall=(fsize-sizeof(np))/sizeof(float)/np;
  q=new float*[nall];

  q[0]=new float[nall*np];
  for (int i=1; i<nall; i++) q[i]=q[0]+i*np;

  fseek(fs, sizeof(np), SEEK_SET);

  fread(q[0], sizeof(float), np*nall, fs);
  fprintf(stderr, "%d tracer fields read in\n", nall);

  fclose(fs);

  //calculate grid size from total number of points:
  ngrid=2*(int) sqrt(np/M_PI/2);

  //read in the dates:
  t=new time_class[nall];

  fs=fopen(datefile, "r");
  //fgets(line, MAXLL, fs);		//current versions of "tracer"
					//do NOT print out a header...

  for (long i=0; i<nall; i++) {
    int ind;
    fgets(line, MAXLL, fs);
    sscanf(line, "%d %s", &ind, tstring);
    t[i].read_string(tstring);
  }
  fclose(fs);

  //if time grids are in reverse order, just flip the whole business around:
  if (t[0] > t[nall-1]) {
    time_class tswap;
    float *qswap;
    for (int i=0; i<nall/1; i++) {
      tswap=t[i];
      t[i]=t[nall-i-1];
      t[nall-i-1]=tswap;
      qswap=q[i];
      q[i]=q[nall-i-1];
      q[nall-i-1]=qswap;
    }
  }

  //check and set ranges for time grids:
  if (flag[4]) {		//use initial date to find start index
    t1=(char *) optarg[4];
    i0=bin_search(t, nall, t1, -1);
  }
  if (flag[5]) {		//use final date to find number of grids
    t2=(char *) optarg[5];
    N=bin_search(t, nall, t2, -1)-i0+1;
  }

  //check to make sure time indices are within bounds:
  if (i0<0) i0=0;
  if (N <= 0 || i0+N>nall) N=nall-i0+1;

  //need measurement window as margin since we are interpolating from the 
  //proxy at the measurement times and locations:
  t1=t[i0];
  t1.add(-window);
  if (t1 < t[0]) {
    t1=t[0];
    t1.add(window);
    i0=ceil(interpolate(t, nall, t1, -1));
  }

  tf=t[N+i0-1];
  tf.add(window);
  if (tf > t[nall-1]) {
    tf=t[nall-1];
    tf.add(-window);
    N=bin_search(t, nall, tf, -1)-i0+1;
  }

  //read in measurements:
  printf("Reading in measurements from, %s\n", measurement_file);
  samp=read_meas(measurement_file, &nsamp, dwid);

  qint=new float[np];
  fs=fopen(outfile, "w");
  fwrite(&np, sizeof(np), 1, fs);
  for (long i=i0; i<N+i0; i++) {
    //calculate lead times:
    tf=t[i];
    t1=tf;
    t1.add(-window);
    t2=tf;
    t2.add(window);

    t1.write_string(tstring);
    fprintf(stderr, "%d: interpolating tracer between %s ", i, tstring);
    t2.write_string(tstring);
    fprintf(stderr, "and %s\n", tstring);

    //select measurements:
    samp_w=select_meas(t1, t2, samp, nsamp, &nsamp_w, hemi);

    if (nsamp_w<=order) {
      fprintf(stderr, "proxy_tracer: insufficient measurements (%d) found in window\n");
      fprintf(stderr, "         skipping...\n");
      err=SAMPLE_COUNT_MISMATCH;
      continue;
    }

    coeff=proxy_tracer(q, ngrid, t, nall, samp_w, nsamp_w, order);
/*
    printf("coefficients: ");
    for (int j=0; j<=order; j++) printf(" %g", coeff[j]);
    printf("\n");
*/
    //calculate interpolated tracer field:
    for (int j=0; j<np; j++) {
      qint[j]=0;
      qval_p=1;
      for (int k=0; k<=order; k++) {
        qint[j]+=coeff[k]*qval_p;
        qval_p*=q[i][j];
      }
    }

    //output final, interpolated field:
    t[i].write_string(tstring);
    printf("%d %s\n", i-i0, tstring);
    //printf("%d: writing vector of length %d\n", i, np);
    fwrite(qint, sizeof(float), np, fs);

    delete [] samp_w;
  }

  delete [] qint;

  //qvec=tracer_multiply(matall+i0, nall-i0, q0);

  //fwrite(qvec[0], sizeof(float), n*(nall-i0), fs);
  fclose(fs);

  delete [] q[0];
  delete [] q;
  delete [] t;
  delete [] samp;

  return err;
}

