#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include "error_codes.h"
#include "parse_command_opts.h"
#include "time_class.h"
#include "peteys_tmpl_lib.h"
#include "read_ascii_all.h"

#include "tracer_anal.h"
#include "ctraj_defaults.h"
#include "av.h"

//*FLAG* -- fixed-length data structure
#define MAXN 10000

using namespace libpetey;
using namespace libsparse;
using namespace ctraj;

int main(int argc, char **argv) {
  FILE *docfs=stderr;

  FILE *fs;
  char *matfile;
  char *datefile;
  char *measurement_file;
  char *outfile;
  char c;
  sparse_matrix dummy;
  size_t ncon;

  sparse_matrix *matall;
  int32_t nall;

  ind_t m, n;	//size of matrix
  int32_t i0, N;		//start point, number of matrices
  int32_t nev, ncv;		//number of eigenvectors, number of Arnoldi vectors

  int32_t nvar;			//for compatibility with libagf

  time_class *t;

  char tstring[TSTRING_LEN];
  char *line;

  //for reading and storing the samples:
  long nsamp;
  meas_data *samp;

  float *q0;
  float **qvec;

  //for calculating lead times:
  float lead;			//length in days of tracer map
  float lead2;		//days between start of tracer map and measurement window
  float window;
  time_class tf;		//date at end of lead time
  time_class tf2;		//date at end of lead2 time
  time_class t1, t2;		//measurement time window
  double l2;			//for truncating matrix array
  int32_t N2;			//size of matrix array to use
  int32_t N3;			//index for output field

  //measurement data within window:
  meas_data *samp_w;
  long nsamp_w;
  int32_t dwid=TFIELD_WIDTH;

  //date range:
  time_class date0;
  time_class datef;

  int flag[20];
  void *optarg[20];

  //northern or southern hemisphere:
  int hemi=0;

  int err=0;

  int cflag=0;
  int kflag=0;

  //for getting sample statistics:
  float nsamp_ave=0;
  int32_t nfield=0;

  i0=0;
  nall=MAXN;
  N=-1;

  ncv=NARNOLDI;
  nev=NEIG;

  optarg[0]=&dwid;
  optarg[3]=&i0;
  optarg[4]=&N;
  optarg[7]=&ncv;
  optarg[8]=&nev;
  optarg[11]=&lead2;
  argc=parse_command_opts(argc, argv, "d-+0NifAvCKl?", 
		"%d%%%d%d%s%s%d%d%%%g%", optarg, flag, OPT_WHITESPACE);
  if (argc<0) exit(21);
  if (flag[1]) hemi=-1; else if (flag[2]) hemi=1;
  if (flag[5]) date0.read_string((char *) optarg[5]);
  if (flag[6]) datef.read_string((char *) optarg[6]);
  cflag=flag[9];
  kflag=flag[10];

  if ((argc<6 || (cflag==0 && argc < 7)) || flag[12]) {
    int err;
    if (flag[12]) {
      docfs=stdout;
      err=0;
    } else {
      docfs=stderr;
      err=INSUFFICIENT_COMMAND_ARGS;
    }
    fprintf(docfs, "\n");
    fprintf(docfs, "usage: pc_proxy_predict [-0 i0|-i t0] [-N N|-f tf] [-A ncv] [-v nev] \n");
    fprintf(docfs, "                  matfile dates measurements lead window outfile\n");
    fprintf(docfs, "\n");
    fprintf(docfs, "where:\n");
    fprintf(docfs, "  matfile      = binary file containing array of matrices representing tracer mapping\n");
    fprintf(docfs, "  dates        = ASCII file containing dates corresponding to each sparse matrix\n");
    fprintf(docfs, "  measurements = ASCII file containing measurements and locations\n");
    fprintf(docfs, "  lead         = lead time in days\n");
    fprintf(docfs, "  window       = measurement window in days\n");
    fprintf(docfs, "  outfile      = binary file containing interpolated tracer field\n");
    fprintf(docfs, "\n");
    ctraj_optargs(docfs, "d0NifAvl-+CK?");
    fprintf(docfs, "\n");
    return err;
  }

  matfile=argv[1];
  datefile=argv[2];
  measurement_file=argv[3];
  sscanf(argv[4], "%f", &lead);
  if (flag[11]==0) lead2=lead;
  sscanf(argv[5], "%f", &window);
  if (cflag==0) outfile=argv[6];

  //read in the array of sparse matrices:
  //if (cflag!=1) {
    fprintf(docfs, "Reading file: %s\n", matfile);
    fs=fopen(matfile, "r");

//*FLAG* -- fixed-length data structure
    matall=new sparse_matrix[nall];

    for (int32_t i=0; i<nall; i++) {
      //printf("%d\n", i);
      ncon=matall[i].read(fs);
      if (ncon==0) {
        nall=i;
        break;
      }
    }
    fprintf(docfs, "%d sparse matrices read in\n", nall);

    //determine the dimensions of the matrics:
    matall[0].dimensions(m, n);
    assert(m==n);			//only square matrices need apply...

    fclose(fs);
  //}

  //read in the dates:
  t=new time_class[nall+1];

  //get time grids:
  fs=fopen(datefile, "r");
  //line=fget_line(fs);		//throw away first grid
  for (int32_t i=0; i<=nall; i++) {
    int32_t ind;
    line=fget_line(fs);
    sscanf(line, "%d %s", &ind, tstring);
    t[i].read_string(tstring);
    delete [] line;
  }
  fclose(fs);

  //if we are specifying dates, then we are interested in the final,
  //output fields:
  if (flag[5]) {
    date0.add(lead2-lead);
    i0=ceil(interpolate(t, nall+1, date0, -1));
  }

  if (N == -1 || N+i0>nall) {
    tf=t[nall-1];
    tf.add(-lead-window);
    N=bin_search_g(t, nall+1, tf, -1)-i0+1;
  }

  if (flag[6]) {
    datef.add(-lead2);
    N=bin_search_g(t, nall+1, datef, -1)-i0+1;
  }
  //***NOTE*** no range-checking...
  
  //read in measurements:
  fprintf(docfs, "Reading in measurements from, %s\n", measurement_file);
  samp=read_meas(measurement_file, &nsamp, dwid);
  //write_meas(samp, nsamp, stdout);

  //fprintf(stderr, "pc_proxy::main: nsamp=%d, nev=%d, ncv=%d\n", nsamp, nev, ncv);
  
  if (cflag == 0) {
    fs=fopen(outfile, "w");
    fwrite(&n, sizeof(n), 1, fs);
  }
  nfield=0;
  for (int32_t i=i0; i<N+i0; i++) {
    int32_t nall_local;		//to save memory...
    //calculate lead times:
    tf=t[i];

    tf.add(lead);
    l2=interpolate(t, nall+1, tf, -1);
    N2=(int32_t) (ceil(l2)-i);		//number of sparse matrix elements

    tf2=t[i];
    tf2.add(lead2);
    t1=tf2;
    t1.add(-window);
    t2=tf2;
    t2.add(window);
    //round to nearest index: (rel. location of measurement window)
    N3=(int32_t) (interpolate(t, nall+1, tf2, -1)-i+0.5);
    //to save memory while doing the interpolation:
    nall_local=bin_search(t, nall+1, t2, -1)+1-i;

    t1.write_string(tstring);
    fprintf(docfs, "Interpolating measurements between %s ", tstring);
    t2.write_string(tstring);
    fprintf(docfs, "and %s\n", tstring);

    //select measurements:
    samp_w=select_meas(t1, t2, samp, nsamp, &nsamp_w, hemi);
    //write_meas(samp_w, nsamp_w, stdout);
    t[N3+i].write_string(tstring);

    if (cflag) {
      //if the -C option is specified, we only output the number of samples per field:
      printf("%d %s: %d\n", nfield, tstring, nsamp_w);
      nfield++;
      nsamp_ave+=nsamp_w;
      delete [] samp_w;
      continue;
    }

    if (nsamp_w < nev) {
      fprintf(stderr, "pc_proxy_predict: %s not included in time series\n", tstring);
      fprintf(stderr, "            fewer data points (%ld) than singular vectors (%d)\n",
			nsamp_w, nev);
      err=OTHER_WARNING;
      continue;
    }

    printf("%d %s\n", nfield, tstring);
    nfield++;

    qvec=pc_proxy(matall+i, t+i, N2, nall_local, samp_w, nsamp_w, nev, ncv, kflag, N3);

    delete [] samp_w;

    //output final, interpolated field:
    //printf("%d: writing vector of length %d for %s\n", i, n, tstring);
    fwrite(qvec[0], sizeof(float), n, fs);

    delete [] qvec[0];
    delete [] qvec;
  }

  //fwrite(qvec[0], sizeof(float), n*(nall-i0), fs);
  if (cflag) {
    printf("Average: %f\n", nsamp_ave/nfield);
  } else {
    fclose(fs);
  }

  //if (cflag!=1) delete [] matall;
  delete [] t;
  delete [] samp;

  return err;
}

