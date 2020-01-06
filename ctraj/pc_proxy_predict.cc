#include <vector>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include "sparse_element.h"

#include "error_codes.h"
#include "parse_command_opts.h"
#include "time_class.h"
#include "peteys_tmpl_lib.h"
#include "read_ascii_all.h"

#include "tracer_anal.h"
#include "ctraj_defaults.h"
#include "av.h"

using namespace std;
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
  size_t ncon;

  sparse_matrix *matall;
  int32_t nall;

  ind_t m, n;	//size of matrix  ****integer types should be fixed...  ****
  int32_t i0, N;		//start point, number of matrices
  int32_t nev, ncv;		//number of eigenvectors, number of Arnoldi vectors

  int32_t nvar;			//for compatibility with libagf

  time_class *t;		//time grids

  char tstring[TSTRING_LEN];
  char *line;

  //for reading and storing the samples:
  long nsamp;
  meas_data *samp;

  //interpolated fields:
  float **qvec;

  //for calculating lead times:
  float int_time;		//integration time in days
  float lead;		//days between start of transport map and measurement window
  float window;			//measurement window in days (+/-)

  //measurement data within window:
  meas_data *samp_w;
  long nsamp_w;

  //width of date field in measurement data:
  int32_t dwid=TFIELD_WIDTH;

  //date range:
  time_class date0;
  time_class datef;

  int flag[20];
  void *optarg[20];

  //northern or southern hemisphere:
  int hemi=0;

  int err=0;

  int cflag=0;	//-Q (query) option: just count measurments in each window
  int kflag=0;		//-K option: constant term included in fit

  //for getting sample statistics:
  float nsamp_ave=0;
  int nsamp_min;
  int nsamp_max=0;
  int32_t nfield=0;

  i0=-1;
  N=-1;

  ncv=NARNOLDI;
  nev=NEIG;

  optarg[0]=&dwid;
  optarg[3]=&i0;
  optarg[4]=&N;
  optarg[7]=&ncv;
  optarg[8]=&nev;
  optarg[11]=&lead;
  argc=parse_command_opts(argc, argv, "t-+0NifAvQKl?", 
		"%d%%%d%d%s%s%d%d%%%g%", optarg, flag, OPT_WHITESPACE);
  if (argc<0) exit(21);
  if (flag[1]) hemi=-1; else if (flag[2]) hemi=1;
  if (flag[5]) date0.read_string((char *) optarg[5]);
  if (flag[6]) datef.read_string((char *) optarg[6]);
  cflag=flag[9];
  kflag=flag[10];

  if ((argc<6 || (cflag==0 && argc < 7)) || flag[12]) {
    if (flag[12]) {
      docfs=stdout;
      err=0;
    } else {
      docfs=stderr;
      err=INSUFFICIENT_COMMAND_ARGS;
    }
    fprintf(docfs, "\n");
    fprintf(docfs, "usage: pc_proxy_predict [-0 i0|-i t0] [-N N|-f tf] [-A ncv] [-v nev] \n");
    fprintf(docfs, "                  matfile dates measurements int_t window outfile\n");
    fprintf(docfs, "\n");
    fprintf(docfs, "where:\n");
    fprintf(docfs, "  matfile      = binary array of sparse matrices representing transport map\n");
    fprintf(docfs, "  dates        = ASCII file of time grids corresponding to above\n");
    fprintf(docfs, "  measurements = ASCII file containing measurements and locations\n");
    fprintf(docfs, "  int_t        = tracer integration time in days\n");
    fprintf(docfs, "  window       = measurement window in days (+/-)\n");
    fprintf(docfs, "  outfile      = binary file of interpolated tracer fields\n");
    fprintf(docfs, "\n");
    ctraj_optargs(docfs, "t0NifAvl-+K?");
    fprintf(docfs, "  -Q   count the number of measurements per interpolate\n");
    fprintf(docfs, "\n");
    return err;
  }

  matfile=argv[1];
  datefile=argv[2];
  measurement_file=argv[3];
  sscanf(argv[4], "%f", &int_time);
  if (flag[11]==0) lead=int_time;
  sscanf(argv[5], "%f", &window);
  if (cflag==0) outfile=argv[6];

  //get time grids:
  fs=fopen(datefile, "r");
  //line=fget_line(fs);		//throw away first grid
  vector<time_class> alldates(0);
  while (feof(fs)==0) {
    int32_t ind;
    int ncon;
    time_class thist;
    line=fget_line(fs);
    if (line==NULL) break;
    ncon=sscanf(line, "%d %s", &ind, tstring);
    if (ncon!=2) break;
    thist.read_string(tstring);
    alldates.push_back(thist);
    delete [] line;
  }
  nall=alldates.size()-1;
  t=alldates.data();
  fclose(fs);

  //if we are specifying dates, then we are interested in the final,
  //output fields:
  if (flag[5]) {
    i0=ceil(interpolate(t, nall+1, date0, -1));
    if (i0<0) {
      fprintf(stderr, "Insufficient date coverage in tracer mapping\n");
      exit(SAMPLE_COUNT_MISMATCH);
    }
  }

  if (i0<0) {
    date0=t[0];
    if (lead-window<0) {
      date0.add(window);
    } else {
      date0.add(lead);
      i0=ceil(interpolate(t, nall+1, date0, -1));
    }
  }

  if (flag[6]) {
    N=bin_search_g(t, nall+1, datef, -1)-i0+1;
  }

  if (N < 0 || N+i0>nall) {
    datef=t[nall];
    if (lead+window>int_time) datef.add(window); else datef.add(-lead+int_time);
    datef.add(-window);
    N=bin_search_g(t, nall+1, datef, -1)-i0+1;
  }

  if (N<0) {
    fprintf(stderr, "Insufficient date coverage in tracer mapping\n");
    exit(SAMPLE_COUNT_MISMATCH);
  }

  //now that we've calculated all these dates and indices, it's time to read
  //in the sparse matrices:
  if (cflag!=1) {
    time_class t0;		//first time grid to load in
    time_class tf;		//last date we need
    int start;
    //we need data starting at this date:
    t0=t[i0];
    if (window>lead) t0.add(-window); else t0.add(-lead);
    start=bin_search_g(t, alldates.size(), t0, -1);
    if (start<0) {
      fprintf(stderr, "Insufficient date coverage in tracer mapping\n");
      exit(SAMPLE_COUNT_MISMATCH);
    }
    //printf("Starting index: %d\n", start);
    fprintf(docfs, "Reading file: %s\n", matfile);
    fs=fopen(matfile, "r");
    //scan ahead to the data that we need:
    for (int i=0; i<start; i++) {
      long nel;				//integer types should be fixed...
      ncon=fread(&m, sizeof(m), 1, fs);
      ncon+=fread(&n, sizeof(n), 1, fs);
      ncon+=fread(&nel, sizeof(nel), 1, fs);
      err=fseek(fs, nel*sizeof(sparse_el<int32_t, float>), SEEK_CUR);
      if (ncon!=3 || err!=0) {
        fprintf(stderr, "pc_proxy_predict: error reading file, %s\n", matfile);
        exit(FILE_READ_ERROR);
      }
    }
    if (m!=n) {
      fprintf(stderr, "pc_proxy_predict: dimension mismatch in sparse matrices\n");
      exit(DIMENSION_MISMATCH);
    }
    //and ending at this date:
    tf=t[i0+N];
    if (lead+window > int_time) tf.add(window); else tf.add(int_time-lead);
    nall=ceil(interpolate(t, alldates.size(), tf))-start+1;
    matall=new sparse_matrix[nall];
    for (int32_t i=0; i<nall; i++) {
      int32_t m1, n1;
      //printf("%d\n", i);
      ncon=matall[i].read(fs);
      if (ncon==0) {
        fprintf(stderr, "pc_proxy_predict: insufficient data found in sparse matrix file, %s\n", matfile);
        exit(SAMPLE_COUNT_MISMATCH);
      }
      matall[i].dimensions(m1, n1);
      if (m1!=m || n1!=n) {
        fprintf(stderr, "pc_proxy_predict: dimension mismatch in sparse matrices\n");
        exit(DIMENSION_MISMATCH);
      }
    }
    fclose(fs);
    fprintf(docfs, "%d sparse matrices read in\n", nall);
    //correct date array and indices:
    i0-=start;
    t+=start;
  }
  
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
  nsamp_min=nsamp;
  long last_ind1=-1;		//speed up binary searches
  long last_ind2=-1;		//speed up binary searches
  for (int32_t i=i0; i<N+i0; i++) {
    int32_t nall_local;		//to save memory...
    int32_t start;		//start grid to pass to pc_proxy
    time_class tgrid;		//current time grid to output
    time_class t1, t2;		//start and end of measurement window
    time_class t0, tf;		//start and end of integration time
    int32_t l1, l2;		//  "

    //current time grid:
    tgrid=t[i];

    //start and end of measurement window:
    t1=tgrid;
    t1.add(-window);
    t2=tgrid;
    t2.add(window);

    //start and end of integration period:
    t0=tgrid;
    t0.add(-lead);
    tf=tgrid;
    tf.add(int_time-lead);
    l1=ceil(interpolate(t, nall+1, t0, last_ind1));
    l2=bin_search_g(t, nall+1, tf, last_ind2);

    //informational message (interpolation interval):
    t1.write_string(tstring);
    fprintf(docfs, "Interpolating measurements between %s ", tstring);
    t2.write_string(tstring);
    fprintf(docfs, "and %s\n", tstring);

    //select measurements:
    samp_w=select_meas(t1, t2, samp, nsamp, &nsamp_w, hemi);
    //write_meas(samp_w, nsamp_w, stdout);

    tgrid.write_string(tstring);
    if (cflag) {
      //if the -C option is specified, we only output the number of samples per field:
      printf("%d %s: %d\n", nfield, tstring, nsamp_w);
      nfield++;
      nsamp_ave+=nsamp_w;
      if (nsamp_w<nsamp_min) nsamp_min=nsamp_w; 
      		else if (nsamp_w>nsamp_max) nsamp_max=nsamp_w;
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

    //earliest time grid we need for interpolation:
    if (window>lead) {
      start=bin_search_g(t, nall+1, t1, last_ind1);
    } else {
      start=l1;
    }

    //to save memory while doing the interpolation:
    //(maximum number of grids we need to pass to pc_proxy)
    if (tf>t2) {
      nall_local=ceil(interpolate(t, nall+1, tf, last_ind2))-start+1;
    } else {
      nall_local=l2-start+1;
    }

    //printf("%d, %d, %d, %d\n", start, i-l1, nall_local, i-start);
    qvec=pc_proxy(matall+start, t+start, i-l1, nall_local, samp_w, nsamp_w, nev, ncv, kflag, i-start);

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
    printf("Min:     %d\n", nsamp_min);
    printf("Max:     %d\n", nsamp_max);
  } else {
    fclose(fs);
  }

  //if (cflag!=1) delete [] matall;
  //delete [] t;
  delete [] samp;

  return err;
}

