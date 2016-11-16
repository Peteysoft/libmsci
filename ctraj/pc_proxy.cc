#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include "error_codes.h"

#include "time_class.h"
#include "peteys_tmpl_lib.h"
#include "parse_command_opts.h"
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

  ind_t m, n;			//size of matrix
  int32_t nvar;			//make sure output file is libagf compatible
  int32_t i0, N;		//start point, end point, number of matrices
  int32_t nev, ncv;		//number of eigenvectors, number of Arnoldi vectors

  time_class *t;

  char tstring[TSTRING_LEN];
  char *line;

  //for reading and storing the samples:
  long nsamp;
  meas_data *samp;

  float *q0;
  float **qvec;

  int wflag=0;

  time_class t0;
  time_class tf;
  time_class ttest;

  short hemi=0;
  int32_t dwid=TFIELD_WIDTH;

  void *optarg[15];
  int flag[15];
  int cflag=0;

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
  argc=parse_command_opts(argc, argv, "tK?0NifAvW", "%d%%%d%d%s%s%d%d%", optarg, flag, OPT_WHITESPACE);
  cflag=flag[1];
  if (argc<0) exit(21);
  wflag=flag[9];
  if (flag[5]) t0.read_string((char *) optarg[5]);
  if (flag[6]) tf.read_string((char *) optarg[6]);

  if (argc<5 || flag[2]) {
    int err;
    if (flag[2]) {
      docfs=stdout;
      err=0;
    } else {
      docfs=stderr;
      err=INSUFFICIENT_COMMAND_ARGS;
    }
    fprintf(docfs, "\n");
    fprintf(docfs, "usage: pc_proxy [-A ncv] [-v nev] [-t dwid] [-W]\n");
    fprintf(docfs, "               [-i t0|-0 i0] [-f tf|-N n]\n");
    fprintf(docfs, "                  matfile dates measurements outfile\n");
    fprintf(docfs, "\n");
    fprintf(docfs, "where:\n");
    fprintf(docfs, "  matfile      = binary file containing array of matrices representing tracer mapping\n");
    fprintf(docfs, "  dates        = ASCII file containing dates corresponding to each sparse matrix\n");
    fprintf(docfs, "  measurements = ASCII file containing measurements and locations\n");
    fprintf(docfs, "  outfile      = binary file containing initial interpolated tracer field\n");
    fprintf(docfs, "\n");
    fprintf(docfs, "options:\n");
    ctraj_optargs(docfs, "t0NifAv-+WK");
    fprintf(docfs, "\n");
    return err;
  }

  matfile=argv[1];
  datefile=argv[2];
  measurement_file=argv[3];
  outfile=argv[4];

  //read in the array of sparse matrices:
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

  //read in the dates:
  t=new time_class[nall+1];

  fs=fopen(datefile, "r");
  //fgets(line, MAXLL, fs);		//no header

  //get time grids:
  //line=fget_line(fs);		//throw away first gird (thus both the map and
				//the resultant tracer are compatible)
  for (int32_t i=0; i<nall+1; i++) {
    int ind;
    line=fget_line(fs);
    sscanf(line, "%d %s", &ind, tstring);
    t[i].read_string(tstring);
    delete [] line;
  }
  fclose(fs);

  if (flag[5]) {
    i0=ceil(interpolate(t, nall+1, t0, -1));
  }

  if (N<0) N=nall-i0;

  if (flag[6]) {
    N=bin_search(t, nall+1, tf, -1)-i0;
  } 

  if (i0 < 0) i0=0;
  if (N > nall-i0) N=nall-i0;

  //read in measurements:
  fprintf(docfs, "Reading in measurements from, %s\n", measurement_file);
  samp=read_meas(measurement_file, &nsamp, dwid);

  if (nsamp < nev) {
    fprintf(docfs, "pc_proxy: error, too few data points (%ld vs. %d SV)\n", nsamp, nev);
    exit(411);
  }

  fprintf(docfs, "Performing interpolation...\n");
  
  fprintf(docfs, "pc_proxy::main: nsamp=%ld, nev=%d, ncv=%d, i0=%d, N=%d\n", nsamp, nev, ncv, i0, N);

  int index;
  if (wflag) index=-1; else index=N;
  qvec=pc_proxy(matall+i0, t+i0, N, nall-i0, samp, nsamp, nev, ncv, cflag, index);

  //printf("Generating final tracer...\n");

  //qvec=tracer_multiply(matall+i0, nall-i0, q0);

  //output final, interpolated initial field:
  fs=fopen(outfile, "w");
  nvar=n;
  if (wflag) {
    fprintf(docfs, "Writing %d vectors of length %d\n", nall-i0+1, n);
    fwrite(&nvar, sizeof(nvar), 1, fs);
    fwrite(qvec[0], sizeof(float), n*(nall-i0+1), fs);
  } else {
    fprintf(docfs, "Writing %d vectors of length %d\n", 1, n);
    fwrite(&nvar, sizeof(nvar), 1, fs);
    fwrite(qvec[0], sizeof(float), n, fs);
  }
  fclose(fs);

  //delete [] q0;

  delete [] qvec[0];
  delete [] qvec;

  delete [] matall;
  delete [] t;
  delete [] samp;
}

