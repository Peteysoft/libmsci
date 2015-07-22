#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>

#include "parse_command_opts.h"

#include "sparse.h"
#include "error_codes.h"

#include "ctraj_defaults.h"

using namespace libpetey;
using namespace libsparse;
using namespace ctraj;

int main(int argc, char *argv[]) {
  sparse_matrix mat;
  char *vect_file;
  char *mat_file;
  FILE *vect_fs;
  FILE *mat_fs;
  ind_t m, n;
  int32_t nvar;
  int32_t i0, i;
  int32_t N=-1;
  float *vect1, *vect2, *swp;
  size_t ncon;
  int errcode;
  void *optarg[3];
  int flags[3];

  float vmax, vmin;

  errcode=0;
  i0=0;

  optarg[0]=&i0;
  optarg[1]=&N;
  argc=parse_command_opts(argc, argv, "0N?", "%d%d%", optarg, flags, OPT_WHITESPACE);
  if (argc < 0) exit(FATAL_COMMAND_OPTION_PARSE_ERROR);
  
  if (argc < 3 || flags[2]) {
    FILE *docfs;
    int err;
    if (flags[2]) {
      docfs=stdout;
      err=0;
    } else {
      docfs=stderr;
      err=INSUFFICIENT_COMMAND_ARGS;
    }
    fprintf(docfs, "\n");
    fprintf(docfs, "  purpose: multiplies a vector with an array of sparse matrices,\n");
    fprintf(docfs, "           renormalizing each step\n");
    fprintf(docfs, "\n");
    fprintf(docfs, "  usage:  tracer_step2_w_renorm [-0 i0] [-N n] infile outfile\n\n");
    fprintf(docfs, "where:\n");
    fprintf(docfs, "  infile  = binary input file containing an array of\n");
    fprintf(docfs, "            sparse matrices (should be square...)\n");
    fprintf(docfs, "  outfile = binary output file containing initial vector\n");
    fprintf(docfs, "            and array of output vectors\n");
    fprintf(docfs, "\n");
    fprintf(docfs, "options:\n");
    ctraj_optargs(docfs, "0N?");
    fprintf(docfs, "\n");
    return err;
  }

  vect_file=argv[2];
  mat_file=argv[1];

  mat_fs=fopen(mat_file, "r");
  mat.read(mat_fs);
  mat.dimensions(m, n);
  //if (m>n) n=m;
  if (m != n) {
    fprintf(stderr, "tracer_step2_w_renorm: matrix must be square: %d x %d\n", m, n);
    exit(401);
  }

  //printf("%d %d\n", m, n);

  vect1=new float[n];
  vect2=new float[n];
  
  vect_fs=fopen(vect_file, "r+");
  fread(&nvar, sizeof(nvar), 1, vect_fs);
  if (nvar != n) {
    fprintf(stderr, "tracer_step2_w_renorm: matrix and vector dimensions do not match: %d vs %d\n", n, nvar);
    exit(401);
  }
  fread(vect2, sizeof(float), n, vect_fs);

  for (i=0; i<i0; i++) {
    if (mat.read(mat_fs)==0) {
      errcode=FILE_READ_ERROR;
      fprintf(stderr, "Error reading file: %s\n", mat_file);
      exit(errcode);
    }
  }

  mat.vect_mult(vect2, vect1);
  fwrite(vect1, sizeof(float), n, vect_fs);

  i=i0;

  while (!feof(mat_fs)) {
    printf("%d\n", i);
    if (flags[1] && i-i0 >= N) break;
    i++;
    if (mat.read(mat_fs)==0) break;
    mat.vect_mult(vect1, vect2);
    //normalize vect2:
    vmin=vect2[0];
    vmax=vect2[0];
    for (ind_t i=1; i<n; i++) {
      if (vect2[i] > vmax) vmax=vect2[i];
      if (vect2[i] < vmin) vmin=vect2[i];
    }
    for (ind_t i=0; i<n; i++) {
      vect2[i]=(vect2[i]-vmin)/(vmax-vmin)*2-1;
    }
    fwrite(vect2, sizeof(float), n, vect_fs);
    swp=vect1;
    vect1=vect2;
    vect2=swp;
  }

  fclose(mat_fs);
  fclose(vect_fs);
  
}

