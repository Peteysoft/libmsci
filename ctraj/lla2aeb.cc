//lon-lat ascii to azimuthal-equidistant binary...
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <typeinfo>

#include "error_codes.h"
#include "parse_command_opts.h"

#include "ctraj_defaults.h"

#include "ctraj_tfield_nd.h"
#include "ctraj_tfield_standard.h"

using namespace libpetey;
using namespace ctraj;

int main(int argc, char **argv) {

  char *infile;
  char *outfile;

  FILE *fs;
  int32_t nvar=0;
  int32_t vtype=0;

  size_t ncon;

  float *dum;

  int32_t *ngrid;

  int32_t nin;
  int32_t nout;

  ctraj_tfield_base<float> *tracer;

  void *optargs[20];
  int flags[20];

  float *qin;
  //the tracer vector (what we are after after all this work...)
  //long k:
  float *qvec;

  int32_t nlon=NLON;
  int32_t nlat=NLAT;
  float parm[2];

  optargs[0]=&nvar;
  optargs[1]=&vtype;
  argc=parse_command_opts(argc, argv, "dV?", "%d%d%", 
		  optargs, flags, OPT_WHITESPACE+2);
  if (argc < 0) {
    fprintf(stderr, "Error parsing command line\n");
    exit(21);
  }
  //printf("%d %g %d %d %s\n", ngrid, rmax, nlon, nlat, mapfile);

  switch (vtype) {
    case (0):
      tracer=new ctraj_tfield_standard<float>();
      break;
    case (1):
      tracer=new ctraj_tfield_nd<float>(2);
      break;
    case (2):
      tracer=new ctraj_tfield_nd<float>(nvar);
      break;
    default:
      tracer=new ctraj_tfield_standard<float>();
      break;
  }

  argc=tracer->setup(argc, argv);
  if (argc < 0) {
    fprintf(stderr, "Error parsing command line\n");
    exit(21);
  }

  if (typeid(*tracer)==typeid(ctraj_tfield_standard<float>)) {
    optargs[0]=&nlon;
    optargs[1]=&nlat;
    argc=parse_command_opts(argc, argv, "xy", "%d%d", 
			optargs, flags, OPT_WHITESPACE);
    if (argc < 0) {
      fprintf(stderr, "Error parsing command line\n");
      exit(21);
    }
    nin=nlon*nlat;
    parm[0]=nlon;
    parm[1]=nlat;
  } else {
    nvar=tracer->ndim();
    ngrid=new int32_t[nvar];
    tracer->get_raw_grids(ngrid);
    nin=1;
    for (int i=0; i<nvar; i++) nin=nin*ngrid[i];
    delete [] ngrid;
  }

  if (argc < 2 || flags[2]) {
    FILE *docfs;
    int err;
    if (flags[6]) {
      docfs=stdout;
      err=0;
    } else {
      docfs=stderr;
      err=INSUFFICIENT_COMMAND_ARGS;
    }
    fprintf(docfs, "\n");
    fprintf(docfs, "Syntax:  lla2aeb [-a] [-r rmax] [-n ngrid] [-x nlon] [-y nlat]\n");
    fprintf(docfs, "                             [infile] outfile\n");
    fprintf(docfs, "\n");
    fprintf(docfs, "options:\n");
    ctraj_optargs(docfs, "nrxydV?", 1);
    if (flags[1]) {
      tracer->help(docfs);
    }
    fprintf(docfs, "\n");
    return err;
  } else if (argc == 2) {
    infile=NULL;
    outfile=argv[1];
  } else if (argc > 2) {
    infile=argv[1];
    outfile=argv[2];
  }

  qin=new float[nin];

  if (infile != NULL) fs=fopen(infile, "r"); else fs=stdin;
  for (int i=0; i<nin; i++) {
    fscanf(fs, "%g", qin+i);
  }
  if (infile != NULL) fclose(fs);

  nout=tracer->nel();
  qvec=tracer->to(qin, parm);
  //printf("%d\n", k);

  fs=fopen(outfile, "w");
  //write size of each record for easier parsing and compatibility with libagf:
  fwrite(&nout, sizeof(nout), 1, fs);
  fwrite(qvec, sizeof(float), nout, fs);

  delete [] qin;
  delete [] qvec;
  delete tracer;

}
