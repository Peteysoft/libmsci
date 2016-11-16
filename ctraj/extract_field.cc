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
  int32_t ndim=0;
  int32_t vtype=0;

  size_t ncon;

  float *dum;

  int32_t *ngrid;

  int32_t index;

  //based on the file alone:
  int nvar;
  int nrec;
  int fsize;

  //calculated:
  int32_t nin;
  int32_t nout;

  float parm[2];

  float missing=0;

  ctraj_tfield_base<float> *tracer;

  void *optargs[20];
  int flags[20];

  float *qout;
  //the tracer vector (what we are after after all this work...)
  //long k:
  float *qvec;

  int32_t nlon=NLON;
  int32_t nlat=NLAT;
  float sl2=SIDELENGTH_Q;

  optargs[0]=&ndim;
  optargs[1]=&vtype;
  optargs[2]=&missing;
  argc=parse_command_opts(argc, argv, "dVa?", "%d%d%g%", 
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
      tracer=new ctraj_tfield_nd<float>(ndim);
      break;
    default:
      tracer=new ctraj_tfield_standard<float>();
      break;
  }

  if (typeid(*tracer)==typeid(ctraj_tfield_standard<float>)) {
    optargs[0]=&nlon;
    optargs[1]=&nlat;
    optargs[2]=&sl2;
    argc=parse_command_opts(argc, argv, "xyr", "%d%d", 
			optargs, flags, OPT_WHITESPACE);
    if (argc < 0) {
      fprintf(stderr, "Error parsing command line\n");
      exit(21);
    }
    nout=nlon*nlat;
    //this part is kind of ugly:
    parm[0]=nlon;
    parm[1]=nlat;
  } else {
    argc=tracer->setup(argc, argv);
    if (argc < 0) {
      fprintf(stderr, "Error parsing command line\n");
      exit(21);
    }
    ndim=tracer->ndim();
    ngrid=new int32_t[ndim];
    tracer->get_raw_grids(ngrid);
    nout=1;
    for (int i=0; i<ndim; i++) nout=nout*ngrid[i];
    delete [] ngrid;
    parm[0]=missing;
  }

  if (argc < 2 || flags[3]) {
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
    fprintf(docfs, "Syntax:  extract_field [-a] [-r rmax] [-n ngrid] [-x nlon] [-y nlat]\n");
    fprintf(docfs, "                             infile [index]\n");
    fprintf(docfs, "\n");
    fprintf(docfs, "Converts fields in azimuthal equidistant gridding, binary format to\n");
    fprintf(docfs, "lon-lat gridding, ASCII format (GMT compatible)\n");
    fprintf(docfs, "\n");
    fprintf(docfs, "where:\n");
    fprintf(docfs, "  infile  =  input file name, binary dump of an array of vectors\n");
    fprintf(docfs, "  index   =  0-based index of field to extract\n");
    fprintf(docfs, "             (- if argument is absent, returns number of records;\n");
    fprintf(docfs, "              - gridding is specified by -x and -y options.)\n\n");
    fprintf(docfs, "options:\n");
    ctraj_optargs(docfs, "anrxydV?", 1);
    if (flags[0]) {
      tracer->help(docfs);
    }
    fprintf(docfs, "\n");
    return err;
  }

  infile=argv[1];
  if (argc>2) sscanf(argv[2], "%d", &index); else index=-1;
  fs=fopen(infile, "r");

  fread(&nvar, sizeof(nvar), 1, fs);

  fseek(fs, 0, SEEK_END);
  fsize=ftell(fs);
  nrec=(fsize-sizeof(nvar))/nvar/sizeof(float);

  if (index<0) {
    printf("%d\n", nrec);
    fclose(fs);
    exit(0);
  }

  if (typeid(*tracer)==typeid(ctraj_tfield_standard<float>)) {
    //initialize the tracer object directly:
    ((ctraj_tfield_standard<float> *) tracer)->init2(nvar, sl2);
  }

  nin=tracer->nel();

  if (nin!=nvar) {
    fprintf(stderr, "extract_field: Something went wrong\n");
    fprintf(stderr, " calculated versus actual number of records don't match (%d vs %d)\n", nin, nvar);
    exit(-1);
  }

  qvec=new float[nin];

  //seek out the position:
  fseek(fs, index*nvar*sizeof(float)+sizeof(nvar), SEEK_SET);
  //read the field:
  fread(qvec, sizeof(float), nvar, fs);

  fclose(fs);


  qout=tracer->from(qvec, parm);
  //printf("%d\n", k);

  //print the converted file to stdout:
  for (int i=0; i<nout; i++) fprintf(stdout, "%g\n", qout[i]);

  delete [] qout;
  delete [] qvec;
  delete tracer;

}

