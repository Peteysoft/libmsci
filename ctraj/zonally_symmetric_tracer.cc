#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <math.h>

#include "parse_command_opts.h"
#include "ctraj_defaults.h"
#include "ctraj_tfield_nd.h"

using namespace libpetey;
using namespace ctraj;

int main(int argc, char **argv) {
  char c;
  long ncon;
  long nlon, nlat;
  float q;

  nlon=NLON;
  nlat=NLAT;

  int pflag=0;

  int32_t vtype=0;
  ctraj_tfield_nd<float> *tracer;

  float missing=0;

  void * optarg[10];
  int flag[10];

  optarg[0]=&nlon;
  optarg[1]=&nlat;
  optarg[4]=&vtype;
  optarg[5]=&missing;


  argc=parse_command_opts(argc, argv, "xyE?Va", "%d%d%%%d%g", optarg, flag, OPT_WHITESPACE);
  if (argc < 0) {
             fprintf(stderr, "zonally_symmetric_tracer: Error parsing command line\n");
             exit(2);
  }
  pflag=flag[2];

  if (vtype==1) {
    tracer=new ctraj_tfield_nd<float>(2);
    tracer->setup(argc, argv);
  } else if (vtype==2) {
    tracer=new ctraj_tfield_nd<float>();
    tracer->setup(argc, argv);
  }

  if (flag[3]) { 
	    fprintf(stdout, "\n");
            fprintf(stdout, "Syntax:  zonally_symmetric tracer [-x nlon] [-y nlat] [-?] [-E]\n");
            fprintf(stdout, "\n");
            fprintf(stdout, "Generates a zonally-symmetric tracer field and prints it to stdout.\n");
            fprintf(stdout, "\n");
            fprintf(stdout, "options:\n");
	    ctraj_optargs(stdout, "xyEVa?");
	    fprintf(stdout, "\n");
	    exit(0);
  }

  if (vtype==1 || vtype==2) {
    int32_t ndim=tracer->ndim();
    int32_t nel=tracer->nel();
    float x[ndim];
    float r2;
    float qvec[nel];
    float *qout;
    int32_t nout;
    float parm[2];

    for (int i=0; i<nel; i++) {
      tracer->get_loc(i, x);
      r2=0;
      for (int j=0; j<ndim; j++) {
        r2+=x[j]*x[j];
      }
      qvec[i]=sqrt(r2);
    }
    //bit stupid:
    parm[0]=missing;
    qout=tracer->from(qvec, parm);
    nout=parm[1];
    for (int i=0; i<nout; i++) printf("%g\n", qout[i]);
    delete tracer;
  } else {
    for (long i=0; i<nlat; i++) {
      for (long j=0; j<nlon; j++) {
        if (pflag) {
          if (j <= nlon/2.) q=2.*j/nlon; else q=2.*(j-nlon)/nlon;
        } else {
          q=2.*i/(nlat-1)-1;
        }
        printf("%g\n", q);
      }
    }
  }

  return 0;

}

