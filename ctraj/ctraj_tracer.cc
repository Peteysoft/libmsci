//
// Semi-Lagrangian tracer integration. Outputs an array of sparse matrices
// representing the tracer dynamics.
//

#include <assert.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
//#include <getopt.h>
#include <stdlib.h>

#include "error_codes.h"

#include "rk_dumb_ts.h"
#include "parse_command_opts.h"

#include "sparse.h"

#include "ctraj_vfield_standard.h"
#include "ctraj_vfield_2d.h"
#include "ctraj_vfield_anal.h"

#include "ctraj_tfield_diffusion.h"
#include "ctraj_tfield_nd.h"

using namespace libpetey;
using namespace libsparse;
using namespace ctraj;

int main(int argc, char *argv[]) {
  //main data file names:
  char *outfile;

  //basic data structures/engine:
  ctraj_vfield_base<float> *vfield;		//velocity field
  ctraj_tfield_base<float> *tracer;		//tracer grid
  sparse<int32_t, float> map;			//transport map

  //composite dataset containing velocity field:
  //time step info:
  ind_type n;				//number of time steps
  int32_t nfine=TSTEP_NFINE;		//number of trajectory time steps
  double tstep_fine;			//time step for trajectory calc.
  float tstep=TSTEP_COARSE;		//time step for tracer field

  //general date variables:
  char date_str[TSTRING_LEN];		//a date as a string

  //time indices:
  double tind1;				//initial time index
  double *tind;				//vector of time indices
  int32_t lt;				//time index as integer

  int32_t ndim=2;			//number of dimensions
  int32_t nwt;				//current number of weights
  int32_t nwtmax;			//maximum number of weights
  int32_t nmap;				//number of points in tracer field

  //stuff for integration:
  float **result;			//integrated values

  void *param[2];		//pass to velocity field to R-K integrator

  az_eq_t<float> *metric;		//properties of metric space

  //interpolation indices:
  int32_t *ind;				//indices of interpolation points
  double *wt;				//weights for each point
  int32_t domain;			//domain of metric space

  FILE *outfun;				//output file unit

  int32_t vtype=0;			//type of velocity field/tracer
  int qflag=0;				//just print out the dates...
  int argc0;				//number of non-optional arguments

  int flag[20];				//argument present or not?
  void *optarg[20];			//returned optional arguments

  //parse the command line arguments:
  tind1=0;
  optarg[0]=&tstep;
  optarg[1]=&nfine;
  optarg[2]=&tind1;
  optarg[3]=&n;
  optarg[4]=&vtype;
  optarg[9]=&ndim;

  argc0=argc;

  argc=parse_command_opts(argc, argv, "hk0NVifQ?L", "%g%d%lg%d%d%s%s%%%d", optarg, flag, OPT_WHITESPACE+2);
  if (argc < 0) {
    fprintf(stderr, "ctraj_tracer: error parsing command line");
    argc=-argc;
  }
  qflag=flag[7];

  switch (vtype) {
    case (0):
      vfield=new ctraj_vfield_standard<float>();
      tracer=new ctraj_tfield_standard<float>();
      break;
    case (1):
      vfield=new ctraj_vfield_2d<float>();
      tracer=new ctraj_tfield_nd<float>(2);
      break;
    case (2):
      vfield=new ctraj_vfield_anal<float>(ndim);
      tracer=new ctraj_tfield_nd<float>(ndim);
      break;
    case (3):
      vfield=new ctraj_vfield_standard<float>();
      tracer=new ctraj_tfield_diffusion<float>();
      break;
    default:
      vfield=new ctraj_vfield_standard<float>();
      tracer=new ctraj_tfield_standard<float>();
      vtype=0;
      break;
  }

  if ((argc <= 2) || flag[8]) {
    FILE *docfs;
    int err;
    if (flag[3]) {
      docfs=stdout;
      err=0;
    } else {
      docfs=stderr;
      err=INSUFFICIENT_COMMAND_ARGS;
    }

    fprintf(docfs, "Purpose: two-dimensional, 'semi-Lagrangian' tracer simulation driven \n");
    fprintf(docfs, "by globally gridded wind fields.  For each time step, outputs a sparse \n");
    fprintf(docfs, "matrix defining the mapping from one tracer field to the next\n");
    fprintf(docfs, "\n");
    fprintf(docfs, "syntax:\n");
    fprintf(docfs, "ctraj_tracer [-Q] [-h dt] [-k nfine] [tracer arguments]\n");
    fprintf(docfs, "               [v-field arguments] outfile\n");
    fprintf(docfs, "\n");
    fprintf(docfs, " where:\n");
    fprintf(docfs, "   outfile   is the binary output file\n");
    fprintf(docfs, "\n");
    fprintf(docfs, "options:\n");
    ctraj_optargs(docfs, "hkif0NVQ?", 1);
    if (flag[4]) {
      fprintf(docfs, "\n");
      vfield->help(docfs);
      fprintf(docfs, "\n");
      tracer->help(docfs);
    } else {
      fprintf(docfs, "\nFor tracer and velocity field options, set -V flag\n");
    }
    fprintf(docfs, "\n");
    return err;
  }

  tracer->setup(argc, argv);
  nwtmax=tracer->nwt();

  vfield->setup(argc, argv);

  if (vtype==0 || vtype==3) {
    //somehow it's important to run this with different-sized Earth's:
    //(maybe this might someday be precise enough to worry about the non-
    //sphericity of the Earth!)
    metric=((ctraj_vfield_standard<float> *) vfield)->get_metric();
    ((ctraj_tfield_standard<float> *) tracer)->set_metric(metric);
  }

  //set integration time:
  if (flag[5]) tind1=vfield->get_tind((char *) optarg[5]);
  if (flag[6]) {
    double tindf=vfield->get_tind((char *) optarg[6]);
    if (tindf<tind1) {
      n=ceil(tind1-tindf);
      if (tstep>0) tstep=-tstep;
    } else {
      n=ceil(vfield->get_tind((char *) optarg[6])-tind1);
    }
    n/=abs(tstep);
  } else if (flag[3]==0) {
    n=(vfield->maxt()-tind1)/abs(tstep);
  }

  outfile=new char[strlen(argv[1])+1];
  strcpy(outfile, argv[1]);

  tind=new double[n+1];
  //figure out its location in relation to the velocity field:
  tind[0]=tind1;
  tstep_fine=-tstep/nfine;

  for (long i=1; i<=n; i++) tind[i]=tind[i-1]+tstep;

  nmap=tracer->nel();

  //initialize the sparse matrix:
  //map.extend(nmap*nwtmax);

  //initialize interpolation variables:
  wt=new double[nwtmax];
  ind=new int32_t[nwtmax];

  //initialize the vector of results for the Runge-Kutta integrations:
  result=new float * [nfine+1];
  result[0]=new float[ndim*(nfine+1)];
  for (int32_t i=1; i<=nfine; i++) result[i]=result[0]+i*ndim;

  //pass velocity field to R-K integrator:
  param[0]=vfield;

  //open the output file and write the headers:
  if (qflag==0) outfun=fopen(outfile, "w");

  //print out first time grid, but do not output a mapping for it:
  vfield->get_t(tind[0], date_str);
  printf("%d %s\n", 0, date_str);

  for (ind_type it=1; it<=n; it++) {
    //get date string and print it out:
    vfield->get_t(tind[it], date_str);
    printf("%d %s\n", it, date_str);

    if (qflag) continue;	//user is only querying the time grids

    //clear sparse matrix:
    map.reset(nmap, nmap);

    for (int32_t i=0; i<nmap; i++) {
      int d0i;
      float x0i[ndim];
      float tw;

      //set initial conditions:
      d0i=tracer->get_loc(i, x0i);
      //pass "domain" to R-K integrator:
      param[1]=&d0i;

      //do a Runge-Kutta integration:
      rk_dumb(tind[it], x0i, ndim, tstep_fine, nfine, result, (void *) param, &ctraj_deriv);

      //get interpolation coefficients:
      nwt=tracer->interpolate(d0i, result[nfine], ind, wt);

      //don't think I've ever used this feature (maybe it should be expunged?)
      tw=0;
      if (nwt < 0) {
        nwt=-nwt;
        //exclude missing indices and weights less than 0:
        for (int32_t j=0; j<nwt; j++) {
          if (ind<0 || wt[j]<0) {
            wt[j]=0;
          } else {
            tw+=wt[j];
          }
        }
        //re-normalize the weights:
        for (int32_t j=0; j<nwt; j++) wt[j]=wt[j]/tw;
      }

      //add interpolation coefficients to matrix:
      tw=0;
      for (int32_t j=0; j<nwt; j++) {
        map.add_el(wt[j], i, ind[j]);
        tw+=wt[j];
      }

      //why was this happening? (maybe it should be expunged too...)
      if (tw > 1.000001) {
        fprintf(stderr, "ctraj_tracer: warning: total of weights at row %d exceeds unity:\n", i);
        fprintf(stderr, "              %g", wt[0]);
        for (int32_t j=1; j<nwt; j++) fprintf(stderr, "+%g", wt[j]);
        fprintf(stderr, "=%g\n", tw);
      }

    }

    //write the date and the tracer field to a file:
    map.write(outfun);
    //map.print(stdout);

  }

  if (qflag==0) fclose(outfun);

  //clean up:
  delete [] outfile;
  delete [] tind;

  delete tracer;
  delete vfield;

  //delete [] x0[0];
  //delete [] x0;
  //delete [] d0;
  delete [] result[0];
  delete [] result;

  delete [] wt;
  delete [] ind;
}
