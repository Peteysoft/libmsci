//
// Returns border samples for the pair of synthetic test classes.
//

#include <stdio.h>
#include <math.h>
#include <string.h>

#include "error_codes.h"
#include "supernewton.h"

#include "agf_lib.h"

#include "sample_class1_obj.h"
#include "sample_class2_obj.h"

#define NSAMP_DEFAULT 100
#define TOL_DEFAULT 0.0001

#define NVAR 2

#define DEFAULT_PRATIO 2

namespace libagf {
  struct cbsc_params {
    sample_class1_obj<real_a> *sc1;
    sample_class2_obj<real_a> *sc2;
    real_a pratio;
  };

  int scsample(void *params1, real_a *x1, real_a *x2) {
    cbsc_params *params=(cbsc_params *) params1;
    params->sc1->sample(x1[0], x1[1]);
    params->sc2->sample(x2[0], x2[1]);
    return 0;
  }

  real_a scrfunc(real_a *x, void *params1, real_a *drdx) {
    cbsc_params *params;	//objects for estimating class prob.
    real_a p1, p2;		//posterior probabilities of each class
    real_a dp1dx, dp1dy;	//derivatives of the prob. of 1st class
    real_a dp2dx, dp2dy;	//  "			2nd  "
    real_a p;			//total probability
    real_a r;			//difference in conditional prob.

    params=(cbsc_params *) params1;

    p1=params->sc1->pdf(x[0], x[1], dp1dx, dp1dy);
    p2=params->sc2->pdf(x[0], x[1], dp2dx, dp2dy);
    p=p1+params->pratio*p2;

    //derivation is not hard, but it's not trivial either
    //--I've done it at least twice!
    drdx[0]=2*params->pratio*(p1*dp2dx-p2*dp1dx)/p/p;
    drdx[1]=2*params->pratio*(p1*dp2dy-p2*dp1dy)/p/p;
  
    r=(params->pratio*p2-p1)/p;

    //printf("%f %f %f\n", *y, params->drdx, params->drdy);
    return r;
  }

}

using namespace libagf;

int main(int argc, char *argv[]) {
  char *brdfile;	//binary file sampling border
  char *grdfile;	//binary file containing gradient vectors

  agf_command_opts opt_args;

  real_a **border;	//border vectors
  real_a **gradient;	//gradient vectors

  FILE *fs;		//file stream

  sample_class1_obj<real_a> *sc1;
  sample_class2_obj<real_a> *sc2;

  cbsc_params params;	//parameters to pass to the minimization function

  dim_ta nvar=NVAR;

  int exit_code;

  //set defaults and parse command line options:
  opt_args.n=NSAMP_DEFAULT;
  opt_args.tol=TOL_DEFAULT;
  opt_args.pratio=DEFAULT_PRATIO;

  exit_code=agf_parse_command_opts(argc, argv, "r:s:t:i:N:X:", &opt_args);
  if (exit_code==FATAL_COMMAND_OPTION_PARSE_ERROR) return exit_code;

  if (argc < 1) {
    printf("\n");
    printf("syntax:  sc_borders [-s n] [-t tol] [-r R] [-X ratio] border\n");
    printf("\n");
    printf("arguments:\n");
    printf("  border   = base name of output files:\n");
    printf("               .brd samples the border;\n");
    printf("               .bgd contains gradient vectors\n");
    printf("\n");
    printf("options:\n");
    printf("  -X ratio = ratio between class size (P(2)/P(1)) (default=%d)\n", DEFAULT_PRATIO); 
    printf("  -s n     = number of times to sample the border (default=%d)\n", NSAMP_DEFAULT);
    printf("  -t tol   = desired tolerance (default=%g)\n", TOL_DEFAULT);
    printf("  -r R     = location of discrimination border (default=0)\n");
    printf("  -i/-N maxiter = maximum number of iterations in supernewton\n");
    printf("             (default=%d)\n", BORDERS_MAXITER);
    printf("\n");
    return INSUFFICIENT_COMMAND_ARGS;
  }

  brdfile=new char[strlen(argv[0])+5];
  strcpy(brdfile, argv[0]);
  strcat(brdfile, ".brd");

  grdfile=new char[strlen(argv[0])+5];
  strcpy(grdfile, argv[0]);
  strcat(grdfile, ".bgd");

  //allocate the arrays for holding the results:
  border=new real_a *[opt_args.n];
  border[0]=new real_a[opt_args.n*nvar];
  for (nel_ta i=1; i<opt_args.n; i++) border[i]=border[0]+i*nvar;

  gradient=new real_a *[opt_args.n];
  gradient[0]=new real_a[opt_args.n*nvar];
  for (nel_ta i=1; i<opt_args.n; i++) gradient[i]=gradient[0]+i*nvar;

  sc1=new sample_class1_obj<real_a>();
  sc2=new sample_class2_obj<real_a>();

  params.sc1=sc1;
  params.sc2=sc2;
  params.pratio=DEFAULT_PRATIO;

  opt_args.n=sample_class_borders(&scrfunc, &scsample, (void *) &params, 
		opt_args.n, nvar, opt_args.tol, agf_global_borders_maxiter,
		border, gradient, opt_args.rthresh);

  //write them to a file:
  fs=fopen(brdfile, "w");
  fwrite(&nvar, sizeof(real_a), 1, fs);
  fwrite(border[0], sizeof(real_a), nvar*opt_args.n, fs);
  fclose(fs);

  fs=fopen(grdfile, "w");
  fwrite(&nvar, sizeof(real_a), 1, fs);
  fwrite(gradient[0], sizeof(real_a), nvar*opt_args.n, fs);
  fclose(fs);

  //clean up:

  //delete character strings containing file names:
  delete[] brdfile;
  delete[] grdfile;

  //delete integer and real_aing point arrays:
  delete [] border[0];
  delete [] border;
  delete [] gradient[0];
  delete [] gradient;

  delete sc1;
  delete sc2;

  return exit_code;

}


