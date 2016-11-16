//
// Copywrite 2004 Peter Mills.  All rights reserved.
//
// Implements a simple two-dimensional trajectory model.
//

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "error_codes.h"
#include "parse_command_opts.h"
#include "rk_dumb_ts.h"
#include "read_ascii_all.h"
#include "ctraj_defaults.h"

#include "ctraj_vfield_base.h"

// syntax:
// model1 vfile initfile outfile
//
// where:
//  	vfile is the binary file containing the velocity field.  Stored in
//  		the format read and written by the "composite_dataset" object.
//
//	initfile is an ascii file containing a list of initial conditions.
//		The initial conditions for each trajectory are stored on
//		a separate line.  The fields are:
//
//		date                    x0 y0 dt nt
//
//		The date is a fixed width field and has the following format:
//		year/month/day-hour:minute:second
//
//		x0 and y0 are floating point fields giving the initial
//		position.
//
//		dt is the non-dimensional time step as a fraction of the
//		time grid spacing of the velocity field.
//
//		nt is the number of time steps.
//
//	outfile is the output file.  The time and position at each time
//		step is simply written to a separate line in an ascii file.
//		The trajectories are written in order, meaning there
//		will be nt+1 consecutive lines for each trajectory.
//
// history:	2004-2-19 PM: formally documented.
// 		- changed file format for velocity fields
// 		- entire velocity field is no longer read in
// 		  (see dependent_swap.cpp)
// 		2004-5-27 PM: updated documentation to reflect changes
//		2005 sometime PM: upgraded to model2...
//
// bugs:	Entire velocity field is read in.  This may cause problems
// 		once it gets too big.  A better solution would be to only
// 		read in those time grids that are needed or read them in
// 		as they are needed.  The likely final "solution" will be
// 		to code the data structures used to store them so that they
// 		can swap from the disk.
//		--> done
//		However, there is still the problem of longer trajectories:
//		for each trajectory the field will be swapped in and out.
//		If all the trajectories start at the same date, it would be
//		better to run them concurrently.
//
// 		No syntax or range checking.  GIGO.

using namespace libpetey;
using namespace ctraj;

int main(int argc, char *argv[]) {
  char **line;		//one line from input file
  char *fg_res;

  int32_t *nt;				//number of time steps
  long ninit;
  double *tstep;				//time step
  char date_str[TSTRING_LEN];		//a date as a string
  int32_t *tind;			//**number of time steps already taken
  double *tind0;
  double tind0min, tind0max;
  double toffs;
  int32_t lt;				//time index as longword

  float ***result;			//final conditions
  int32_t **hemi;
  
  int32_t checksteps;
  int32_t nt_coarse;		//number of coarse steps
  int32_t ss, srem;

  ctraj_vfield_base<float> *vfield;		//nicely encapsulates (most of) the process
  void *param[2];

  FILE *initfun;		//initialization file unit
  FILE *outfun;			//output file unit

  int32_t ndim;
  float *x0;

  //read from the command line:
  int32_t twid=TFIELD_WIDTH;
  int32_t nrk=TSTEP_NFINE;
  double ind1;
  double tstep1=1;

  //read from the initialization file:
  float checkint=1;
  int32_t nt1;

  //type of velocity field:
  int32_t vtype=0;

  int flag[20];
  void *optarg[20];

  //defaults:
  ind1=-1;
  nt_coarse=0;

  optarg[0]=&nrk;
  optarg[1]=&twid;
  optarg[2]=&checkint;

  argc=parse_command_opts(argc, argv, "kth?", "%d%d%g%d%", optarg, flag, OPT_WHITESPACE+2);
  if (argc < 0) {
    fprintf(stderr, "ctraj_model: error parsing command line\n");
    argc=-argc;
  }
  vfield=ctraj_loader(argc, argv, ind1, nt_coarse);
  if (argc < 0) {
    fprintf(stderr, "ctraj_model: error parsing command line\n");
    argc=-argc;
  }

  if (argc <= 3 || flag[3]) {
    FILE *docfs;
    int err;
    if (flag[7]) {
      docfs=stdout;
      err=0;
    } else {
      docfs=stderr;
      err=INSUFFICIENT_COMMAND_ARGS;
    }

    fprintf(docfs, "\n");
    fprintf(docfs, "Runs the trajectory model over one hemisphere for a series of initial conditions\n");
    fprintf(docfs, "\n");
    fprintf(docfs, "syntax:\n");
    fprintf(docfs, "	ctraj_model [v-field-arguments] infile outfile [checkint]\n");
    fprintf(docfs, "\n");
    fprintf(docfs, "where:\n");
    fprintf(docfs, "    v-field-arguments are arguments supplied to the velocity field loader\n");
    fprintf(docfs, "            the -V option controls the type of velocity field\n");
    fprintf(docfs, "    infile  Ascii file containing initial conditions\n");
    fprintf(docfs, "            The format of each one-column record is as follows:\n");
    fprintf(docfs, "\n");
    fprintf(docfs, "             date                    x0 y0 dt nt\n\n");
    fprintf(docfs, "             date is a fixed width (%d) field and has the following format:\n", TFIELD_WIDTH);
    fprintf(docfs, "                  year/month/day-hour:minute:second\n");
    fprintf(docfs, "\n");
    fprintf(docfs, "             x0   floating point fields giving the initial position\n");
    fprintf(docfs, "             y0   \n");
    fprintf(docfs, "\n");
    fprintf(docfs, "             dt   is the non-dimensional time step as a fraction of the\n");
    fprintf(docfs, "                  time grid spacing of the velocity field.\n");
    fprintf(docfs, "\n");
    fprintf(docfs, "             nt   is the number of time steps.");
    fprintf(docfs, "\n");
    fprintf(docfs, "    outfile  Ascii file containing final, integrated trajectories.\n");
    fprintf(docfs, "             The time and position at each time\n");
    fprintf(docfs, "             step is simply written to a separate line in an ascii file.\n");
    fprintf(docfs, "             The trajectories are written in order, meaning there\n");
    fprintf(docfs, "             will be nt+1 consecutive lines for each trajectory.\n");
    fprintf(docfs, "\n");
    fprintf(docfs, "    checkint Interval to check for switching hemispheres\n");
    fprintf(docfs, "\n");
    fprintf(docfs, "options:\n");
    ctraj_optargs(docfs, "Vhif0Ntk", 1);
    if (vfield!=NULL) vfield->help(docfs);
    return err;
  }

  //read the initialization file:
  line=read_ascii_all(argv[1], &ninit);

  nt=new int32_t[ninit];
  tind0=new double[ninit];
  tind=new int32_t[ninit];
  tstep=new double[ninit];
  hemi=new int32_t*[ninit];

  result=new float **[ninit];
  ndim=vfield->ndim();

  //collect initial conditions:
  tind0max=0;
  x0=new float[ndim];
  for (long i=0; i<ninit; i++) {
    int strptr=twid+1;
    int nchar;
    for (int32_t j=0; j<ndim; j++) {
      sscanf(line[i]+strptr, "%g%n", x0+j, &nchar);
      strptr+=nchar;
    }
    sscanf(line[i]+strptr, "%lg %d", &tstep1, &nt1);
    if (flag[2]==0) tstep[i]=tstep1; else tstep[i]=checkint/nrk;
    if (ind1==-1) {
      line[i][twid]='\0';
      tind0[i]=vfield->get_tind(line[i]);
    } else {
      tind0[i]=ind1;
    }
    if (nt_coarse==0) {
      nt[i]=nt1;
    } else {
      nt[i]=nt_coarse*nrk;
      printf("nt[%d]=%d\n", i, nt[i]);
    }
    result[i]=new float *[nt[i]+1];
    result[i][0]=new float[(nt[i]+1)*ndim];
    for (int32_t j=1; j<=nt[i]; j++) result[i][j]=result[i][0]+ndim*j;

    hemi[i]=new int32_t[nt[i]+1];

    for (int32_t j=0; j<ndim; j++) result[i][0][j]=x0[j];

    hemi[i][0]=vfield->absolute(-1, result[i][0]);

    if (tind0[i] < tind0min || i==0) tind0min=tind0[i];
    if (tind0[i]+tstep[i]*nt[i] > tind0max) tind0max=tind0[i]+tstep[i]*nt[i];
    tind[i]=0;
    delete [] line[i];
  }
  delete [] line;
  delete [] x0;

  //perform the integration:
  param[0]=vfield;		//first parameter if velocity field
  for (double i=tind0min; i<=tind0max; ) {
    i+=checkint;
    printf("%lg:", i);
    for (long j=0; j<ninit; j++) {
      //vfield->get_t(tind[j]*tstep[j]+tind0[j], date_str);
      //printf(" %s", date_str);

      //the algorithm is designed so that if there is any time overlap between
      //trajectories, then they are integrated simultaneously:
      if (i>tind[j]*tstep[j]+tind0[j] && tind[j]<nt[j]) {
        nt1=ceil((i-tind0[j])/tstep[j])-tind[j];
        if (nt1+tind[j]>nt[j]) nt1=nt[j]-tind[j];

        //if "hemi" (domain index) is less than 0, then point has gone out-
        //-of-bounds: stop integration
        if (hemi[j][tind[j]] < 0) {
	  for (int k=1; k<=nt1; k++) {
            for (int dim=0; dim<ndim; dim++) result[j][tind[j]+k][dim]=result[j][tind[j]][dim];
            hemi[j][tind[j]+k]=hemi[j][tind[j]];
          }
          continue;
        }

        //do the integration:
        param[1]=hemi[j]+tind[j];		//second parameter is the domain
        //Runge-Kutta integration, fixed time-step, default derivative
        //calculator:
        rk_dumb(tind0[j]+tstep[j]*tind[j], result[j][tind[j]], (long) ndim, 
			tstep[j], (long) nt1, result[j]+tind[j], 
			(void *) param, &ctraj_deriv);
        //all fine time-steps within single coarse step have same domain:
	for (int k=1; k<=nt1; k++) hemi[j][tind[j]+k]=hemi[j][tind[j]];
        //check if trajectory has strayed from domain:
        hemi[j][tind[j]+nt1]=vfield->fix(hemi[j][tind[j]+nt1], 
			tind0[j]+tstep[j]*tind[j], result[j][tind[j]+nt1]);
	tind[j]+=nt1;
        //printf("x");
      }
    }
    printf("\n");
  }

  outfun=fopen(argv[2], "w");

  for (long i=0; i<ninit; i++) {
    //write it to the output file:
    for (long j=0; j<=nt[i]; j++) {
      //compute the date and transform to a pretty ascii format:
      ind1=tind0[i]+j*tstep[i];
      vfield->get_t(ind1, date_str);
      //convert to global coordinates:
      vfield->absolute(hemi[i][j], result[i][j]);
      //print to output file:
      fprintf(outfun, "%23s", date_str);
      for (int32_t k=0; k<ndim; k++) fprintf(outfun, " %14.7g", result[i][j][k]);
      fprintf(outfun, "\n");
    }
    delete [] result[i][0];
    delete [] result[i];
    delete [] hemi[i];
  }
  //close the open file pointers:
  fclose(outfun);

  delete [] result;
  delete [] hemi;
  delete [] tind0;
  delete [] tind;
  delete [] tstep;
  delete [] nt;

  delete vfield;
}

