#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <typeinfo>

#include "error_codes.h"
#include "parse_command_opts.h"
#include "ctraj_defaults.h"

#include "az_eq_t.h"
#include "ctraj_vfield_standard.h"
#include "ctraj_vfield_2d.h"
#include "ctraj_boundary_list.h"

//syntax: contour2 initfile
//
//where:
//	initfile:	name of the initialisation file
//
//format of the initialisation file:
//	nfile:	Binary file containing N hemisphere velocity field
//	sfile:	Binary file containing S hemisphere velocity field
//	infile:	Binary file containing initial contour
//	outfile:	Binary file to hold results
//	t0:	Date to start integration
//	tstep:	Time step (coarse)
//	n:	Number of steps
//	nrk:	Number of Runge-Kutta steps between each time step
//	thresh_arc:	Maximum radians of arc between each pair of points
//	min_spac:	Minimum spacing between each pair of points
//	max_spac:	Maximum spacing between each pair of points

#define MAXLL 200

using namespace libpetey;
using namespace ctraj;

int main(int argc, char *argv[]) {
  char *initfile;
  FILE *initfs;

  //the basic "engine":
  ctraj_boundary_list<double,float> *contour;		//for advecting the contour

  ctraj_vfield_base<float> *vfield;		//velocity field
  metric_base<float> *metric;			//metric

  char *infile;
  FILE *infs;
  char *outfile;
  FILE *outfs;
  float x, y, r;
  float tarc=MAXARC;		//maximum degrees of arc
  float minspac=DSMIN;		//minimum point spacing
  float maxspac=DSMAX;		//maximum point spacing
  double ts=TSTEP_COARSE;
  int32_t n;			//number of coarse time steps
  int32_t nrec;			//number of records in output file
  int32_t nrk=TSTEP_NFINE;	//number of fine time steps
  int32_t npts;
  int32_t nc;			//number of contours in input file
  int32_t magic=MAGIC;
  char *t0=NULL;
  char *tf=NULL;
  double tcur;
  char tstring[MAXLL];
  size_t nread;
  interpol_index ind1;
  int oflag;
  int dflag=0;			//use dates in '.bev' file
  int32_t maxnpt=-1;
  int32_t writeint=WRITE_INT;

  int flag[20];
  void *optarg[20];

  int argc0=argc;

  //defaults:
  ind1=-1;
  n=-1;

  optarg[0]=&tarc;
  optarg[1]=&minspac;
  optarg[2]=&maxspac;
  optarg[3]=&ts;
  optarg[4]=&nrk;
  optarg[5]=&maxnpt;
  optarg[6]=&writeint;

  argc=parse_command_opts(argc, argv, "cmshkMOo?if", "%g%g%g%lg%d%d%d%%%s%s", optarg, flag, OPT_WHITESPACE+2);
  if (argc < 0) {
    fprintf(stderr, "ctraj_contour: error parsing command line\n");
    argc=-argc;
  }
  oflag=flag[7];
  vfield=ctraj_loader(argc, argv, ind1, n);
  if (argc < 0) {
    fprintf(stderr, "ctraj_contour: error parsing command line\n");
    argc=-argc;
  }

  if (argc <= 3 || flag[8]) {
    FILE *docfs;
    int err;
    if (flag[8]) {
      docfs=stdout;
      err=0;
    } else {
      docfs=stderr;
      err=INSUFFICIENT_COMMAND_ARGS;
    }
    fprintf(docfs, "\n");
    fprintf(docfs, "Runs the contour advection program.\n");
    fprintf(docfs, "\n");
    fprintf(docfs, "syntax: ctraj_contour [-h tstep] [-k nfine] [-c arc] [-m dsmin] [-s dsmax]\n");
    fprintf(docfs, "                 [-i t0] [-f tf] [vfield-arguments] infile outfile\n");
    fprintf(docfs, "\n");
    fprintf(docfs, "where:\n");
    fprintf(docfs, "  vfield-arguments are arguments supplied to the velocity field loader\n");
    fprintf(docfs, "  infile   initial contour\n");
    fprintf(docfs, "  outfile  simulation results\n");
    fprintf(docfs, "\n");
    fprintf(docfs, "options:\n");
    ctraj_optargs(docfs, "Vcmshk0NifOoM?", 1);
    fprintf(docfs, "\n");
    fprintf(docfs, "Format of binary file:\n");
    fprintf(docfs, "  - 4 byte 'magic' number\n");
    fprintf(docfs, "  - 4 byte integer header containing number of records\n");
    fprintf(docfs, "  Each record has the following layout:\n");
    fprintf(docfs, "    - 2 byte integer year\n");
    fprintf(docfs, "    - 2 byte integer month\n");
    fprintf(docfs, "    - 2 byte integer day\n");
    fprintf(docfs, "    - 2 byte integer hour\n");
    fprintf(docfs, "    - 2 byte integer minute\n");
    fprintf(docfs, "    - 4 byte float second\n");
    fprintf(docfs, "    - 4 byte integer number of points in contour\n");
    fprintf(docfs, "      all longitudes are stored contiguously as floats\n");
    fprintf(docfs, "      followed by latitudes, NOT as lon-lat pairs\n");
    fprintf(docfs, "\n");
    fprintf(docfs, "  Contour advection files may be queried using:\n");
    fprintf(docfs, "  > bev2xy\n");
    fprintf(docfs, "\n");
    fprintf(docfs, "  and created using:\n");
    fprintf(docfs, "  > c0\n");
    fprintf(docfs, "  > xy2bev\n");
    fprintf(docfs, "  e.g.:\n");
    fprintf(docfs, "  > c0 45 | xy2bev c45deg.bev\n");
    fprintf(docfs, "\n");
    if (vfield!=NULL) vfield->help(docfs);
    return err;
  }

  if (vfield==NULL) {
    fprintf(stderr, "ctraj_contour: unable to load velocity field\n");
    return FILE_READ_ERROR;
  }

  if (vfield->ndim() != 2) {
    fprintf(stderr, "ctraj_contour: contour advection only works for two-dimensional velocity fields\n");
    return PARAMETER_OUT_OF_RANGE;
  }

  //for (int i=0; i<argc0; i++) printf("%s ", argv[i]);
  //printf("\n");

  //parse the argument list:
  infile=argv[1];
  outfile=argv[2];

  //determine the metric:
  if (typeid(*vfield)==typeid(ctraj_vfield_standard<float>)) {
    metric=((ctraj_vfield_standard<float> *) vfield)->get_metric();
    dflag=1;
  } else if (typeid(*vfield)==typeid(ctraj_vfield_2d<float>)) {
    metric=new Cart_t<float>(2);
    dflag=1;
  } else {
    metric=new Cart_t<float>(2);
  }

  //convert degrees to radians:
  tarc=tarc*M_PI/180;

  //initialize the contour:
  if (maxspac <= minspac) maxspac=10*minspac;
  contour=new ctraj_boundary_list<double, float>(vfield, metric, tarc, minspac, maxspac, dflag);
  //contour->set_parm(ts, nrk, tarc, minspac, maxspac);

  //get the initial contour:
  infs=fopen(infile, "r");
  if (infs==NULL) {
    fprintf(stderr, "Unable to open file, %s, for reading\n", infile);
    exit(UNABLE_TO_OPEN_FILE_FOR_READING);
  }
  fread(&magic, sizeof(magic), 1, infs);
  if (magic != MAGIC) {
    printf("Wrong file type: %s\n", infile);
    return 2;
  }
  fread(&nc, sizeof(nc), 1, infs);
  contour->read(infs);
  fclose(infs);

  //figure out dates:
  if (flag[9]) t0=(char *) optarg[9];
  if (t0!=NULL) ind1=vfield->get_tind(t0);
  //unless set otherwise, we start at the date supplied with the initial
  //contour:
  if (ind1<0) ind1=contour->gett();

  //figure out integration interval:
  if (n<0) {
    if (flag[10]) {
      tcur=vfield->get_tind((char *) optarg[10]);
      n=abs((ind1-tcur)/ts);
      if (tcur < ind1 && ts > 0) ts=-ts;
    } else {
      if (ts>0) {
        n=(vfield->maxt()-ind1)/ts;
      } else {
        //count down from 0-index:
        n=-ind1/ts;
      }
    }
    n++;
  }

  printf("Time step: %f\n", ts);
  printf("Number of steps: %d\n", n);
  printf("Number of intermediate steps: %d\n", nrk);
  printf("Threshold arc: %f\n", tarc);
  printf("Minimum spacing: %f\n", minspac);
  printf("Maximum spacing: %f\n", maxspac);

  if (oflag) contour->wrap_off();

  //set the start time:
  contour->sett(ind1);

  //open the output file and write the headers:
  outfs=fopen(outfile, "w");
  fwrite(&magic, sizeof(magic), 1, outfs);
  nrec=n/writeint+1;
  n++;
  fwrite(&nrec, sizeof(nrec), 1, outfs);
//  contour->print(stdout);
  npts=contour->fix();
  contour->write(outfs);

  for (int32_t i=1; i<n; i++) {
    tcur=contour->advance(ts/nrk, nrk);
    vfield->get_t(tcur, tstring);
    printf("%4d %s %d\n", i, tstring, npts);
//    contour->print(stdout);
    if (i % writeint == 0) contour->write(outfs);
    if (maxnpt > 0 && npts >= maxnpt) {
      nrec=i/writeint+1;
      break;
    }
    fflush(outfs);
    npts=contour->fix();
  }

  if (maxnpt > 0) {
    fseek(outfs, sizeof(magic), SEEK_SET);
    fwrite(&nrec, sizeof(nrec), 1, outfs);
  }

  fclose(outfs);

  if (typeid(*vfield)!=typeid(ctraj_vfield_standard<float>)) delete metric;
  delete vfield;
  delete contour;

  return 0;

}
