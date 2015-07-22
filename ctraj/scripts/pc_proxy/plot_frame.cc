//lon-lat ascii to azimuthal-equidistant binary...
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>

#include "simple_temp.h"
#include "dependent_temp.h"

#include "coordtran.h"
#include "field_tran.h"

#define MAX_COMMAND_LEN 200

int main(int argc, char **argv) {

  char *infile;
  char *outfile;
  char *cfile=NULL;
  long index;

  int cwflag=0;		//write to contour file
  int gflag=0;		//use geometric progression
  int bflag=0;		//have bottom
  int tflag=0;		//have top
  float bottom, top;	//range for contours...

  char command[MAX_COMMAND_LEN];	//command for calling GMS

  FILE *fs;

  long ncon;

  float *dum;

  float rmax;		//side length/2
  ind_type ngrid;		//grid points per side
  ind_type nlon, nlat;	//lon-lat grid
  char c;

  //datasets:
  simple<float> *lon;
  simple<float> *lat;

  long fsize;		//file size
  long recsize;		//record size
  long nrec;		//number of records

  ngrid=NGRID_Q;
  rmax=SIDELENGTH_Q;

  nlon=NLON;
  nlat=NLAT;

  while ((c = getopt(argc, argv, "r:x:y:n:b:t:g")) != -1) {
    switch (c) {
      case ('n'):
        ncon=sscanf(optarg, "%d", &ngrid);
        if (ncon != 1) {
          fprintf(stderr, "Warning: garbled command option argument: -n %s", optarg);
          exit(2);
        }
        break;
      case ('r'):
        ncon=sscanf(optarg, "%f", &rmax);
        if (ncon != 1) {
          fprintf(stderr, "Warning: garbled command option argument: -r %s", optarg);
          exit(2);
        }
        break;
      case ('a'):
	aflag=1;
	break;
      case ('x'):
        ncon=sscanf(optarg, "%d", &nlon);
        if (ncon != 1) {
          fprintf(stderr, "Warning: garbled command option argument: -x %s", optarg);
          exit(2);
        }
	break;
      case ('y'):
        ncon=sscanf(optarg, "%d", &nlat);
        if (ncon != 1) {
          fprintf(stderr, "Warning: garbled command option argument: -y %s", optarg);
          exit(2);
        }
	break;
      case ('b'):
        ncon=sscanf(optarg, "%f", &bottom);
        if (ncon != 1) {
          fprintf(stderr, "Warning: garbled command option argument: -n %s", optarg);
          exit(2);
        }
	cwflag=1;
	bflag=1;
        break;
      case ('t'):
        ncon=sscanf(optarg, "%f", &top);
        if (ncon != 1) {
          fprintf(stderr, "Warning: garbled command option argument: -n %s", optarg);
          exit(2);
        }
	cwflag=1;
	tflag=1;
        break;
      case ('c'):
	ncon=strlen(optarg);
        if (ncon == 0) {
          fprintf(stderr, "Warning: garbled command option argument: -n %s", optarg);
          exit(2);
        }
	cfile=new char[ncon+1];
	strcpy(cfile, optarg);
        break;
      case ('g'):
	gflag=1;
	cwflag=1;
	break;
      case ('?'):
             fprintf(stderr, "Unknown option: -%c -- ignored\n", optopt);
             break;
      default:
             fprintf(stderr, "Error parsing command line\n");
             exit(2);
             break;
    }
  }
  argc-=optind;
  argv+=optind;

  if (argc < 1) {
    fprintf(stdout, "\n");
    fprintf(stdout, "Syntax:  plot_frame [-a] [-r rmax] [-n ngrid] [-x nlon] [-y nlat]\n");
    fprintf(stdout, "                    [-b bottom] [-t top] [-c cfile] [-g]\n");
    fprintf(stdout, "                    infile index [outfile]\n");
    fprintf(stdout, "\n");
    fprintf(stdout, "Plots a single tracer field.**\n");
    fprintf(stdout, "\n");
    fprintf(stdout, "where:\n");
    fprintf(stdout, "  infile  =  input file name, binary dump of an array of vectors\n");
    fprintf(stdout, "  index   =  0-based index of field to extract\n");
    fprintf(stdout, "\n");
    fprintf(stdout, "For input data (azimuthal equidistant coords):\n");
    fprintf(stdout, "  rmax    =  (side-length)/2 [%7.0f]\n", rmax);
    fprintf(stdout, "  ngrid   =  grid points per side [%d]\n", ngrid);
    fprintf(stdout, "For output data (lon-lat coords):\n");
    fprintf(stdout, "  nlon    =  number of longitude grids [%d]\n", nlon);
    fprintf(stdout, "  nlat    =  number of latitude grids [%d]\n\n", nlat);
    fprintf(stdout, "For setting contours:\n");
    fprintf(stdout, "  bottom  =  bottom contour\n");
    fprintf(stdout, "  top     =  top contour\n");
    fprintf(stdout, "  cfile   =  GMT-compatible file containing contours\n");
    fprintf(stdout, "  -g      =  signifies geometric progression\n");
    fprintf("\n");
    fprintf("** Note: require Generic Mapping Tools (GMT) package.\n");
    return 1;
  }

  infile=argv[0];
  if (argc>1) sscanf(argv[1], "%d", &index); else index=-1;
  outfile=argv[2];

  fs=fopen(infile, "r");

  fseek(fs, 0, SEEK_END);
  fsize=ftell(fs);
  recsize=ntot*2*sizeof(float);
  nrec=fsize/recsize;

  if (fsize%recsize!=0) {
    fprintf(stderr, "Not an even number of records:\n");
    fprintf(stderr, "     %ld mod %ld = %ld\n", fsize, recsize, fsize%recsize);
  }

  if (index < 0 || index >=nrec) {
    printf("%ld\n", nrec);
    exit(0);
  }

  fseek(fs, index*recsize, SEEK_SET); //should do some range-checking here...

  fread(vec, sizeof(float), ntot*2, fs);
    fread(readq_n, sizeof(float), ntot, fs);
    //map this onto a rectangular grid:
    for (ind_type i=0; i<ngrid; i++) {
      for (ind_type j=0; j<ngrid; j++) {
        map->get(mapval, i, j);
	//(this is not the most efficient way of doing things...)
	if (mapval > -1) {
          qn->cel(readq_n[mapval], i, j);
	  qs->cel(readq_s[mapval], i, j);
	  //printf("%d %f\n", mapval, readq_s[mapval]);
	}
      }
    }
  } else {
    qs->read(fs);
    qn->read(fs);
  }
  fclose(fs);

  //create the output grids:
  nlats=nlat/2;

  lon=new simple<float>(0., 360.-360./nlon, nlon);
  lat=new simple<float>(-90., 90., nlat);

  //lon->print();
  //lat->print();

  c1=new dependent<double>(lon, lat);
  c2=new dependent<double>(lon, lat);
  intcoeff2(lon, lat, xgrid, ygrid, c1, c2, -1);

  //Southern Hemisphere:
  //perform the interpolation and write the result to the output file:
  for (long j=0; j<nlats; j++) {
    for (long i=0; i<nlon; i++) {
      c1->get(c1val, i, j);
      c2->get(c2val, i, j);
      qs->interpol(val, c1val, c2val);
      //printf("%d %d\n", i, j);
      //printf("%lg %lg\n", c1val, c2val);
      //printf("%f\n", val);
      printf("%g\n", val);
    }
  }

  //Northern Hemisphere:
  //create the output grids and interpolation coefficients:
  intcoeff2(lon, lat, xgrid, ygrid, c1, c2, 1);

  //perform the interpolation and write the result to the output file:
  for (long j=nlats; j<nlat; j++) {
    for (long i=0; i<nlon; i++) {
      //printf("%d %d\n", i, j);
      c1->get(c1val, i, j);
      c2->get(c2val, i, j);
      qn->interpol(val, c1val, c2val);
      printf("%g\n", val);
    }
  }

  delete c1;
  delete c2;

  delete qn;
  delete qs;
  delete lon;
  delete lat;

  if (aflag != 1) delete map;

  delete xgrid;
  delete ygrid;

  delete [] readq_n;
  delete [] readq_s;

}
