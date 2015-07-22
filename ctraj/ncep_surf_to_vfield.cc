#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "netcdfcpp.h"

#include "ctraj_defaults.h"

#include "time_class.h"
#include "simple_temp.h"
#include "dependent_swap.h"
#include "composite_dataset.h"

#include "default_grid.h"
#include "coordtran.h"
#include "ncep_util.h"

#define FNAMELEN 40

//syntax: ecmwf_sigma_to_vfield t1 t2 level outfile [rmax [ngrid]]
//where:
//  t1		Date to start import
//  t2		Date to finish
//  level	Sigma level (+ indicates N. hemisphere, - S.)
//  outfile	Output file
//  rmax	(Side-length)/2
//  ngrid	Number of grid points per side

int main(int argc, char *argv[]) {
  char *outfile;		//output file name
  FILE *outfs;			//output file stream

  time_class t1, t2;
  int yr1, yr2;
  int it10, it2f;

  ind_type nt;			//number of time grids
  ind_type nlat;
  ind_type nlon;

  float rmax;			//(side-length)/2
  ind_type ngrid;			//number of grid points per side

  float val;
  float xval, yval, r;		//for holding x and y grid values
  float tval, pval;		//temperature, pressure
  float qval;			//final value
  double c1val, c2val;		//for holding interpolation coeffs.
  int hemi;			//hemisphere
  double tint;

  //variables for input files:
  int year;

  char *path;
  char *tfile;
  char *pfile;
  char fname[FNAMELEN];

  NcFile *nc_T;
  NcFile *nc_p;

  composite_dataset all(1);	//holds the data to write to the output file
  long loc;

  simple<float> *xgrid;		//x grid
  simple<float> *ygrid;		//y grid
  simple<time_class> *tgrid;		//time grid
  simple<time_class> *tdum;		//time grid

  simple<float> *lon;		//ecmwf lon. grid
  simple<float> *lat;		//ecmwf lat. grid

  dependent<float> *pp;		//ncep u-field
  dependent<float> *tt;		//ncep temperature

  dependent_swap<float> *sfc;	//transformed surface field

  dependent<interpol_index> *c1;	//interpolation coefficients (longitude)
  dependent<interpol_index> *c2;	//interpolation coefficients (latitude)

  char c;
  int ncon;
  char tstring[30];

  long nlen;

  int tflag=0;
  int wflag=0;

  ngrid=NGRID;
  rmax=SIDELENGTH;

  path=new char [3];
  strcpy(path, "./");

  hemi=1;

  while ((c = getopt(argc, argv, "tsr:n:p:")) != -1) {
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
          fprintf(stderr, "Warning: garbled command option argument: -n %s", optarg);
          exit(2);
        }
        break;
      case ('p'):
	nlen=strlen(optarg);
	delete [] path;
        path=new char[nlen+2];
	strcpy(path, optarg);
	if (path[nlen-1] != '/') strcat(path, "/");
        break;
      case ('s'):
        hemi=-1;
        break;
      case('t'):
	tflag=1;
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
  
  if (argc < 3) {
    fprintf(stdout, "\n");
    fprintf(stdout, "Syntax:  ncep_surf_to_vfield [-r rmax] [-n ngrid] [-p path]\n");
    fprintf(stdout, "         [-t] [-s] t1 t2 outfile\n");
    fprintf(stdout, "\n");
    fprintf(stdout, "where:\n");
    fprintf(stdout, "  t1	=	date to begin import\n");
    fprintf(stdout, "  t2	=	last date to import\n");
    fprintf(stdout, "  outfile	=  name of output file\n");
    fprintf(stdout, "  rmax	=  (side-length)/2 [=%7.0f]\n", rmax);
    fprintf(stdout, "  ngrid	=  number of grid points per side [=%d]\n", ngrid);
    fprintf(stdout, "  path    =   directory containing input files [%s]\n", path);
    fprintf(stdout, "\n");
    fprintf(stdout, "flags:\n");
    //fprintf(stdout, "  -w      =       import vertical velocities\n", path);
    fprintf(stdout, "  -s      =  Southern hemisphere\n", path);
    fprintf(stdout, "  -t      =  use theta levels (default is pressure)\n", path);
    return 1;
  }

  //parse the command line arguments:
  t1.read_string(argv[0]);
  t2.read_string(argv[1]);

  //hemi is used as a flag for the command-line option
  //but we still keep the negative level as another means
  outfile=argv[2];
  if (argc>3) {
    sscanf(argv[3], "%f", &rmax);
    if (argc>4) {
      sscanf(argv[4], "%d", &ngrid);
    }
  }

  char tstr1[25], tstr2[25];
  t1.write_string(tstr1);
  t2.write_string(tstr2);
  printf("Converting data from %s to %s\n", tstr1, tstr2);
  printf("Writing to file: %s\n", outfile);

  yr1=t1.year();
  yr2=t2.year();

  tfile=new char[strlen(path)+FNAMELEN];
  pfile=new char[strlen(path)+FNAMELEN];

  //open files to get the header data:
  printf("Reading time grids...\n");
  strcpy(pfile, path);
  strcat(pfile, "pres.sfc");
  tgrid=get_ncep_tgrid(pfile, t1, t2, it10, it2f);

  //read in the grids for the NCEP data:
  printf("Reading ncep grids...\n");
  sprintf(pfile+strlen(path)+8, ".%4.4d.nc", yr1);
  nc_p=new NcFile(pfile);
  get_ncep_dim(nc_p, "lon", lon);
  get_ncep_dim(nc_p, "lat", lat);
  delete nc_p;

  nlon=lon->nel()-1;
  nlat=lat->nel();
  nt=tgrid->nel();
  printf("%d longitude x %d latitude\n", 
		  nlon, nlat);

  printf("Longitude grid:\n");
  lon->print();
  printf("Latitude grid:\n");
  lat->print();

  printf("Time grid:\n");
  tgrid->print();
  
  //x- and y-grids:
  printf("Generating interpolation coefficients\n");
  xgrid=new simple<float>(-rmax, rmax, ngrid);
  ygrid=new simple<float>(-rmax, rmax, ngrid);

  c1=new dependent<interpol_index>(xgrid, ygrid);
  c2=new dependent<interpol_index>(xgrid, ygrid);

  //calculate the interpolation coefficients:
  intcoeff(lon, lat, xgrid, ygrid, c1, c2, hemi);

  //create the datasets to hold the ecmwf data:
  pp=new dependent<float>(lon, lat);
  if (tflag) tt=new dependent<float>(lon, lat);  

  //create the output data and group it into the "all" variable:
  printf("Creating data structure to hold output data...\n");
  loc=all.add_var("xgrid");
  all.cvar(loc, (dataset *) xgrid);
  loc=all.add_var("ygrid");
  all.cvar(loc, (dataset *) ygrid);
  loc=all.add_var("time");
  all.cvar(loc, (dataset *) tgrid);

  sfc=new dependent_swap<float>(xgrid, ygrid, tgrid);

  loc=all.add_var("sfc");
  all.cvar(loc, (dataset *) sfc);
	  
  //note: this whole process is rather awkward
  //...there must be a better way of doing it...

  //start to read and write the data:
  outfs=fopen(outfile, "w");
  all.write(outfs);

  //let's make this really un-ambigous:
  //there's the time index for the input files
  //and the time index for the output files...
  ind_type time_ind=0;		//this is the time index for the output file...

  for (int year=yr1; year<=yr2; year++) {

    //create the file names:
    if (tflag) {
      strcpy(tfile, path);
      sprintf(fname, "air.sfc.%4.4d.nc", year);
      strcat(tfile, fname);
      printf("Opening file, %s...\n", tfile);
      nc_T=new NcFile(tfile);
    }
    strcpy(pfile, path);
    sprintf(fname, "pres.sfc.%4.4d.nc", year);
    strcat(pfile, fname);

    //open all the files:
    printf("Opening file, %s...\n", pfile);
    nc_p=new NcFile(pfile);

    ind_type it1=0;
    ind_type it2=ncep_nt(nc_p)-1;

    printf("Found %d time grids\n", it2);

    if (year==yr1) it1=it10;
    if (year==yr2) it2=it2f;

    //this for loop counts out the time grids in the input files...
    for (ind_type it=it1; it<=it2; it++) {
      time_class ttmpa, ttmpb;

      tgrid->get(ttmpa, time_ind);
      ttmpa.write_string(tstring);
      printf("%d Reading in grid %s...\n", time_ind, tstring);

      //time interval should always be the same value, 
      //but we calculate it anyways:
      if (time_ind == 0) {
        tgrid->get(ttmpb, (ind_type) 1);
        tint=ttmpb.diff(ttmpa);
      } else if (time_ind == nt-1) {
        tgrid->get(ttmpb, nt-2);
        tint=ttmpa.diff(ttmpb);
      } else {
        tgrid->get(ttmpa, time_ind-1);
        tgrid->get(ttmpb, time_ind+1);
        tint=ttmpb.diff(ttmpa)/2;
      }

      //read in the data from the ncep files:
      //printf("Interpolating to desired level...\n");
      if (tflag) {
      }

      //printf("velocity conversion: %f\n", vunit);
      for (ind_type j=0; j<ngrid; j++) {
        for (ind_type i=0; i<ngrid; i++) {
          //do the interpolation:
          c1->get(c1val, i, j);
          c2->get(c2val, i, j);
          pp->interpol(pval, c1val, c2val);

          if (tflag) {
            tt->interpol(tval, c1val, c2val);
            qval=tval*pow(P0/pval, KAPPA);
          } else {
            qval=pval;
          }

          sfc->cel(qval, i, j, 0, time_ind);
	}
      }
      time_ind++;
    }

    if (tflag) delete nc_T;
    delete nc_p;

  }

  //finish:
  printf("Deleting:\n");
  printf("p & t\n");
  delete sfc;

  //must not be moved above deletion of variables u and v...
  printf("Closing outfile\n");
  fclose(outfs);

  printf("Deleting:\n");
  printf("arrays\n");
  delete [] path;
  delete [] tfile;
  delete [] pfile;

  printf("pp & tt\n");
  delete pp;
  if (tflag) delete tt;

  printf("c1 & c2\n");
  delete c1;
  delete c2;

  printf("grids\n");
  delete xgrid;
  delete ygrid;
  delete tgrid;

  printf("lon & lat\n");
  delete lon;
  delete lat;
}

