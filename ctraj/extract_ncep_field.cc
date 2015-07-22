#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "netcdfcpp.h"

#include "error_codes.h"
#include "time_class.h"
#include "parse_command_opts.h"
#include "simple_temp.h"
#include "dependent_swap.h"

#include "ctraj_defaults.h"
#include "ncep_util.h"

using namespace libpetey;
using namespace datasets;
using namespace ctraj;

#define FNAMELEN 14

int main(int argc, char *argv[]) {
  FILE *docfs=stdout;
  char *var;			//output file name

  float lev;			//sigma level
  time_class date;

  int yr1, yr2;
  long it1, it2;		//these are long for compatibility with netcdf

  ind_type nt;			//number of time grids
  ind_type nlat;
  ind_type nlon;
  ind_type nz;

  double tind;
  double c3val;			//for holding interpolation coeffs.

  //variables for input files:
  int year;

  char *path;
  char *tfile;
  char *vfile;
  char fname[FNAMELEN];

  NcFile *nc_T;
  NcFile *nc_V;

  simple<time_class> *tgrid;		//time grid
  simple<float> *lon;		//ecmwf lon. grid
  simple<float> *lat;		//ecmwf lat. grid
  simple<float> *pres;		//pressure grid

  dependent<float> *vv;		//ncep field
  dependent<float> *vv1;		//ncep field

  dependent<double> *c3;	//vertical interpolation coeffs.

  float val, val1;

  char c;
  size_t ncon;
  char tstring[30];

  size_t nlen;

  int tflag=0;
  int wflag=0;

  //for parsing command options:
  void *optargs[20];
  int flags[20];

  char *command;
  FILE *ps;
  int yr1a, yr2a;

  int err=0;

  path=new char [0];
  strcpy(path, "./");

  //parse command options:
  argc=parse_command_opts(argc, argv, "pT?H", "%s%%%", 
		  optargs, flags, OPT_WHITESPACE);
  if (argc < 0) exit(21);

  if (flags[0]) {
    delete [] path;
    nlen=strlen((char *) optargs[0]);
    path=new char[nlen+2];
    strcpy(path, (char *) optargs[0]);
    if (path[nlen-1] != '/') strcat(path, "/");
  }
  tflag=flags[1];

  if (flags[2] || argc < 4) {
    int err;
    if (flags[2]) {
      docfs=stdout;
      err=0;
    } else {
      docfs=stderr;
      err=INSUFFICIENT_COMMAND_ARGS;
    }
    fprintf(docfs, "\n");
    fprintf(docfs, "Syntax:  extract_ncep_field [-T] [-p path] variable level date\n");
    fprintf(docfs, "\n");
    fprintf(docfs, "where:\n");
    fprintf(docfs, "  variable =  variable name\n");
    fprintf(docfs, "  level    =  level to extract\n");
    fprintf(docfs, "  date     =  date to extract\n");
    fprintf(docfs, "\n");
    fprintf(docfs, "options:\n");
    ctraj_optargs(docfs, "pT?");
    fprintf(docfs, "\n");
    return err;
  }

  var=argv[1];
  sscanf(argv[2], "%f", &lev);
  date=argv[3];

  tfile=new char[strlen(path)+FNAMELEN];
  vfile=new char[strlen(path)+FNAMELEN];

  //open files to get the header data:
  //assume first argument is a variable name:
  year=date.year();
  sprintf(vfile, "%s%s.%4.4d.nc", path, var, year);
  sprintf(tfile, "%sair.%4.4d.nc", path, year);

  //read in the grids for the NCEP data:
  //fprintf(docfs, "Reading ncep grids from file %s...\n", vfile);
  nc_V=new NcFile(vfile);

  if (nc_V==NULL) {
    fprintf(docfs, "extract_ncep_field: unable to open input file: %s\n", vfile);
    exit(UNABLE_TO_OPEN_FILE_FOR_READING);
  }
  get_ncep_grid(nc_V, lon, lat, pres, tgrid);

  nlon=lon->nel()-1;
  nlat=lat->nel();
  nz=pres->nel();
  nt=tgrid->nel();

  if (flags[3]) printf("%d %d\n", nlon, nlat);
  //fprintf(docfs, "%d longitude x %d latitude x %d level\n", 
  //		  nlon, nlat, nz);

  //fprintf(docfs, "Longitude grid:\n");
  //lon->print();
  //fprintf(docfs, "Latitude grid:\n");
  //lat->print();
  //fprintf(docfs, "NCEP pressure grid:\n");
  //pres->print();

  //fprintf(docfs, "Time grid:\n");
  //tgrid->print();
  
  //vertical interpolation coeff:
  if (tflag==0) {
    c3val=nz-pres->interp(lev)-1;
  } else {
    c3=new dependent<interpol_index>(lon, lat);
  }

  //create the datasets to hold the ncep data:
  vv=new dependent<float>(lon, lat);
  vv1=new dependent<float>(lon, lat);

  //fprintf(docfs, "Opening file, %s...\n", vfile);
  nc_V=new NcFile(vfile);

  tind=tgrid->interp(date);
  if (tind < 0 || tind > nt-1) {
    fprintf(stderr, "extract_ncep_field: date not contained in file\n");
    exit(PARAMETER_OUT_OF_RANGE);
  }

  it1=tind;
  it2=it1+1;

  //fprintf(docfs, "Reading in grid %lg...\n", tind);

  //read in the data from the ncep files:
  //printf("Interpolating to desired level...\n");
  if (tflag) {
    //fprintf(docfs, "Opening file, %s...\n", tfile);
    if (strcmp(tfile, vfile)==0) {
      nc_T=nc_V;
    } else {
      nc_T=new NcFile(tfile);
    }
    //fprintf(docfs, "Calculating vertical interpolation coefficients\n");
    get_ncep_theta_interp(nc_T, lev/pow(P0, KAPPA), it1, c3);
    //fprintf(docfs, "Interpolating to desired theta level\n");
    get_ncep_theta_level(nc_V, var, it1, c3, vv);
    if (it1!=tind) {
      //fprintf(docfs, "Calculating vertical interpolation coefficients\n");
      get_ncep_theta_interp(nc_T, lev/pow(P0, KAPPA), it2, c3);
      //fprintf(docfs, "Interpolating to desired theta level\n");
      get_ncep_theta_level(nc_V, var, it2, c3, vv1);
    }
    if (nc_T!=nc_V) delete nc_T;
  } else {
    get_ncep_tz(nc_V, var, it1, c3val, vv);
    if (it1!=tind) get_ncep_tz(nc_V, var, it2, c3val, vv1);
  }

  if (it1==tind) {
    for (int j=0; j<nlat; j++) {
      for (int i=0; i<nlon; i++) {
        vv->get(val, i, j);
        printf("%g\n", val);
      }
    }
  } else {
    for (int j=0; j<nlat; j++) {
      for (int i=0; i<nlon; i++) {
        vv->get(val, i, j);
        vv1->get(val1, i, j);
        printf("%g\n", val+(val1-val)*(tind-it1));
      }
    }
  }

  delete nc_V;

  //finish:

  //fprintf(docfs, "Deleting:\n");
  //fprintf(docfs, "arrays\n");
  delete [] path;
  delete [] tfile;
  delete [] vfile;

  //fprintf(docfs, "uu & vv\n");
  delete vv;
  delete vv1;

  //fprintf(docfs, "c3\n");
  if (tflag) delete c3;

  //fprintf(docfs, "grids\n");
  delete tgrid;

  //fprintf(docfs, "lon & lat & pres\n");
  delete lon;
  delete lat;
  delete pres;

  return err;
}

