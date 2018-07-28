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
#include "composite_dataset.h"

#include "ctraj_vfield_standard.h"
#include "coordtran.h"
#include "ncep_util.h"

#define FNAMELEN 14


using namespace libpetey;
using namespace datasets;
using namespace ctraj;

int main(int argc, char *argv[]) {
  FILE *docfs=stdout;
  char *outfile[2];		//output file name

  float lev;			//sigma level

  time_class t1, t2;
  int yr1, yr2;
  long it10, it2f;		//these are long for compatibility with netcdf

  ind_type nt;			//number of time grids
  ind_type nlat;
  ind_type nlon;
  ind_type nz;

  float rmax;			//(side-length)/2
  ind_type ngrid;			//number of grid points per side

  float wval;			//vertical velocity
  double c3val;			//for holding interpolation coeffs.

  //variables for input files:
  int year;

  char *path;
  char *tfile;
  char *ufile;
  char *vfile;
  char *wfile;
  char fname[FNAMELEN];

  NcFile *nc_T;
  NcFile *nc_u;
  NcFile *nc_v;
  NcFile *nc_w;

  simple<time_class> *tgrid;		//time grid
  simple<time_class> *tdum;		//time grid

  simple<float> *lon;		//ecmwf lon. grid
  simple<float> *lat;		//ecmwf lat. grid
  simple<float> *pres;		//pressure grid

  dependent<float> *uu;		//ncep u-field
  dependent<float> *vv;		//ncep v-field

  dependent<double> *c3;

  char c;
  size_t ncon;
  char tstring[30];

  size_t nlen;

  int tflag=0;
  int wflag=0;

  //for parsing command options:
  void *optargs[20];
  int flags[20];
  int32_t page_size=0;

  char *command;
  FILE *ps;
  int yr1a, yr2a;

  int err=0;

  ctraj_vfield_standard<float> vconv;

  ngrid=NGRID;
  rmax=SIDELENGTH;

  optargs[2]=(void *) &ngrid;
  optargs[4]=(void *) &rmax;
  optargs[7]=&page_size;

  path=new char [3];
  strcpy(path, "./");

  //parse command options:
  argc=parse_command_opts(argc, argv, "ifnprTwB?", "%s%s%d%s%g%%%ld%", 
		  optargs, flags, OPT_WHITESPACE);
  if (argc < 0) exit(21);

  if (flags[0]) t1.read_string((char *) optargs[0]);
  if (flags[1]) t2.read_string((char *) optargs[1]);
  if (flags[3]) {
    delete [] path;
    nlen=strlen((char *) optargs[3]);
    path=new char[nlen+2];
    strcpy(path, (char *) optargs[3]);
    if (path[nlen-1] != '/') strcat(path, "/");
  }
  tflag=flags[5];
  wflag=flags[6];

  //help page print block:
  if (flags[8] || argc < 4) {
    int err;
    if (flags[8]) {
      docfs=stdout;
      err=0;
    } else {
      docfs=stderr;
      err=INSUFFICIENT_COMMAND_ARGS;
    }
    fprintf(docfs, "\n");
    fprintf(docfs, "Syntax:  nc2ds [-r rmax] [-n ngrid] [-p path]\n");
    fprintf(docfs, "               [-T] [-w] [-i t1] [-f t2] \n");
    fprintf(docfs, "              level Sfile Nfile\n");
    fprintf(docfs, "\n");
    fprintf(docfs, "where:\n");
    fprintf(docfs, "  level    =  level to import (sign determines hemisphere)\n");
    fprintf(docfs, "  Sfile    =  name of S. hemi. output file\n");
    fprintf(docfs, "  Nfile    =  name of N. hemi. output file\n");
    fprintf(docfs, "\n");
    fprintf(docfs, "options:\n");
    ctraj_optargs(docfs, "ifnprTwB?");
    fprintf(docfs, "\n");
    return err;
  }

  //parse the command line arguments:
  //sscanf(argv[0], "%d", &year);
  sscanf(argv[1], "%f", &lev);

  outfile[0]=argv[2];
  outfile[1]=argv[3];

  //this is a bit convoluted:
  //should be more complete, but if you supply no dates, what do you expect??
  if (flags[0] != 1 || flags[1] != 1) {
    command=new char[strlen(path)+25];
    sprintf(command, "%s %s %s", "ls ", path, "| grep air\\.....\\.nc$");
    fprintf(docfs, "%s\n", command);
    system(command);
    //air temperature files:
    ps=popen(command, "r");
    if (fgets(fname, FNAMELEN, ps)==NULL) {
      fprintf(stderr, "nc2ds: No air temperature data found\n");
      exit(101);
    }
    sscanf(fname+4, "%d", &year);
    yr1=year;
    yr2=year;
    while (feof(ps) == 0) {
      fgets(fname, FNAMELEN, ps);
      sscanf(fname+4, "%d", &year);
      if (year<yr1) yr1=year;
      if (year>yr2) yr2=year;
    }
    fclose(ps);
    //u-wind files:
    sprintf(command, "%s %s %s", "ls ", path, "| grep uwnd\\.....\\.nc$");
    ps=popen(command, "r");
    if (fgets(fname, FNAMELEN, ps)==NULL) {
      fprintf(stderr, "nc2ds: No zonal wind data found\n");
      exit(101);
    }
    sscanf(fname+5, "%d", &year);
    yr1a=year;
    yr2a=year;
    while (feof(ps) == 0) {
      fgets(fname, FNAMELEN, ps);
      sscanf(fname+5, "%d", &year);
      if (year<yr1a) yr1a=year;
      if (year>yr2a) yr2a=year;
    }
    fclose(ps);
    if (yr1a>yr1) yr1=yr1a;
    if (yr2a<yr2) yr2=yr2a;
    //v-wind files:
    sprintf(command, "%s %s %s", "ls ", path, "| grep vwnd\\.....\\.nc$");
    ps=popen(command, "r");
    if (fgets(fname, FNAMELEN, ps)==NULL) {
      fprintf(stderr, "nc2ds: No meridional wind data found\n");
      exit(101);
    }
    sscanf(fname+5, "%d", &year);
    yr1a=year;
    yr2a=year;
    while (feof(ps) == 0) {
      fgets(fname, FNAMELEN, ps);
      sscanf(fname+5, "%d", &year);
      if (year<yr1a) yr1a=year;
      if (year>yr2a) yr2a=year;
    }
    fclose(ps);
    if (yr1a>yr1) yr1=yr1a;
    if (yr2a<yr2) yr2=yr2a;
    if (flags[0]==0) t1.init(yr1, 1, 1, 0, 0, 0);
    if (flags[2]==0) t2.init(yr2, 12, 31, 23, 59, 59);
  }

  char tstr1[25], tstr2[25];
  t1.write_string(tstr1);
  t2.write_string(tstr2);
  if (t2 < t1) {
    fprintf(stderr, "nc2ds: start date (%s) must be before end data (%s)\n", tstr1, tstr2);
    exit(420);
  }
  fprintf(docfs, "Converting data from %s to %s at level %f\n", tstr1, tstr2, lev);
  fprintf(docfs, "sidelength = %g; ngrid = %d; page size = %ld\n", rmax, ngrid, page_size);
  fprintf(docfs, "Writing to files: %s %s\n", outfile[0], outfile[1]);

  yr1=t1.year();
  yr2=t2.year();

  tfile=new char[strlen(path)+FNAMELEN];
  ufile=new char[strlen(path)+FNAMELEN];
  vfile=new char[strlen(path)+FNAMELEN];

  //open files to get the header data:
  fprintf(docfs, "Reading time grids...\n");
  strcpy(tfile, path);
  strcat(tfile, "air");
  tgrid=get_ncep_tgrid(tfile, t1, t2, it10, it2f);
  sprintf(tfile+strlen(path)+3, ".%4.4d.nc", yr1);

  //read in the grids for the NCEP data:
  fprintf(docfs, "Reading ncep grids...\n");
  nc_T=new NcFile(tfile);
  get_ncep_grid(nc_T, lon, lat, pres, tdum);
  delete tdum;		//don't need this...
  delete nc_T;

  nlon=lon->nel()-1;
  nlat=lat->nel();
  nz=pres->nel();
  nt=tgrid->nel();
  fprintf(docfs, "%d longitude x %d latitude x %d level\n", 
		  nlon, nlat, nz);

  fprintf(docfs, "Longitude grid:\n");
  lon->print();
  fprintf(docfs, "Latitude grid:\n");
  lat->print();
  fprintf(docfs, "NCEP pressure grid:\n");
  pres->print();

  fprintf(docfs, "Time grid:\n");
  tgrid->print();
  
  //x- and y-grids:
  fprintf(docfs, "Generating interpolation coefficients\n");
  vconv.init(outfile, page_size, "w");
  vconv.set_outgrid(lev, ngrid, rmax);
  vconv.set_tgrid(*tgrid);
  vconv.set_ingrid(*lon, *lat);

  //vertical interpolation coeff:
  if (tflag==0) {
    c3val=nz-pres->interp(lev)-1;
  } else {
    c3=new dependent<interpol_index>(lon, lat);
  }

  //create the datasets to hold the ecmwf data:
  uu=new dependent<float>(lon, lat);  
  vv=new dependent<float>(lon, lat);

  //if dataset is set up for paging, write out grids before filling in
  //velocity data:
  if (page_size!=-1) vconv.write();

  //let's make this really un-ambigous:
  //there's the time index for the input files
  //and the time index for the output files...
  ind_type time_ind=0;		//this is the time index for the output file...

  for (int year=yr1; year<=yr2; year++) {

    //create the file names:
    strcpy(tfile, path);
    sprintf(fname, "air.%4.4d.nc", year);
    strcat(tfile, fname);
    strcpy(ufile, path);
    sprintf(fname, "uwnd.%4.4d.nc", year);
    strcat(ufile, fname);
    strcpy(vfile, path);
    sprintf(fname, "vwnd.%4.4d.nc", year);
    strcat(vfile, fname);

    //open all the files:
    fprintf(docfs, "Opening file, %s...\n", tfile);
    nc_T=new NcFile(tfile);
    fprintf(docfs, "Opening file, %s...\n", ufile);
    nc_u=new NcFile(ufile);
    fprintf(docfs, "Opening file, %s...\n", vfile);
    nc_v=new NcFile(vfile);

    if (wflag) {
      strcpy(wfile, path);
      sprintf(fname, "omega.%4.4d.nc", year);
      fprintf(docfs, "Opening file, %s...\n", wfile);
      nc_w=new NcFile(wfile);
    }

    ind_type it1=0;
    ind_type it2=ncep_nt(nc_T)-1;

    fprintf(docfs, "Found %d time grids\n", it2);

    if (year==yr1) it1=it10;
    if (year==yr2) it2=it2f;

    //this for loop counts out the time grids in the input files...
    for (ind_type it=it1; it<=it2; it++) {
      time_class ttmpa;

      tgrid->get(ttmpa, time_ind);
      ttmpa.write_string(tstring);
      fprintf(docfs, "%d Reading in grid %s...\n", time_ind, tstring);

      //read in the data from the ncep files:
      //printf("Interpolating to desired level...\n");
      if (tflag) {
        fprintf(docfs, "Calculating vertical interpolation coefficients\n");
        get_ncep_theta_interp(nc_T, lev/pow(P0, KAPPA), it, c3);
        if (wflag) {
          //get_dthdp(nc_T, lev/pow(P0, KAPPA), it, c3, dthdp);
          //get_ncep_theta_level(nc_w, "omega", it, c3, ww1);
          //ww=ww1*dthdp;
        }
        fprintf(docfs, "Interpolating to desired theta level\n");
        get_ncep_theta_level(nc_u, "uwnd", it, c3, uu);
        get_ncep_theta_level(nc_v, "vwnd", it, c3, vv);
      } else {
        get_ncep_tz(nc_u, "uwnd", it, c3val, uu);
        get_ncep_tz(nc_v, "vwnd", it, c3val, vv);
      }

      vconv.add_field(time_ind, uu, vv);

      /*
      printf("u\n");
      float val;
      for (int j=0; j<nlat; j++) {
        for (int i=0; i<nlon; i++) {
          uu->get(val, i, j);
          printf("%g\n", val);
        }
      }

      printf("v\n");
      for (int j=0; j<nlat; j++) {
        for (int i=0; i<nlon; i++) {
          vv->get(val, i, j);
          printf("%g\n", val);
        }
      }
      */

      time_ind++;
    }

    delete nc_T;
    delete nc_u;
    delete nc_v;

    if (wflag) delete nc_w;
  }

  //finish:
  //if dataset is not set up for paging, it is held completely in RAM and
  //must be written at the end:
  if (page_size==-1) vconv.write();

  //must not be moved above deletion of variables u and v...
  fprintf(docfs, "Deleting:\n");
  fprintf(docfs, "arrays\n");
  delete [] path;
  delete [] tfile;
  delete [] ufile;
  delete [] vfile;

  fprintf(docfs, "uu & vv\n");
  delete uu;
  delete vv;

  fprintf(docfs, "c3\n");
  if (tflag) delete c3;

  if (wflag) {
    //delete w;
    //delete ww;
    if (tflag) {
      //delete dthdp;
      //delete ww1;
    }
  }

  fprintf(docfs, "grids\n");
  delete tgrid;

  fprintf(docfs, "lon & lat & pres\n");
  delete lon;
  delete lat;
  delete pres;

  return err;
}

