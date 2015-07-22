#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "error_codes.h"
#include "time_class.h"
#include "parse_command_opts.h"
#include "simple_temp.h"
#include "dependent_swap.h"
#include "composite_dataset.h"

#include "ctraj_defaults.h"
#include "ctraj_vfield_standard.h"

using namespace libpetey;
using namespace datasets;
using namespace ctraj;

#define FNAMELEN 14

int main(int argc, char *argv[]) {
  FILE *docfs=stdout;
  FILE *infs=stdin;
  char *outfile[2];		//output file name
  char mode[3]="w";

  float lev;			//sigma level
  float val;

  time_class t1, t2, dt;

  ind_type nt;			//number of time grids
  ind_type nlat=NLAT;
  ind_type nlon=NLON;
  ind_type nz;
  ind_type dim[2];
  ind_type j0, jf, dj;

  float rmax;			//(side-length)/2
  ind_type ngrid;			//number of grid points per side

  //variables for input files:
  ind_type time_ind=0;		//this is the time index for the output file...

  char *path;

  simple<time_class> *tgrid;		//time grid

  simple<float> *lon;		//ecmwf lon. grid
  simple<float> *lat;		//ecmwf lat. grid

  dependent<float> *uu;		//ncep u-field
  dependent<float> *vv;		//ncep v-field

  char c;
  size_t ncon;
  char tstring[30];

  size_t nlen;

  int tflag=0;
  int wflag=0;

  //for parsing command options:
  void *optargs[30];
  int flags[30];
  int64_t page_size=0;

  char *command;
  FILE *ps;
  int yr1a, yr2a;

  int err=0;

  ctraj_vfield_standard<float> vconv;

  ngrid=NGRID;
  rmax=SIDELENGTH;

  optargs[2]=(void *) &ngrid;
  optargs[4]=(void *) &rmax;
  optargs[5]=&lev;
  optargs[7]=&page_size;
  optargs[8]=&nlon;
  optargs[9]=&nlat;
  optargs[11]=&nt;

  path=new char [3];
  strcpy(path, "./");

  //parse command options:
  argc=parse_command_opts(argc, argv, "ifnprzwBxyhN?RX", "%s%s%d%s%g%g%%ld%d%d%s%d%%%", 
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

  if (flags[12] || argc < 4) {
    if (flags[12]) {
      docfs=stdout;
      err=0;
    } else {
      docfs=stderr;
      err=INSUFFICIENT_COMMAND_ARGS;
    }
    fprintf(docfs, "\n");
    fprintf(docfs, "Purpose: converts global wind data in lon-lat gridding to a native\n");
    fprintf(docfs, "binary format.  Takes as input an ASCII stream of numbers, cycling\n");
    fprintf(docfs, "through longitude first, then latitute, with the zonal wind first\n");
    fprintf(docfs, "then the meridional.  By specifying the time grid, this will insert\n");
    fprintf(docfs, "a single, 2D field into the file.  Parameters (input grid, output\n");
    fprintf(docfs, "grid, time grid, vertical level) need only be specified once,\n");
    fprintf(docfs, "when the file is first created.\n");
    fprintf(docfs, "\n");
    fprintf(docfs, "Syntax:  vdata2ds [-r rmax] [-n ngrid] [-z level]\n");
    fprintf(docfs, "               [-w] [-i t1] [-f t2] [-h dt] [-N nt-1] \n");
    fprintf(docfs, "               tind Sfile Nfile\n");
    fprintf(docfs, "\n");
    fprintf(docfs, "where:\n");
    fprintf(docfs, "  tind     =  time index\n");
    fprintf(docfs, "  Sfile    =  name of S. hemi. data file\n");
    fprintf(docfs, "  Nfile    =  name of N. hemi. data file\n");
    fprintf(docfs, "\n");
    fprintf(docfs, "options:\n");
    ctraj_optargs(docfs, "ifnprwBxyN?");
    fprintf(docfs, "  -h   time grid interval [0/0/0-6:0:0.0]\n");
    fprintf(docfs, "  -R   reverse latitude grids\n");
    fprintf(docfs, "  -X   x- (longitude) major\n");
    fprintf(docfs, "  -z   vertical level\n");
    fprintf(docfs, "\n");
    return err;
  }

  //parse the command line arguments:
  //sscanf(argv[0], "%d", &year);
  sscanf(argv[1], "%d", &time_ind);

  outfile[0]=argv[2];
  outfile[1]=argv[3];

  tgrid=NULL;

  //logic is a bit screwy, but I think I've covered all bases...
  if (flags[0]) {
    if (flags[10]) {
      dt.read_string((char *) optargs[10]);
      if (flags[11]) {
        t2=t1+dt*(double) nt;
        tgrid=new simple<time_class>(t1, t2, nt+1);
      }
      if (flags[1]) {
        nt=(t2-t1)/dt;
        tgrid=new simple<time_class>(t1, t2, nt+1);
      }
    } else if (flags[1]) {
      if (flags[11]) {
        tgrid=new simple<time_class>(t1, t2, nt+1);
      } else {
        nt=(t2-t1)/(time_class) "0/0/0-6";
        tgrid=new simple<time_class>(t1, t2, nt+1);
      }
    }
  }

  //if the user specifies the time grid, then we create a new file,
  //otherwise we assume it exists already....
  if (tgrid==NULL) {
    strcpy(mode, "r+");
  }

  lon=NULL;
  lat=NULL;
  if (strcmp(mode, "w")==0 || (flags[8] && flags[9])) {
    lon=new simple<float>(0., 360., nlon+1);
    lat=new simple<float>(-90., 90., nlat);
  }

  char tstr1[25], tstr2[25];
  t1.write_string(tstr1);
  t2.write_string(tstr2);
  if (t2 < t1) {
    fprintf(stderr, "vdata2ds: start date (%s) must be before end data (%s)\n", tstr1, tstr2);
    exit(420);
  }
  fprintf(docfs, "Converting data from %s to %s at level %f\n", tstr1, tstr2, lev);
  fprintf(docfs, "sidelength = %g; ngrid = %d\n", rmax, ngrid);
  fprintf(docfs, "Writing to files: %s %s\n", outfile[0], outfile[1]);

  if (tgrid != NULL ) {
    fprintf(docfs, "Time grid:\n");
    tgrid->print();
  }
  
  //x- and y-grids:
  vconv.init(outfile, page_size, mode);
  if (strcmp(mode, "w")==0) {
    vconv.set_outgrid(lev, ngrid, rmax);
    vconv.set_tgrid(*tgrid);
  }
  if (lon!=NULL) {
    vconv.set_ingrid(*lon, *lat);
  }

  //create the datasets to hold the input data:
  uu=vconv.empty_field();  
  vv=vconv.empty_field();  

  if (strcmp(mode, "w")==0) vconv.write(1);

  uu->get_dim(dim);

  if (flags[13]) {
    //reverse latitude grids:
    j0=dim[1]-1;
    jf=-1;
    dj=-1;
  } else {
    j0=0;
    jf=dim[1];
    dj=1;
  }

  //subtract 1 from longitude grid because they haven't been "wrapped" yet:
  if (flags[14]) {
    //longitude major:
    for (int32_t i=0; i<dim[0]-1; i++) {
      for (int32_t j=j0; j!=jf; j+=dj) {
        fscanf(infs, "%g", &val);
        uu->cel(val, i, j);
      }
    }
    for (int32_t i=0; i<dim[0]-1; i++) {
      for (int32_t j=j0; j!=jf; j+=dj) {
        fscanf(infs, "%g", &val);
        vv->cel(val, i, j);
      }
    }
  } else {
    //latitude major:
    printf("u\n");
    for (int32_t j=j0; j!=jf; j+=dj) {
      for (int32_t i=0; i<dim[0]-1; i++) {
        fscanf(infs, "%g", &val);
        uu->cel(val, i, j);
        //printf("%g\n", val);
      }
    }
    printf("v\n");
    for (int32_t j=j0; j!=jf; j+=dj) {
      for (int32_t i=0; i<dim[0]-1; i++) {
        fscanf(infs, "%g", &val);
        vv->cel(val, i, j);
        //printf("%g\n", val);
      }
    }
  }

  //here's where we "wrap" the longitude grids:
  for (int32_t j=0; j<dim[1]; j++) {
    uu->get(val, 0, j);
    uu->cel(val, dim[0]-1, j);
    vv->get(val, 0, j);
    vv->cel(val, dim[0]-1, j);
  }

  vconv.add_field(time_ind, uu, vv);

  //finish:

  fprintf(docfs, "Deleting:\n");
  fprintf(docfs, "arrays\n");

  fprintf(docfs, "uu & vv\n");
  delete uu;
  delete vv;

  if (wflag) {
    //delete w;
    //delete ww;
  }

  fprintf(docfs, "grids\n");
  delete tgrid;

  fprintf(docfs, "lon & lat\n");
  delete lon;
  delete lat;

  return err;
}

