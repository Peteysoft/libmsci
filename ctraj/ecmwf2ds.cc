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
#include "coordtran.h"
#include "ncep_util.h"
#include "ctraj_vfield_standard.h"

using namespace libpetey;
using namespace datasets;
using namespace ctraj;

#define FNAMELEN 14

int get_ecmwf_lev(const char *tmpfile, 
		const char *var,		//requested variable
		const char *levtype,		//level type (={pt, pl})
		float level,			//level
		time_class date,
		dependent<float> *ff, 		//returned field
		int dflag);			//delete temporary files?

int main(int argc, char *argv[]) {
  char *outfile[2];		//output file name

  char *path;

  float lev;			//vertical level
  char levtype[3]="pl";

  time_class t1, t2;
  int yr1, yr2;
  double dy1, dy2;
  long it10, it2f;		//these are long for compatibility with netcdf

  ind_type nt;			//number of time grids
  ind_type nlat=ECMWF_NLAT;
  ind_type nlon=ECMWF_NLON;

  float rmax;			//(side-length)/2
  ind_type ngrid;			//number of grid points per side

  float val;
  float xval, yval, r;		//for holding x and y grid values
  float uval, vval;		//for holding velocity vectors
  float utval, vtval;		//for holding transformed v vectors
  float wval;			//vertical velocity
  double c1val, c2val, c3val;	//for holding interpolation coeffs.
  double tint;
  float vunit;			//unit conversion for velocity

  //variables for input files:
  int year;

  char fname[FNAMELEN];

  FILE *p1u;
  FILE *p2u;

  simple<time_class> *tgrid;		//time grid
  simple<time_class> *tdum;		//time grid

  simple<float> *lon;		//ecmwf lon. grid
  simple<float> *lat;		//ecmwf lat. grid

  dependent<float> *uu;		//ncep u-field
  dependent<float> *vv;		//ncep v-field
  dependent<float> *ww;		//ncep w-field
  dependent<float> *ww1;	//ncep w-field (untransformed)
  dependent<float> *dthdp;	//derivative of potential temp. vs. pressure

  dependent_swap<float> *u;	//transformed u-field
  dependent_swap<float> *v;	//transformed v-field
  dependent_swap<float> *w;	//transformed w-field

  dependent<interpol_index> *c3;	//interpolation coefficients (vertical)

  char c;
  size_t ncon;
  char tstring[30];

  size_t nlen;

  int tflag=0;
  int wflag=0;
  int dflag;

  //for parsing command options:
  void *optargs[20];
  int flags[20];
  int64_t page_size=0;

  FILE *ps;

  ctraj_vfield_standard<float> vconv;

  int err=0;

  ngrid=NGRID;
  rmax=SIDELENGTH;

  optargs[2]=(void *) &ngrid;
  optargs[4]=(void *) &rmax;
  optargs[3]=&page_size;

  path=new char [3];
  strcpy(path, "./");

  //parse command options:
  argc=parse_command_opts(argc, argv, "ifnBrTwp?", "%s%s%d%ld%g%%%s%", 
		  optargs, flags, OPT_WHITESPACE);
  if (argc < 0) exit(21);

  if (flags[0]) t1.read_string((char *) optargs[0]);
  if (flags[1]) t2.read_string((char *) optargs[1]);
  wflag=flags[6];
  if (flags[5]) {
    strcpy(levtype, "pt");
    if (wflag) {
      fprintf(stderr, "ecmwf2ds: Vertical velocities not available on theta levels\n");
      wflag=0;
      err=10;
    }
  }

  if (flags[7]) {
    delete [] path;
    nlen=strlen((char *) optargs[9]);
    path=new char[nlen+2];
    strcpy(path, (char *) optargs[9]);
    if (path[nlen-1] != '/') strcat(path, "/");
    //dflag=2*(1-flags[11]);
    dflag=1;
  } else {
    //dflag=1+2*(1-flags[11]);
    dflag=3;
  }

  if (flags[8] || argc != 4) {
    FILE *docfs;
    int err;
    if (flags[10]) {
      docfs=stdout;
      err=0;
    } else {
      docfs=stderr;
      err=INSUFFICIENT_COMMAND_ARGS;
    }
    fprintf(docfs, "\n");
    fprintf(docfs, "Syntax:  ecmwf2ds [-r rmax] [-n ngrid] [-?] \n");
    fprintf(docfs, "               [-T] [--] [-+] [-w] [-i t1] [-f t2] \n");
    fprintf(docfs, "              level outfile\n");
    fprintf(docfs, "\n");
    fprintf(docfs, "where:\n");
    fprintf(docfs, "  level    =  level to import\n");
    fprintf(docfs, "  Sfile  =  name of output file for Southern hemisphere\n");
    fprintf(docfs, "  Nfile  =  name of output file for Northern hemisphere\n");
    fprintf(docfs, "\n");
    fprintf(docfs, "options:\n");
    ctraj_optargs(docfs, "ifnrTwBp?");
    fprintf(docfs, "if path (-p) is specified, grib files are saved for future use.\n");
    fprintf(docfs, "\n");
    return err;
  }

  if (flags[0]==0) {
    fprintf(stderr, "Initial date argument (-i) required\n");
    exit(1);
  }
  if (flags[1]==0) {
    fprintf(stderr, "Final date argument (-f) required\n");
    exit(1);
  }

  //parse the command line arguments:
  //sscanf(argv[0], "%d", &year);
  sscanf(argv[1], "%f", &lev);
  outfile[0]=argv[2];
  outfile[1]=argv[3];

  char tstr1[25], tstr2[25];
  t1.write_string(tstr1);
  t2.write_string(tstr2);
  if (t2 < t1) {
    fprintf(stderr, "ecmwf2ds: start date (%s) must be before end data (%s)\n", tstr1, tstr2);
    exit(420);
  }
  printf("Converting data from %s to %s at level %f\n", tstr1, tstr2, lev);
  printf("sidelength = %g; ngrid = %d\n", rmax, ngrid);
  printf("Writing to files: %s %s\n", outfile[0], outfile[1]);

  printf("Generating time grids...\n");
  yr1=t1.year();
  yr2=t2.year();
  dy1=t1.doy();
  dy2=t2.doy();
  t1.add(floor(dy1*4)/4-dy1);
  t2.add(ceil(dy2*4)/4-dy2);
  nt=((double) t2 - (double) t1)/0.25+1;
  tgrid=new simple<time_class>(t1, t2, nt);
  vconv.init(outfile, page_size, "w");
  vconv.set_outgrid(lev, ngrid, rmax);
  vconv.set_tgrid(*tgrid);

  //read in the grids for the NCEP data:
  printf("Generating ECMWF grids...\n");
  lon=new simple<float>(0., 360., nlon+1);
  lat=new simple<float>(-180., 180., nlat);
  vconv.set_ingrid(*lon, *lat);

  printf("%d longitude x %d latitude\n", 
		  nlon, nlat);

  printf("Longitude grid:\n");
  lon->print();
  printf("Latitude grid:\n");
  lat->print();

  printf("Time grid:\n");
  tgrid->print();
  
  //create the datasets to hold the ecmwf data:
  uu=new dependent<float>(lon, lat);  
  vv=new dependent<float>(lon, lat);

  vconv.write();

  if (wflag) {
    ww=new dependent<float>(lon, lat);
    dthdp=new dependent<float>(lon, lat);
  }

  for (int it=0; it<nt; it++) {
    time_class ttmpa;

    tgrid->get(ttmpa, it);

    //read in the data from the ncep files:
    //printf("Interpolating to desired level...\n");
    //temporary files:
    err=get_ecmwf_lev(path, "131.128", levtype, lev, ttmpa, uu, dflag);
    if (err != 0) exit(err);
    err=get_ecmwf_lev(path, "132.128", levtype, lev, ttmpa, vv, dflag);
    if (err != 0) exit(err);
    if (wflag) {
      err=get_ecmwf_lev(path, "135.128", levtype, lev, ttmpa, ww, dflag);
      if (err != 0) exit(err);
    }
    vconv.add_field(it, uu, vv);
	
  }

  //finish:

  printf("uu & vv\n");
  delete uu;
  delete vv;

  if (wflag) {
    delete ww;
    if (tflag) {
      delete dthdp;
      delete ww1;
    }
  }

  printf("grids\n");
  delete tgrid;

  printf("lon & lat & pres\n");
  delete lon;
  delete lat;

  delete [] path;

  return err;
}

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

float ctraj_glob_ecmwf_pl[37]={1, 2, 3, 5, 7, 10, 20, 30, 50, 70, 100, 125,
		150, 175, 200, 225, 250, 300, 350, 400, 450, 500, 550, 600, 650,
		700, 750, 775, 800, 825, 850, 875, 900, 925, 950,  975, 1000};
long ctraj_glob_ecmwf_npl=37;

long ctraj_glob_ecmwf_npt=15;
float ctraj_glob_ecmwf_pt[15]={265, 275, 285, 300, 315, 330, 350, 370, 395, 430,
		475, 530, 600, 700, 850};

int fetch_ecmwf_field(const char *fbase,
		const char *var,
		const char *levtype,
		float level,
		time_class date,
		int dflag=0)
{
  char *command;
  
  //date fields:
  short year, mon, day, hour, min;
  float second;

  int err;
  struct stat filestat;

  date.get_fields(year, mon, day, hour, min, second);

  command=new char[200+strlen(fbase)];
  sprintf(command, "%s.grib", fbase);
  stat(command, &filestat);

  if (!S_ISREG(filestat.st_mode)) {

    sprintf(command, "get_ecmwf.py %4.4d%2.2d%2.2d %3d %s %5.0f %s \"%s.grib\"", 
		    year, mon, day, hour, levtype, level, var, fbase);
    printf("%s\n", command);
    err=system(command);
    if (err != 0) goto finish;
  }

  sprintf(command, "wgrib -d all -nh -o \"%s.dat\" \"%s.grib\"", 
		      fbase, fbase);
  printf("%s\n", command);
  err=system(command);

  if (dflag) {
    sprintf(command, "rm -f \"%s.grib\"", fbase);
    printf("%s\n", command);
    system(command);
  }

finish:
  delete [] command;

  return err;
}

int get_ecmwf_lev(const char *path,		//path for storing temporary file...
		const char *var,		//requested variable
		const char *levtype,		//level type (={pt, pl})
		float level,			//level
		time_class date,
		dependent<float> *ff,		//returned field
		int dflag)
{
  FILE *p1u;
  FILE *p2u;

  char *command;
  char *fbase1;
  char *fbase2;

  //date fields:
  short year, mon, day, hour, min;
  float second;

  float *lev;
  long nlev;

  double zint;
  long znum;
  double frac;

  long n=ECMWF_NLON*ECMWF_NLAT;

  float data1[n];
  float data2[n];

  int err=0;
  int nread;

  int ddflag=dflag % 2;		//delete .dat file on exit
  int dgflag=dflag/2;		//delete grib file on exit

  //format precision specifiers:
  int p1, p2;

  if (strcmp(levtype, "pl") == 0) {
    lev=ctraj_glob_ecmwf_pl;
    nlev=ctraj_glob_ecmwf_npl;
  } else if (strcmp(levtype, "pt") == 0) {
    lev=ctraj_glob_ecmwf_pt;
    nlev=ctraj_glob_ecmwf_npt;
  } else {
    return 420;
  }

  zint=interpolate(lev, nlev, level, -1);
  if (zint < 0 || zint > nlev-1) {
    fprintf(stderr, "ecmwf2ds/get_ecmwf_lev: level (%f) out of range [%f, %f]\n", level, lev[0], lev[nlev-1]);
    return 420;
  }
   
  date.get_fields(year, mon, day, hour, min, second);
  znum=(long) zint;
  frac=zint-znum;
  command=new char[200+strlen(path)];
  fbase1=new char[200+strlen(path)];
  //create name of temporary files:
  p1=log10(lev[znum])+1;
  if ((int) lev[znum] == lev[znum]) p2=0; else p2=2;
  sprintf(fbase1, "%s%4.4d%2.2d%2.2d%2.2d.%s.%s.%*.*f", 
		    path, year, mon, day, hour, var, levtype, p1, p2, lev[znum]);
  //check to see if the file exists already:
  sprintf(command, "%s.dat", fbase1);
  if (ddflag) {
    p2u=NULL;
  } else {
    p2u=fopen(command, "r");
  }
  if (p2u == NULL) {
    fetch_ecmwf_field(fbase1, var, levtype, lev[znum], date, dgflag);
    p2u=fopen(command, "r");
    if (p2u == NULL) {
      fprintf(stderr, "get_ecmwf_lev: could not open data file: %s.dat\n", fbase1);
      err=101;
      goto finish;
    }
  }

  nread=fread(data1, sizeof(float), n, p2u);
  fclose(p2u);
  if (nread != n) {
    fprintf(stderr, "get_ecmwf_lev: insufficient records in data file: %s.dat\n", fbase1);
    err=111;
    goto finish;
  }

  if (frac == 0) {
    int k=0;
    //we have to reverse the latitudes:
    for (int i=ECMWF_NLAT-1; i>=0; i--) {
      for (int j=0; j<ECMWF_NLON; j++) {
        ff->cel(data1[k], j, i);
	//printf("%f\n", data1[k]);
        k++;		//do this the easy way...
      }
      //printf("\n");
    }
    //printf("\n");
  } else {

    fbase2=new char[200+strlen(path)];
    //create name of temporary files:
    p1=log10(lev[znum+1])+1;
    if ((int) lev[znum+1] == lev[znum+1]) p2=0; else p2=2;
    sprintf(fbase2, "%s%4.4d%2.2d%2.2d%2.2d.%s.%s.%*.*f", 
		    path, year, mon, day, hour, var, levtype, p1, p2, lev[znum+1]);
    //check to see if the file exists already:
    sprintf(command, "%s.dat", fbase2);
    if (ddflag) {
      p2u=NULL;
    } else {
      p2u=fopen(command, "r");
    }
    if (p2u == NULL) {
      fetch_ecmwf_field(fbase2, var, levtype, lev[znum+1], date, dgflag);
      p2u=fopen(command, "r");
      if (p2u == NULL) {
        fprintf(stderr, "get_ecmwf_lev: could not open data file: %s.dat\n", fbase2);
        err=101;
        goto finish;
      }
    }

    nread=fread(data2, sizeof(float), n, p2u);
    fclose(p2u);
    if (nread != n) {
      fprintf(stderr, "get_ecmwf_lev: insufficient records in data file: %s.dat\n", fbase1);
      err=111;
      goto finish;
    }

    int k=0;
    //we have to reverse the latitudes:
    for (int i=ECMWF_NLAT-1; i>=0; i--) {
      for (int j=0; j<ECMWF_NLON; j++) {
        ff->cel(frac*data2[k]+(1-frac)*data1[k], j, i);
	//printf("%f ", frac*data2[k]+(1-frac)*data1[k]);
        k++;		//do this the easy way...
      }
      //printf("\n");
    }
    //printf("\n");

  }

  for (int i=0; i<ECMWF_NLAT; i++) {
    float val;
    ff->get(val, 0, i);
    ff->cel(val, ECMWF_NLON, i);
  }

finish:
  if (ddflag) {
    sprintf(command, "rm -f \"%s.dat\"", fbase1);
    system(command);
    if (frac != 0) {
      sprintf(command, "rm -f \"%s.dat\"", fbase2);
      system(command);
    }
  }

  delete [] command;
  delete [] fbase1;
  if (frac != 0) delete [] fbase2;

  return err;

}

