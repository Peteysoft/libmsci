//
// Copywrite 2004 Peter Mills.  All rights reserved.
//

#include <assert.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <getopt.h>
#include <stdlib.h>

#include "time_class.h"
#include "peteys_tmpl_lib.h"
#include "simple_temp.h"
#include "dependent_swap.h"
#include "composite_dataset.h"

#include "dependent_intc.h"

#include "tcoord_defs.h"
#include "coordtran.h"
#include "default_grid.h"

#include "sparse.h"

//global datasets containing the velocity field:
simple<float> *x_grid_n;
simple<float> *y_grid_n;
simple<float> *x_grid_s;
simple<float> *y_grid_s;
//simple<float> *level;
simple<time_class> *time_grid;

dependent_swap<float> *u_n;
dependent_swap<float> *v_n;
dependent_swap<float> *u_s;
dependent_swap<float> *v_s;

//northern hemisphere derivatives:
void derivs_n(float tind, float xvec[], float dxdt[], long n) {
  interpol_index xint, yint;

  xint=x_grid_n->interp(xvec[0]);
  yint=y_grid_n->interp(xvec[1]);
  //printf("%f %f\n", xint, yint);
  u_n->interpol(dxdt[0], xint, yint, 0, (interpol_index) tind);
  v_n->interpol(dxdt[1], xint, yint, 0, (interpol_index) tind);

  //printf("Velocity at t= %f : (%f, %f)\n", tind, dxdt[0], dxdt[1]);

}

//southern hemisphere derivatives:
void derivs_s(float tind, float xvec[], float dxdt[], long n) {
  interpol_index xint, yint;

  xint=x_grid_s->interp(xvec[0]);
  yint=y_grid_s->interp(xvec[1]);
  //printf("%f %f\n", xint, yint);
  u_s->interpol(dxdt[0], xint, yint, 0, (interpol_index) tind);
  v_s->interpol(dxdt[1], xint, yint, 0, (interpol_index) tind);

  //printf("Velocity at t= %f : (%f, %f)\n", tind, dxdt[0], dxdt[1]);

}

#define MAXLL 200
#define DATE_WIDTH 23

#define RMAX 10000

//syntax:
//tracer4_step1 initfile
//
// where initfile is the initialization file
//
// format of initfile:
//
//nfile		Binary file containing velocity field in the N. hemi.
//sfile		Binary file containing the velocity field in the S. hemi.
//outfile	Output file
//mapfile	File containing the mapping from the two tracer fields to a vector
//t0		Date to start integration
//nt		Number of (Eulerian) time steps
//tstep		Non-dimensional time (Eulerian, i.e. coarse) step
//nfine		Number of Runge-Kutta steps between each time step
//np		Number of grid points per side

int main(int argc, char *argv[]) {
  char *initfile;
  FILE *initfs;
  char line[MAXLL];

  //main data file names:
  char nfile[MAXLL];
  char sfile[MAXLL];
  char outfile[MAXLL];

  char c;
  long ncon;

  //composite dataset containing velocity field:
  composite_dataset ndata;
  composite_dataset sdata;
  FILE *nfs;				//file stream for velocity field
  FILE *sfs;
  long loc, dum;			//location of the datasets (members of vdata)

  //time step info:
  ind_type n;				//number of time steps
  int32_t nfine=4;				//number of trajectory time steps
  float tstep_fine;			//time step for trajectory calc.
  float tstep=1;			//time step for tracer field

  //general date variables:
  time_class date1, date2, date3;		//a date
  char date_str[DATE_WIDTH];		//a date as a string

  //time indices:
  double *tind;				//vector of time indices
  ind_type lt;				//time index as integer

  //stuff for integration:
  float x0[2];				//initial conditions
  float **result;			//integrated values

  FILE *outfun;				//output file unit

  ind_type nlon=NLON, nlat=NLAT;
  simple<float> *longrid;		//tracer x grid
  simple<float> *latgrid;		//tracer y grid

  //dummy tracer field for calculating interpolation coefficients:
  dependent_intc *tracer;

  //sparse matrix class instance for holding the matrix of interpolation coeffs:
  sparse_matrix map;

  //for checking the consistancy of the two files:
  simple<float> *check1;
  simple<float> *check2;
  simple<time_class> *check3;

  //intermediate values:
  float r;	//radius
  float dx;
  sub_1d_type ind, k;
  sub_1d_type ind0, ind1, ind2, ind3;

  int qflag=0;		//just print out the dates...

  if (argc < 2) {
    printf("Purpose: two-dimensional, 'semi-Lagrangian' tracer simulation driven \n");
    printf("by globally gridded wind fields.  For each time step, outputs a sparse \n");
    printf("matrix defining the mapping from one tracer field to the next\n");
    printf("\n");
    printf("syntax:\n");
    printf("tracerLL [-Q] [-x nlon] [-y nlat] [-t nt] [-i t0] [-h dt] [-f nfine]\n");
    printf("		  sfile nfile outfile\n");
    printf("\n");
    printf(" where:\n");
    printf("\n");
    printf("sfile:	Binary file containing velocity field in the N. hemi.\n");
    printf("nfile:	Binary file containing the velocity field in the S. hemi.\n");
    printf("outfile:	Output file\n");
    printf("\n");
    printf("t0:		Date to start integration\n");
    printf("		(default is bottom time grid in velocity field)\n");
    printf("nt:		Number of (Eulerian) time steps\n");
    printf("		(default is number of time grids in velocity field)\n");
    printf("h:		Non-dimensional time (Eulerian, i.e. coarse) step [%f]\n", tstep);
    printf("nfine:	Number of Runge-Kutta steps between each time step [%d]\n", nfine);
    printf("nlon:	Number of longitude grids [%d]\n", nlon);
    printf("nlat:	Number of latitude grids [%d]\n", nlat);
    return 1;
  }

  while ((c = getopt(argc, argv, "Qx:y:i:t:h:f:")) != -1) {
    switch (c) {
      case ('x'):
        ncon=sscanf(optarg, "%d", &nlon);
        if (ncon != 1) {
          fprintf(stderr, "Warning: garbled command option argument: -x %s", optarg);
          exit(2);
        }
        break;
      case ('y'):
        ncon=sscanf(optarg, "%f", &nlat);
        if (ncon != 1) {
          fprintf(stderr, "Warning: garbled command option argument: -y %s", optarg);
          exit(2);
        }
        break;
      case ('Q'):
	qflag=1;
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

  initfile=argv[0];

  //parse the initialisation file:
  initfs=fopen(initfile, "r");
  fgets(line, MAXLL, initfs);
  sscanf(line, "%s", nfile);
  fgets(line, MAXLL, initfs);
  sscanf(line, "%s", sfile);
  fgets(line, MAXLL, initfs);
  sscanf(line, "%s", outfile);
  fgets(line, MAXLL, initfs);
  sscanf(line, "%s", mapfile);
  fgets(line, MAXLL, initfs);
  date1=line;
  fgets(line, MAXLL, initfs);
  sscanf(line, "%g", &tstep);
  fgets(line, MAXLL, initfs);
  sscanf(line, "%d", &n);
  fgets(line, MAXLL, initfs);
  sscanf(line, "%d", &nfine);
  fgets(line, MAXLL, initfs);
  if (np <= 0) sscanf(line, "%d", &np);

  printf("tstep: %f, number of time steps: %d, Runge-Kutta steps: %d, grid size: %d\n",
  		tstep, n, nfine, np);

  //get the N. hemi. velocity fields:
  nfs=fopen(nfile, "r");
  ndata.read(nfs);
  loc=ndata.search_var("xgrid", dum);
  x_grid_n=(simple<float> *) ndata.get_var(loc);
  loc=ndata.search_var("ygrid", dum);
  y_grid_n=(simple<float> *) ndata.get_var(loc);
  loc=ndata.search_var("time", dum);
  time_grid=(simple<time_class> *) ndata.get_var(loc);
  loc=ndata.search_var("u", dum);
  u_n=(dependent_swap<float> *) ndata.get_var(loc);
  loc=ndata.search_var("v", dum);
  v_n=(dependent_swap<float> *) ndata.get_var(loc);

  //get the S. hemi. velocity fields:
  sfs=fopen(sfile, "r");
  sdata.read(sfs);
  loc=sdata.search_var("xgrid", dum);
  x_grid_s=(simple<float> *) sdata.get_var(loc);
  loc=sdata.search_var("ygrid", dum);
  y_grid_s=(simple<float> *) sdata.get_var(loc);
  loc=sdata.search_var("u", dum);
  u_s=(dependent_swap<float> *) sdata.get_var(loc);
  loc=sdata.search_var("v", dum);
  v_s=(dependent_swap<float> *) sdata.get_var(loc);

  //check the consistency of the two fields:
  loc=ndata.search_var("zgrid", dum);
  check1=(simple<float> *) sdata.get_var(loc);
  loc=sdata.search_var("zgrid", dum);
  check2=(simple<float> *) sdata.get_var(loc);
  assert(*check1 == *check2);
  loc=sdata.search_var("time", dum);
  check3=(simple<time_class> *) sdata.get_var(loc);
  assert(*time_grid==*check3);

  //figure out the initial time index:
  tind=new double[n+1];
  //figure out its location in relation to the velocity field:
  tind[0]=time_grid->interp(date1);
  tstep_fine=-tstep/nfine;

  for (long i=1; i<=n; i++) tind[i]=tind[i-1]+tstep;

  //initialize the dummy tracer field:
  //tracer gridding:
  /* really, what the fuck IS this shit???
  qgrid=new float [np];
  dx=2.0*SIDELENGTH_Q/(np-4.0);
  maxr=SIDELENGTH_Q+sqrt(2)*dx;
  for (long i=0; i<np; i++) {
    qgrid[i]=dx*(i-1.5)-RMAX;
  }
  x_qgrid=new simple<float>(qgrid, np, 0);
  y_qgrid=new simple<float>(qgrid, np, 0);*/

  x_qgrid=new simple<float>(-sidelength, sidelength, np);
  y_qgrid=new simple<float>(-sidelength, sidelength, np);
  //x_qgrid->print();

  tracer=new dependent_intc(x_qgrid, y_qgrid);

  map_map=new dependent<sub_1d_type>(x_qgrid, y_qgrid);

  nocorner_map(x_qgrid, y_qgrid, map_map, nmap);

  //generate the mapping of a hemispherical tracer field to a vector
  //(omitting everything below the equator):
  mfs=fopen(mapfile, "w");
  for (ind_type i=0; i<np; i++) {
    for (ind_type j=0; j<np; j++) {
      float x, y;
      map_map->get(ind, i, j);
      x_qgrid->get(x, i);
      y_qgrid->get(y, j);
      fprintf(mfs, "%8d %8d %8d %12.4f %12.4f\n", i, j, ind, x, y);
    }
  }
  fclose(mfs);

  //initialize the sparse matrix:
  map.extend(nmap*8);

  //initialize the vector of results for the Runge-Kutta integrations:
  result=new float * [nfine+1];
  for (long i=0; i<=nfine; i++) result[i]=new float[2];

  //open the output file and write the headers:
  outfun=fopen(outfile, "w");
//  fwrite(&maxr, 1, sizeof(maxr), outfun);
//  fwrite(&np, 1, sizeof(np), outfun);
//  fwrite(&n, 1, sizeof(n), outfun);

  for (long it=0; it<n; it++) {
    float x0[2];		//initial cond. for traj.
    float xf, yf;		//final cond.
    float val;			//interpolated value
    interpol_index xyind[2];	//interpolation indices
    sub_1d_type row_sub;
    sub_1d_type sub[4];		//1d subscripts
    double weight[4];		//interpolation coeffs
    short hemi;			//hemisphere

    //get the date:
    lt=(ind_type) tind[it];
    time_grid->get(date1, lt);
    time_grid->get(date2, lt+1);
    date3=date2-date1;
    date2=date3*(tind[it]-(float) lt);
    date3=date1+date2;
    date3.write_string(date_str);

    printf("%d %s\n", it, date_str);

    if (qflag) continue;

    map.reset(nmap*2, nmap*2);

    for (ind_type i=0; i<np; i++) {
      for (ind_type j=0; j<np; j++) {
        //check to see that the point falls within the useable regions:
	map_map->get(row_sub, i, j);
	if (row_sub != -1) {
          //printf("%d\n", row_sub);
          //use the current grid point as the initial condition:
          x_qgrid->get(x0[0], i);
          y_qgrid->get(x0[1], j);
//	  printf("(%f, %f)\n", x0[0], x0[1]);

          //Southern hemisphere:
          //do a Runge-Kutta integration:
          rk_dumb((float) tind[it+1], x0, 2L, tstep_fine, nfine, result, &derivs_s);
          //get the final position:
          xf=result[nfine][0];
          yf=result[nfine][1];

	  //"fix" it:
	  hemi=-1;
	  tcoord_fix(xf, yf, hemi);
          hemi=(hemi+1)/2;

	  //printf("Initial: (%f, %f); final: (%f, %f)\n", x0[0], x0[1], xf, yf);
          //find the interpolated value in the previous tracer field:
          xyind[0]=x_qgrid->interp(xf);
          xyind[1]=y_qgrid->interp(yf);
	  //printf("(%f, %f)\n", xyind[0], xyind[1]);

          tracer->interpol_coeff(xyind, sub, weight);
	  /*printf("%d %d %d %d\n", sub[0], sub[1], sub[2], sub[3]);
	  printf("(%d %d)\n", sub[0]/np, sub[0] % np);
	  printf("(%d %d)\n", sub[1]/np, sub[1] % np);
	  printf("(%d %d)\n", sub[2]/np, sub[2] % np);
	  printf("(%d %d)\n", sub[3]/np, sub[3] % np);*/
	  //map_map->get(ind0, sub[0]%np, sub[0]/np);
	  //printf("%d\n", ind0);

	  map_map->get_1d(ind0, sub[0]);
	  //printf("%d\n", ind0);
	  map_map->get_1d(ind1, sub[1]);
	  map_map->get_1d(ind2, sub[2]);
	  map_map->get_1d(ind3, sub[3]);

	  assert(ind0 != -1);
	  assert(ind1 != -1);
	  assert(ind2 != -1);
	  assert(ind3 != -1);

          map.add_el(weight[0], row_sub, ind0+hemi*nmap);
          map.add_el(weight[1], row_sub, ind1+hemi*nmap);
          map.add_el(weight[2], row_sub, ind2+hemi*nmap);
          map.add_el(weight[3], row_sub, ind3+hemi*nmap);

          //Northern hemisphere:
          //do a Runge-Kutta integration:
          rk_dumb((float) tind[it+1], x0, 2L, tstep_fine, nfine, result, &derivs_n);
          //get the final position:
          xf=result[nfine][0];
          yf=result[nfine][1];

	  //"fix" it:
	  hemi=1;
	  tcoord_fix(xf, yf, hemi);
	  hemi=(hemi+1)/2;

//	  printf("Initial: (%f, %f); final: (%f, %f)\n", x0[0], x0[1], xf, yf);
          //find the interpolated value in the previous tracer field:
          xyind[0]=x_qgrid->interp(xf);
          xyind[1]=y_qgrid->interp(yf);
	  //printf("(%f, %f)\n", xyind[0], xyind[1]);

          tracer->interpol_coeff(xyind, sub, weight);
	  /*printf("%d %d %d %d\n", sub[0], sub[1], sub[2], sub[3]);
	  printf("(%d %d)\n", sub[0]/np, sub[0] % np);
	  printf("(%d %d)\n", sub[1]/np, sub[1] % np);
	  printf("(%d %d)\n", sub[2]/np, sub[2] % np);
	  printf("(%d %d)\n", sub[3]/np, sub[3] % np);*/

	  row_sub+=nmap;
	  map_map->get_1d(ind0, sub[0]);
	  map_map->get_1d(ind1, sub[1]);
	  map_map->get_1d(ind2, sub[2]);
	  map_map->get_1d(ind3, sub[3]);

	  assert(ind0 != -1);
	  assert(ind1 != -1);
	  assert(ind2 != -1);
	  assert(ind3 != -1);

          map.add_el(weight[0], row_sub, ind0+hemi*nmap);
          map.add_el(weight[1], row_sub, ind1+hemi*nmap);
          map.add_el(weight[2], row_sub, ind2+hemi*nmap);
          map.add_el(weight[3], row_sub, ind3+hemi*nmap);
	}
      }
    }

    //write the date and the tracer field to a file:
//    fwrite(date_str, 1, sizeof(char)*DATE_WIDTH, outfun);
    map.write(outfun);
    //map.print(stdout);

  }

  fclose(outfun);
  fclose(nfs);
  fclose(sfs);

  //clean up:
  delete [] tind;
  for (long i=0; i<=nfine; i++) delete [] result[i];
  delete [] result;

  /*
  printf("u_n\n");
  delete u_n;
  printf("v_n\n");
  delete v_n;
  printf("x_grid_n\n");
  delete x_grid_n;
  printf("y_grid_n\n");
  delete y_grid_n;

  printf("u_s\n");
  delete u_s;
  printf("v_s\n");
  delete v_s;
  printf("x_grid_s\n");
  delete x_grid_s;
  printf("y_grid_s\n");
  delete y_grid_s;

  printf("time_grid\n");
  delete time_grid;
  */

  //printf("tracer\n");
  delete tracer;
  //printf("map_map\n");
  delete map_map;

  //printf("x_qgrid\n");
  delete x_qgrid;
  //printf("y_qgrid\n");
  delete y_qgrid;

}
