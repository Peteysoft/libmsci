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
  char *initfile=NULL;
  FILE *initfs;

  char line[MAXLL];
  char *tok;

  //main data file names:
  char **nfile;
  char **sfile;
  char outfile[MAXLL];

  char *zfile;

  char c;
  long ncon;

  traj_3D_ojb *straj;
  traj_3D_obj *ntraj;

  //vertical grids:
  long nz;
  float *zval;
  simple<float> *zgrid;

  //time step info:
  long n;				//number of time steps
  long nfine=4;				//number of trajectory time steps
  float tstep_fine;			//time step for trajectory calc.
  float tstep=1;			//time step for tracer field

  //general date variables:
  time_class date1, date2, date3;		//a date
  char date_str[DATE_WIDTH];		//a date as a string

  //time indices:
  double tind1, tind2;
  double *tind;				//vector of time indices
  ind_type lt;				//time index as integer

  //stuff for integration:
  float x0[3];				//initial conditions
  float **result;			//integrated values

  FILE *outfun;				//output file unit

  //grid info:
  ind_type np=NGRID_Q;			//number of x and y tracer grids
  float sidelength=SIDELENGTH_Q;	//actually sidelength/2...
  ind_type nz_q;			//number of levels for tracer

  float *qgrid;
  simple<float> *x_qgrid;		//tracer x grid
  simple<float> *y_qgrid;		//tracer y grid
  simple<float> *z_qgrid;		//tracer z grid

  //dummy tracer field for calculating interpolation coefficients:
  dependent_intc *tracer;

  //sparse matrix class instance for holding the matrix of interpolation coeffs:
  sparse_matrix map;

  //intermediate values:
  float r;	//radius
  float dx;
  sub_1d_type ind, k;
  sub_1d_type ind0, ind1, ind2, ind3;

  int qflag=0;		//just print out the dates...
  
  zfile=NULL;

  while ((c = getopt(argc, argv, "Qr:n:h:k:")) != -1) {
    switch (c) {
      case ('n'):
        ncon=sscanf(optarg, "%d", &np);
        if (ncon != 1) {
          fprintf(stderr, "Warning: garbled command option argument: -n %s", optarg);
          exit(2);
        }
        break;
      case ('r'):
        ncon=sscanf(optarg, "%f", &sidelength);
        if (ncon != 1) {
          fprintf(stderr, "Warning: garbled command option argument: -r %s", optarg);
          exit(2);
        }
        break;
      case ('h'):
        ncon=sscanf(optarg, "%f", &tstep);
        if (ncon != 1) {
          fprintf(stderr, "Warning: garbled command option argument: -h %s", optarg);
          exit(2);
        }
        break;
      case ('k'):
        ncon=sscanf(optarg, "%d", &nfine);
        if (ncon != 1) {
          fprintf(stderr, "Warning: garbled command option argument: -k %s", optarg);
          exit(2);
        }
        break;
      case ('v'):
        zfile=new char [strlen(optarg)+1];
        strcpy(zfile, optarg);
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

  if (argc != 4) {
    printf("Purpose: two-dimensional, 'semi-Lagrangian' tracer simulation driven \n");
    printf("by globally gridded wind fields.  For each time step, outputs a sparse \n");
    printf("matrix defining the mapping from one tracer field to the next\n");
    printf("\n");
    printf("syntax:\n");
    printf("tracer_3D [-Q] [-r sidelengh/2] [-n ngrid] [-h dt] [-k nfine] [-v zfile]\n");
    printf("               initfile t0 tf outfile\n");
    printf("\n");
    printf(" where:\n");
    printf("   initfile  is the initialization file\n");
    printf("   t0        is the date to start integration\n");
    printf("   tf        is the date to finish integration\n");
    printf("\n");
    return 1;
  }

  strcpy(initfile, argv[0]);
  date1=argv[1];
  date2=argv[2];
  strcpy(outfile, argv[3]);

  initfs=fopen(initfile, "r");
  fgets(line, MAXLL, initfs);
  sscanf(line, "%g", &nz);
  zval=new float[nz];
  nfile=new char *[nz+1];
  sfile=new char *[nz+1];
  for (int i=0; i<nz; i++) {
    fgets(line, MAXLL, initfs);
    tok=strsep(&line, " ");
    sscanf(tok, "%f", zval+i);
    tok=strsep(&line, " ");
    nfile[i]=new char[strlen(tok)];
    strcpy(nfile[i], tok);
    tok=strsep(&line, " ");
    sfile[i]=new char[strlen(tok)];
    strcpy(sfile[i], tok);
  }
  fgets(line, MAXLL, initfs);
  tok=strsep(&line, " ");
  nfile[nz]=new char[strlen(tok)];
  strcpy(nfile[i], tok);
  tok=strsep(&line, " ");
  sfile[nz]=new char[strlen(tok)];
  strcpy(sfile[i], tok);
  fclose(initfs);

  zgrid=new simple<float>(zval, nz);
  delete [] zval;

  ntraj=new traj_3D_obj(nfile, zgrid);
  straj=new traj_3D_obj(sfile, zgrid);

  if (zfile!=NULL) {
    initfs=fopen(zfile, "r");
    fgets(line, MAXLL, initfs);
    sscanf(line, "%f", &nz_q);
    zval=new float[nz_q];

    for (int i=0; i<nz_q; i++)
      fgets(line, MAXLL, initfs);
      sscanf(line, "%f", zval+i);
    }
    z_qgrid=new simple<float>(zval, nz_q);
    delete [] zval;
  } else {
    z_qgrid=zgrid;
    nz_q=z_qgrid->nel();
  }

  //figure out the initial time index:
  tind1=time_grid->interp(date1);
  if (argc == 5) {
    tind2=time_grid->interp(date2);
    n=ceil((tind2-tind1)/tstep)+1;
  }

  printf("tstep: %f, number of time steps: %d, Runge-Kutta steps: %d, grid size: %d\n",
  		tstep, n, nfine, np);

  tind=new double[n+1];
  //figure out its location in relation to the velocity field:
  tind[0]=tind1;
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

  tracer=new dependent_intc(x_qgrid, y_qgrid, z_qgrid);

  map_map=new dependent<sub_1d_type>(x_qgrid, y_qgrid);

  nocorner_map(x_qgrid, y_qgrid, map_map, nmap);

  //initialize the sparse matrix:
  map.extend(nmap*16*nz_q);

  //open the output file and write the headers:
  outfun=fopen(outfile, "w");
//  fwrite(&maxr, 1, sizeof(maxr), outfun);
//  fwrite(&np, 1, sizeof(np), outfun);
//  fwrite(&n, 1, sizeof(n), outfun);

  for (long it=0; it<n; it++) {
    float x0[3];		//initial cond. for traj.
    float xf[3];		//final cond.
    float val;			//interpolated value
    interpol_index xyzind[3];	//interpolation indices
    sub_1d_type row_sub;
    sub_1d_type sub[8];		//1d subscripts
    double weight[8];		//interpolation coeffs
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

    map.reset(nmap*2*nz_q, nmap*2*nz_q);

    for (ind_type k=0; k<nz_q; k++) {
      z_qgrid->get(x0[2], k);
      for (ind_type i=0; i<np; i++) {
        x_qgrid->get(x0[0], i);
        for (ind_type j=0; j<np; j++) {
          //check to see that the point falls within the useable regions:
          map_map->get(row_sub, i, j);
          if (row_sub != -1) {
            //printf("%d\n", row_sub);
            //use the current grid point as the initial condition:
            y_qgrid->get(x0[1], j);
//          printf("(%f, %f)\n", x0[0], x0[1]);

            //Southern hemisphere:
            //do a Runge-Kutta integration:
            straj->integrate(x0, xf, tind[it+1], tstep_fine, nfine);

            //"fix" it:
            hemi=-1;
            tcoord_fix(xf[0], xf[1], hemi);
            hemi=(hemi+1)/2;

            //printf("Initial: (%f, %f); final: (%f, %f)\n", x0[0], x0[1], xf, yf);
            //find the interpolated value in the previous tracer field:
            xyzind[0]=x_qgrid->interp(xf);
            xyzind[1]=y_qgrid->interp(yf);
            xyzind[2]=z_qgrid->interp(zf);
            //printf("(%f, %f)\n", xyind[0], xyind[1]);

            tracer->interpol_coeff(xyzind, sub, weight);
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

//          printf("Initial: (%f, %f); final: (%f, %f)\n", x0[0], x0[1], xf, yf);
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
