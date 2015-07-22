#include <assert.h>

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_multifit.h>

#include "time_class.h"

#include "dependent_intc.h"
#include "coordtran.h"
#include "tcoord_defs.h"
#include "ctraj_defaults.h"
#include "sparse_array.h"
#include "av.h"

#include "ctraj_vfield_base.h"
#include "ctraj_tfield_base.h"

#include "meas_data.h"

//*FLAG*

// *** BE CAREFULL! ***

// THIS FILE IS AN ATTEMPTED UPDATE TO tracer_anal.cc

// DO NOT CONFUSE THE TWO

// CURRENTLY NOT USED ANYWHERE...

//*FLAG*

float ** tracer_multiply(sparse_matrix *a, int32_t N, float *v0) {
  ind_t m, n;
  float **tracer;

  a[0].dimensions(m, n);

  tracer=new float *[N+1];
  tracer[0]=new float[(N+1)*n];

  for (int32_t i=0; i<n; i++) tracer[0][i]=v0[i];

  for (int32_t i=1; i<=N; i++) {
    tracer[i]=tracer[0]+i*n;
    a[i].vect_mult(tracer[i-1], tracer[i]);
  }

  return tracer;
}

float ** tracer_multiply_renorm(sparse_matrix *a, int32_t N, float *v0) {
  ind_t m, n;
  float **tracer;
  float vmin, vmax;

  a[0].dimensions(m, n);

  tracer=new float *[N+1];
  tracer[0]=new float[(N+1)*n];

  for (int32_t i=0; i<n; i++) tracer[0][i]=v0[i];

  for (int32_t i=1; i<=N; i++) {
    tracer[i]=tracer[0]+i*n;
    vmin=tracer[i-1][0];
    vmax=tracer[i-1][0];
    for (ind_type j=1; j<n; j++) {
      if (tracer[i-1][j] > vmax) vmax=tracer[i-1][j];
      if (tracer[i-1][j] < vmin) vmin=tracer[i-1][j];
    }
    for (ind_t j=0; j<n; j++) {
      tracer[i-1][j]=(tracer[i-1][j]-vmin)/(vmax-vmin)*2-1;
    }
    a[i].vect_mult(tracer[i-1], tracer[i]);
  }

  return tracer;
}

//assumes standard side-length of 10000 km
template <class real>
void tracer_interp(ctraj_tfield_base<real> *trac, real **qall, time_class *t, int32_t nt, int32_t np, meas_data *data, long nsamp) {

  real xs, ys;
  int32_t domain;

  //for generating the interpolation coefficients:
  int32_t nint;
  int32_t *sub;
  double *intc;
  double tind;
  double frac;
  int tindi;

  //for solving the coefficient matrix:
  float q1, q2;

  nint=trac->nwt();
  sub=new int32_t[nint];
  intc=new double[nint];

  for (long i=0; i<nsamp; i++) {
    //float norm;		//normalization coefficient
    			//(in case of point near border)
			//should never be the case...
    //get the interpolation coefficients:
    if (data[i].lat < 0) hemi=-1; else hemi=1;
    tcoord2_2lonlat(xs, ys, -1, hemi, data[i].lon, data[i].lat);
    float lon, lat;
    tcoord2_2lonlat(xs, ys, 1, hemi, lon, lat);
    //printf("Transformed coords: (%f, %f) -> (%f, %f)\n", data[i].lon, data[i].lat, xs, ys);
    //printf("Transformed coords: (%f, %f) -> (%f, %f)\n", lon, lat, xs, ys);
    int_ind[0]=xgrid->interp(xs);
    int_ind[1]=ygrid->interp(ys);
    qdum->interpol_coeff(int_ind, sub, intc);
    //printf("Interpolation coefficients: %g, %g, %g, %g\n", intc[0], intc[1], intc[2], intc[3]);
    tind=interpolate(t, nt, data[i].t, last_ind);
    tindi=(int32_t) tind;
    frac=tind-(double) tindi;
    //printf("Time index: %lg\n", tind);
    //printf("1-d subscripts: %d, %d, %d, %d\n", sub[0], sub[1], sub[2], sub[3]);
    map->get_1d(sub[0], sub[0]);
    map->get_1d(sub[1], sub[1]);
    map->get_1d(sub[2], sub[2]);
    map->get_1d(sub[3], sub[3]);
    //printf("1-d subscripts: %d, %d, %d, %d\n", sub[0], sub[1], sub[2], sub[3]);

    q1=0;
    q2=0;
    //norm=0;
    for (int k=0; k<4; k++) {
      if (sub[k]>=0) {
        q1+=qall[tindi][sub[k]+(hemi+1)*nmap/2]*intc[k];
        q2+=qall[tindi+1][sub[k]+(hemi+1)*nmap/2]*intc[k];
	//norm+=int_ind[k];
      }
    }
    //printf("q1 %g; q2 %g\n", q1, q2);

    data[i].q=(q1*(1-frac)+q2*frac);//norm;
  }

  delete qdum;
  delete map;
  delete xgrid;
  delete ygrid;

}

float ** pc_proxy(sparse_matrix *matall, time_class *t, int32_t N, int32_t nall, meas_data *samp, long nsamp, int32_t nev, int32_t ncv, int cflag=0) {

  ind_t m, n;		//size of matrix
  int32_t i0;		//start point, number of matrices

  float **v;		//eigenvectors
  float *eval;		//eigenvalues

  //the array of "right" singular vectors
  float ***qall;

  //for generating the interpolation coefficients:
  float xs, ys;			//sample transformed coords
  short hemi;
  ind_type np;
  float sidelength=SIDELENGTH_Q;
  dependent_intc *qdum;
  interpol_index int_ind[3];
  sub_1d_type sub[8];
  interpol_index intc[8];
  interpol_index tind;
  interpol_index frac;
  int tindi;
  simple<float> *xgrid, *ygrid;
  dependent<ind_type> *map;
  sub_1d_type nmap;

  sparse_array<int32_t, float> *a_mat;
  sparse_array<int32_t, float> at_mat;
  sparse_array<int32_t, float> ata_mat;

  if (cflag != 0) cflag=1;

  //for solving the coefficient matrix:
  gsl_matrix *a;
  gsl_vector *x;
  gsl_vector *b;		//coefficients
  float **qvec;
  float q1, q2;
  gsl_multifit_linear_workspace *work;
  gsl_matrix *cov;
  double chisq;

  i0=0;

  //get matrix dimensions:
  matall[0].dimensions(m, n);
  assert(m==n);

  //prepare for SVD:
  //sa_dir=1;
  if (N+i0 > nall) N=nall-i0;
  /*nsparse=2*N;
  sparse_array=new sparse_matrix[N*2];
  for (ind_type i=0; i<N; i++) {
    sparse_array[i]=matall[i+i0];
    sparse_array[nsparse-i-1]=sparse_array[i];
    sparse_array[nsparse-i-1].transpose();
  }*/

  a_mat=new sparse_array<int32_t, float>(matall+i0, N);
  at_mat=*a_mat;
  at_mat.transpose();
  a_mat->mat_mult(at_mat, ata_mat);

  //perform the SVD:
  v=new float*[ncv];
  v[0]=new float[n*ncv];
  for (int i=1; i<ncv; i++) v[i]=v[0]+i*n;
  eval=new float[ncv];

  //call the fortran program:
  FORTRAN_FUNC(svd_sparse_array)(&n, &nev, &ncv, v[0], eval, (void *) (&ata_mat));

  //calculate number of gridpoints from number of tracer points
  np=2*(int) sqrt(n/M_PI/2);

  //printf("ngrid (estimate)=%d\n", np);
   
  xgrid=new simple<float>(-sidelength, sidelength, np);
  ygrid=new simple<float>(-sidelength, sidelength, np);

  qdum=new dependent_intc(xgrid, ygrid);
  map=new dependent<ind_type>(xgrid, ygrid);

  //get the "mapping"
  nocorner_map(xgrid, ygrid, map, nmap);
  assert(2*nmap==n);

  //first we multiply through to get all the tracers:
  qall=new float **[nev];
  for (int32_t i=0; i<nev; i++) {
    qall[i]=tracer_multiply(matall+i0, nall, v[i]);

  }

  //create the matrix:
  //printf("nsamp=%d, nev=%d\n", nsamp, nev);
  a=gsl_matrix_alloc(nsamp, nev+cflag);
  b=gsl_vector_alloc(nsamp);
  for (long i=0; i<nsamp; i++) {
    //we add a constant term to the matrix:
    if (cflag) gsl_matrix_set(a, i, nev, 1);

    //get the interpolation coefficients:
    if (samp[i].lat < 0) hemi=-1; else hemi=1;
    tcoord2_2lonlat(xs, ys, -1, hemi, samp[i].lon, samp[i].lat);
    int_ind[0]=xgrid->interp(xs);
    int_ind[1]=ygrid->interp(ys);
    qdum->interpol_coeff(int_ind, sub, intc);
    tind=interpolate(t, nall, samp[i].t, -1)-i0;
    tindi=(int) tind;
    frac=tind-(double) tindi;
    map->get_1d(sub[0], sub[0]);
    map->get_1d(sub[1], sub[1]);
    map->get_1d(sub[2], sub[2]);
    map->get_1d(sub[3], sub[3]);

    for (int32_t j=0; j<nev; j++) {
      q1=0;
      q2=0;
      for (int k=0; k<4; k++) {
        q1+=qall[j][tindi][sub[k]+(hemi+1)*nmap/2]*intc[k];
        q2+=qall[j][tindi+1][sub[k]+(hemi+1)*nmap/2]*intc[k];
      }
      //printf("%g ", q1*(1-frac)+q2*frac);
      gsl_matrix_set(a, i, j, q1*(1-frac)+q2*frac);
    }
    //printf("\n");

    gsl_vector_set(b, i, samp[i].q); 
  }

  work=gsl_multifit_linear_alloc(nsamp, nev+cflag);
  x=gsl_vector_alloc(nev+cflag);
  cov=gsl_matrix_alloc(nev+cflag, nev+cflag);

  gsl_multifit_linear(a, b, x, cov, &chisq, work);

  //gsl_vector_fprintf(stdout, x, "%g");

  //output final, interpolated field:
  //if there is a constant term, we only output one field:
  //at the lead time...
  if (cflag) {
    qvec=new float *[1];
    qvec[0]=new float[n];
    for (int32_t k=0; k<n; k++) qvec[0][k]=gsl_vector_get(x, nev);
    for (int32_t i=0; i<nev; i++) {
      for (ind_t k=0; k<n; k++) qvec[0][k]+=gsl_vector_get(x, i)*qall[i][N][k];
    }
  } else {

    qvec=new float *[nall];
    qvec[0]=new float[nall*n];
    for (int32_t k=0; k<n; k++) qvec[0][k]=0;
    for (int32_t i=1; i<nall; i++) {
      qvec[i]=qvec[0]+i*n;
      for (ind_t k=0; k<n; k++) qvec[i][k]=0;
    }

    for (int32_t i=0; i<nev; i++) {
      for (int32_t j=0; j<nall; j++) {
        for (ind_t k=0; k<n; k++) qvec[j][k]+=gsl_vector_get(x, i)*qall[i][j][k];
      }
    }
  }

  for (int32_t i=0; i<nev; i++) {
    delete [] qall[i][0];
    delete [] qall[i];
  }
  delete [] qall;

  delete qdum;
  delete map;
  delete xgrid;
  delete ygrid;

  //delete [] sparse_array;
  //sparse_array=NULL;

  delete [] v[0];
  delete [] v;

  //there probably should be a use for these:
  delete [] eval;

  //we can't forget to delete the gsl data structures:
  gsl_matrix_free(a);
  gsl_vector_free(b);
  gsl_vector_free(x);
  gsl_matrix_free(cov);
  gsl_multifit_linear_free(work);

  delete a_mat;

  return qvec;

}

//1. returns only the coefficients
//2. replaces all the measurements in samp with interpolated values...
float * proxy_tracer(float **tall, int32_t np, time_class *t, int32_t nall, meas_data *samp, long nsamp, int32_t order) {

  gsl_matrix *q1;		//we copy the interpolated tracer values here transformed by the basis functions...
  gsl_vector *q2;		//measured tracer values here...

  gsl_vector *k;		//regression coefficients
  float *k2;			//coefficients as regular vector

  gsl_matrix *cov;
  gsl_multifit_linear_workspace *work;
  double chisq;
  
  q1=gsl_matrix_alloc(nsamp, order+1);
  q2=gsl_vector_alloc(nsamp);

  //copy measured values to a gsl vector:
  for (long i=0; i<nsamp; i++) gsl_vector_set(q2, i, samp[i].q);

  //interpolated tracer to measurement locations:
  tracer_interp(tall, t, nall, np, samp, nsamp);

  //fill the matrix with interpolated tracer values transformed by basis functions:
  for (long i=0; i<nsamp; i++) {
    gsl_matrix_set(q1, i, 0, 1);
    for (int32_t j=1; j<=order; j++) {
      gsl_matrix_set(q1, i, j, gsl_matrix_get(q1, i, j-1)*samp[i].q);
    }
  }

  //perform the regression analysis:
  k=gsl_vector_alloc(order+1);
  work=gsl_multifit_linear_alloc(nsamp, order+1);
  cov=gsl_matrix_alloc(order+1, order+1);
  
  gsl_multifit_linear(q1, q2, k, cov, &chisq, work);

  //copy the result into a regular C++ vector:
  k2=new float[order];
  for (int i=0; i<order; i++) k2[i]=gsl_vector_get(k, i);

  gsl_matrix_free(q1);
  gsl_vector_free(q2);
  gsl_vector_free(k);
  gsl_matrix_free(cov);
  gsl_multifit_linear_free(work);

  return k2;

}

//creates the sparse matrix mapping for tracer simulation:
template <class real>
void tracer_map(ctraj_tfield_base<real> * tracer, 
		ctraj_vfield_base<real> *vfield,
		sparse<int32_t, real> *map,
                real t, real dt, int32_t nt) {
  int32_t n;		//number of grid points
  int32_t ndim;		//number of dimensions:
  int32_t nvert;	//number of vertices:
  //grid point positions:
  real **x;
  int32_t *d;
  //for performing the integration:
  real **result;
  void *param[2];
  //interpolation weights and indices:
  double *wt;
  int32_t *ind;
  int32_t dnew;

  //important parameters:
  n=tracer->nel();
  ndim=tracer->ndim();
  nvert=tracer->nvert();

  //get the grid locations:
  x=new real * [n];
  x[0]=new real[n*ndim];
  d=new int32_t[n];                 //initial domain

  //get the initial conditions:
  for (int32_t i=0; i<n; i++) {
    x[i]=x[0]+i*ndim;
    d[i]=tracer->get_loc(i, x[i]);
  }

  //initialize interpolation variables:
  wt=new double[nvert];
  ind=new int32_t[nvert];

  //initialize the vector of results for the Runge-Kutta integrations:
  result=new float * [nt+1];
  result[0]=new float[ndim*(nt+1)];
  for (int32_t i=1; i<=nt; i++) result[i]=result[0]+i*ndim;

  param[0]=vfield;

  map.reset(n, n);

  for (int32_t i=0; i<n; i++) {
    //set initial conditions:
    param[1]=d+i;

    //do a Runge-Kutta integration:
    rk_dumb(t, x[i], ndim, dt, nt, result, (void *) param, &ctraj_deriv);

    dnew=tracer->interpolate(d[i], result[nt], ind, wt);

    if (dnew < 0) {
      map.add_el(1, i, i);
    } else {
      for (int32_t j=0; j<nvert; j++) {
        map.add_el(wt[j], i, ind[j]);
      }
    }
  }

  delete [] x[0];
  delete [] x;
  delete [] d;
  delete [] result[0];
  delete [] result;
  delete [] wt;
  delete [] ind;

}

  
