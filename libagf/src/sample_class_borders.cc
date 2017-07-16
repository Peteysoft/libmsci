//
// This software is released under the following terms:
//
// 1. No commercial use.
// 2. Copies and derivative works are free to use and modify.
// 3. Attribution must be given to all contributors of both original and derivative works.
//
// Authors:
//
// 2017-07-16 Peter Mills: added license information 
//

//
// Given a function, r(x), and a set of pairs of samples, {(y_i, z_i)}, returns 
// a set of points, {b_j}, and a set of vectors, {v_j}, for which:
//
// r(b_j) = 0
// v_j = grad(r(x)) | x=b_j
// for every j there exists an i s.t. b_j = f y_i + (1 - f) z_i | 0 <= f <= 1
//

#include <stdio.h>
#include <math.h>
#include <stdint.h>
#include <assert.h>

#include "supernewton.h"
#include "agf_lib.h"

#include "tree_lg.h"
#include "peteys_tmpl_lib.h"

using namespace std;
using namespace libpetey;

namespace libagf {

  //a set of program designed to find the class borders
  //and use these to perform classifications...

  //this structure holds the parameters to pass to the 1-D function 
  //to be minimized (section of a line between two points of opposite sign):
  template <class real>
  struct bf_params {
    real *vec1;			//the first vector
    real *v;			//the line between the two
  
    real *vec;		//to avoid having to allocate and reallocate it...
    real *grd;		//  "
    //--also for passing back to the main routine...

    real D;		//same fucking thing as in the other "param" structure

    real rthresh;

    //for calculating the conditional probabilities:
    real (*rfunc) (real *, void *, real *);	//cond. prob. (2-D)
    void *param;			//parameters for calculating the cond. prob.

    int dumflag;

  };

  //**note: these functions are compatible with the GSL...
  template <class real>
  void bfind(real t, 			//line parameter
		void *params1, 		//parameters passed from minim. caller
		real *d1, 		//difference in cond. probabilities
		real *dddt1) { 		//deriv. of diff. in cond. prob. along line
    bf_params<real> *params=(bf_params<real> *) params1;
    real d, dddt;

    /*	//this is numerical derivative calculation for validation...
    real dt=0.001;
    real da, db;
    real dddt0;

    if (params->dumflag) {
      params->dumflag=0;
      bfind(t+dt, params1, &da, &dddt0);
      bfind(t-dt, params1, &db, &dddt0);
      params->dumflag=1;
    }
    */

    //calculate the current position along the line:
    for (dim_ta i=0; i<params->D; i++) {
      params->vec[i]=params->vec1[i]+t*params->v[i];
    }

    //calculate the difference between the two conditional probabilities
    //and the gradient vector:
    d=(* params->rfunc) (params->vec, params->param, params->grd);
    *d1 = d-params->rthresh;

    //calculate the change in d wrt the line parameter:
    dddt=0;

    for (dim_ta i=0; i<params->D; i++) dddt+=params->grd[i]*params->v[i];

    //if (params->dumflag) printf("r=%g; drdt=%g drdt (num.)=%g\n", d, dddt, (da-db)/dt/2);
    *dddt1=dddt;
  }

  //returns a pair of random samples from each class
  //return value is an error code: negative value indicates error
  template <class real>
  int oppositesample(void *param1, 			//paramters: sample locations etc.
			real *x1, 		//sample from first class
			real *x2) {		//sample from second class
    bordparam<real> *param=(bordparam<real> *) param1;
    nel_ta s1, s2;

    do {
      if (param->sind1->nel() >= param->maxsample) {
        //param->sind1->print(stdout);
        //fprintf(stderr, "agfsample: no more classes to sample\n");
        return -1;
      }
      s1=(nel_ta) (gsl_rng_uniform(param->rann)*param->ind);
      s2=(nel_ta) (gsl_rng_uniform(param->rann)*(param->n-param->ind));
      for (dim_ta i=0; i<param->D; i++) {
        x1[i]=param->train[s1][i];
        x2[i]=param->train[s2+param->ind][i];
      }
    } while (param->sind1->add_member((int64_t) s2*(int64_t) param->ind+(int64_t) s1) < 0);

    return 0;
  }

  template <class real>
  int oppositesample_small(void *param1, real *x1, real *x2) {
    bordparam<real> *param=(bordparam<real> *) param1;
    nel_ta s1, s2;
  
    if (param->s >= param->maxsample) {
      fprintf(stderr, "agfsample_small: no more classes to sample\n");
      return -1;
    }

    s1=param->sind2[param->s] % param->n;
    s2=param->sind2[param->s]/param->n;

    //printf("oppositesample_small: s1=%d; s2=%d\n", s1, s2);
    param->s++;
    for (dim_ta i=0; i<param->D; i++) {
      x1[i]=param->train[s1][i];
      x2[i]=param->train[s2][i];
    }

    return 0;
  }

  //initialize stuff for sampling routines primarily:
  template <class real>
  int bordparam_init(bordparam<real> *param,
			real **train,
			dim_ta D,
			nel_ta n,
			nel_ta clind,
			int smallflag) {
    double *ran_vec;

    //not much to this...
    param->D = D;
    param->train = train;
    param->ind=clind;
    param->n=n;
    
    //random number generator:
    param->rann=gsl_rng_alloc(agf_gsl_rng_type);
    gsl_rng_set(param->rann, seed_from_clock());

    if (smallflag) {
      //maximum possible number of samples, random # for each, sort the #s:
      param->maxsample=0;
      param->sind2=new int64_t[n*(n-1)/2];
      for (nel_ta i=0; i<n; i++) {
        for (nel_ta j=0; j<i; j++) {
          param->sind2[param->maxsample]=i*n+j;
	  param->maxsample++;
	}
      }
      //printf("%d %ld\n", n*(n-1)/2, param->maxsample);
      //assert(param->maxsample==n*(n-1)/2);
      ran_vec=new double[param->maxsample];
      for (nel_ta i=0; i<param->maxsample; i++) {
        ran_vec[i]=gsl_rng_uniform(param->rann);
      }
      long *sind3=heapsort(ran_vec, param->maxsample);
      map_vector_inplace(param->sind2, sind3, param->maxsample);
      param->s=0;
      param->sind1=NULL;

      delete [] ran_vec;
      delete [] sind3;
    } else {
      param->sind1=new tree_lg<int64_t>;
      param->sind2=NULL;
      param->maxsample=clind*(n-clind);
    }

    return 0;
  }

  template <class real>
  void bordparam_clean(bordparam<real> *param) {
    if (param->sind1!=NULL) delete param->sind1;
    if (param->sind2!=NULL) delete [] param->sind2;
    gsl_rng_free(param->rann);
  }

  template <class real>
  void zero_bord_diag(bord_diag<real> *diag) {
    //zero diagnostics:
    diag->min_iter=-1;
    diag->max_iter=0;
    diag->total_iter=0;
    diag->min_tol=1;			//maximum allowable
    diag->max_tol=0;
    diag->total_tol=0;
    diag->nfail=0;
    diag->nbis=0;
  }

  template <class real>
  void print_bord_diag(FILE *fs, bord_diag<real> *diag, nel_ta nfound) {
    //print out the diagnostics:
    fprintf(fs, "\n\n");
    fprintf(fs, "diagnostic parameter          %8s   %8s   %8s\n", 
		  "min", "max", "average");
    fprintf(fs, "\n");
    fprintf(fs, "iterations in supernewton:**   %8d   %8d %10.3g\n",  
		  diag->min_iter, diag->max_iter, (float) diag->total_iter/(float) nfound);
    fprintf(fs, "tolerance of samples:          %10.3g %10.3g %10.3g\n",  
		  diag->min_tol, diag->max_tol, diag->total_tol/nfound);
    fprintf(fs, "\n");
    fprintf(fs, "** total number of iterations:     %d\n", diag->total_iter);
    fprintf(fs, "   number of bisection steps:      %d\n", diag->nbis);
    fprintf(fs, "   number of convergence failures: %d\n", diag->nfail);
    fprintf(fs, "\n");
  }

  //train a binary classifier in which the Bayesian borders are discretely sampled:
  template <class real>
  nel_ta single_sample_cb(
		//returns difference in conditional prob. plus derivatives:
		real (*rfunc) (real *, void *, real *),
                //returns a pair of random samples (one for each class):
                int (*sample) (void *, real *, real *),
                void *param,		//these are just along for the ride
                dim_ta D,		//number of dimensions
                real tol,		//desired tolerance
		iter_ta maxit,		//max no. of iterations in root-finder
                real *border,		//returned border samples
                real *gradient,		//return border gradients
		bord_diag<real> *diag,	//diagnostics
                real rthresh)		//location of Bayesian border
  {
    real d1, d2;		//difference between conditional probabilities
  				//of two classes
    real *grad1, *grad2;	//gradients of initial brackets
    real dddt1, dddt2;		//derivatives at initial brackets

    real t0, t2;		//line parameters: at the root, at the second vector
    supernewton_stat err;	//error code returned from root-finding routine
    bf_params<real> bfparam;	//structure of parameters to pass function to minimize

    void (* funcd) (real, void *, real *, real *);

    gsl_error_handler_t * old_handler;

    funcd=&bfind<real>;

    //when you strip it down to its essentials, there really isn't much to this algorithm:

    //allocate space:
    bfparam.vec1=new real[D];
    bfparam.v=new real[D];

    bfparam.vec=new real[D];
    bfparam.grd=new real[D];

    //parameters we need in the function to be zeroed:
    bfparam.dumflag=1;
    bfparam.D=D;
    bfparam.rthresh=rthresh;
    bfparam.rfunc=rfunc;

    /*FLAG*/
    //I think this is pretty clear:
    //three levels of voided "param" structures....   AWESOME
    //(really need to fix this or at least clean up the naming...)
    bfparam.param=param;

    //allocate more space:
    grad1=new real[D];
    grad2=new real[D];

    //change the GSL error handler so that an interpolation failure
    //in supernewton will not interrupt execution:
    old_handler=gsl_set_error_handler(&agf_gsl_handler);

    //find a single sample of the class borders:
    do {
      if ((*sample) (param, bfparam.vec1, bfparam.vec)<0) {
        err.code=OUT_OF_DATA;
        goto finish;
      }
      d1=(*rfunc) (bfparam.vec1, param, grad1)-rthresh;
      d2=(*rfunc) (bfparam.vec, param, grad2)-rthresh;
      //printf("d1=%g; d2=%g\n", d1, d2);
    } while (d1*d2>=0);

    //find the parametric representation of the line between the two vectors:
    for (dim_ta j=0; j<D; j++) {
      bfparam.v[j]=bfparam.vec[j]-bfparam.vec1[j];
    }
    t2=1;

    //to avoid redundant function calls, we pass the conditional prob.
    //and their derivatives to the root-finding routine:
    dddt1=0;
    dddt2=0;
    for (dim_ta j=0; j<D; j++) {
      dddt1+=bfparam.v[j]*grad1[j];
      dddt2+=bfparam.v[j]*grad2[j];
    }

    //find the root of the difference between the 
    //conditional probabilities along
    //the line between the two vectors in order to find the class border:
    t0=supernewton(funcd, (void *)&bfparam, (real) 0., t2, tol, (real) 0., 
		maxit, &err, d1, dddt1, d2, dddt2);
    if (err.code != 0) {
      diag->nfail++;
      goto finish;
    }

    if (err.niter < diag->min_iter) diag->min_iter=err.niter;
    else if (err.niter > diag->max_iter) diag->max_iter=err.niter;
    diag->total_iter+=err.niter;
    diag->nbis+=err.nbis;

    //the class border and the gradient vector
    //should be sitting in the "bfparam" structure:
    for (dim_ta j=0; j<D; j++) {
      border[j]=bfparam.vec[j];
      gradient[j]=bfparam.grd[j];
    }
    //printf("d1=%f\n", d1);
    //for (long k=0; k<D; k++) printf("%f ", gradient[i][k]);
    //printf("\n");
    
    //calculate more diagnostics:
    d1=fabs(d1);
    if (d1 < diag->min_tol) diag->min_tol=d1;
    else if (d1 > diag->max_tol) diag->max_tol=d1;
    diag->total_tol+=d1;

    finish:
      //delete integer and floating point arrays:
      delete[] bfparam.vec1;
      delete[] bfparam.v;
      delete [] bfparam.vec;
      delete [] bfparam.grd;
      delete [] grad1;
      delete [] grad2;

      //set error handler back to previous one:
      gsl_set_error_handler(old_handler);

    return err.code;

  }

  //train a binary classifier in which the Bayesian borders are discretely sampled:
  template <class real>
  nel_ta sample_class_borders(
		//returns difference in conditional prob. plus derivatives:
		real (*rfunc) (real *, void *, real *),
                //returns a pair of random samples (one for each class):
                int (*sample) (void *, real *, real *),
                void *param,		//these are just along for the ride
                nel_ta n,		//number of time to sample
                dim_ta D,		//number of dimensions
                real tol,		//desired tolerance
		iter_ta maxit,		//max no. of iterations in root-finder
                real **border,		//returned border samples
                real **gradient,	//return border gradients
                real rthresh)		//location of Bayesian border
  {
    bord_diag<real> diag;
    nel_ta nfound=0;
    int err;

    zero_bord_diag<real>(&diag);
    printf("%7d of %7d vectors found: %5.1f%%", 0, n, 0.);
    for (nel_ta i=0; i<n; i++) {
      for (int j=0; j<40; j++) printf("\b");
      printf("%7d of %7d vectors found: %5.1f%%", i+1, n, (100.*(i+1))/n);
      fflush(stdout);
      err=single_sample_cb(rfunc, sample, param, D, tol, maxit, border[i], gradient[i], &diag, rthresh);
      if (err==0) {
        nfound++;
      } else if (err==OUT_OF_DATA) {
        fprintf(stderr, "sample_class_borders: ran out of samples (%d found), returning\n", nfound);
	break;
      } else {
        i--;
      }
    }
    print_bord_diag<real>(stdout, &diag, nfound);
    return nfound;
  }

  //function that uses borders and border gradients to perform a
  //classification:
  template <class real>
  real border_classify0(real **brd, 		//border samples
			real **grd, 		//gradient vectors
			dim_ta D, 		//number of dimensions
			nel_ta n, 		//number of samples
			real *x, 		//test point
			nel_ta &k,		//index of nearest sample
			real &d2m) {		//distance to sample

    real d2;		//distance squared
    real p;		//dot product

    //find the nearest border sample:
    d2m=metric2(x, brd[0], D);
    k=0;
    for (nel_ta i=1; i<n; i++) {
      d2=metric2(x, brd[i], D);
      if (d2 < d2m) {
        d2m=d2;
        k=i;
      }
    }

    //calculate the dot product of the difference with the gradient
    //--the sign is the class
    p=0;
    for (dim_ta j=0; j<D; j++) {
      p+=grd[k][j]*(x[j]-brd[k][j]);
    }
    //printf("d1=%g; d2=%g\n", sqrt(d2m), p/sqrt(dg));

    return p;

  }

  template <class real>
  int test_oppositesample(nel_ta n1, nel_ta n2) {
    bordparam<real> param;
    int maxsample=n1*n2;
    dim_ta D=1;
    real *data[n1+n2];
    int found[maxsample];
    int cnt;
    real sample1, sample2;
    int err=0;
    real sep;

    data[0]=new real [n1+n2];
    for (nel_ta i=0; i<n1+n2; i++) {
      data[i]=data[0]+i;
      data[i][0]=i;
    }
    for (int i=0; i<maxsample; i++) {
      found[i]=-1;
    }

    bordparam_init(&param, data, D, n1+n2, n1, 0);

    cnt=0;
    while (oppositesample((void *) &param, &sample1, &sample2)>=0) {
      long ind1, ind2, k;
      ind1=bin_search(data[0], n1, sample1, -1);
      if (ind1<0 || ind1>=n1) {
        fprintf(stderr, "test_oppositesample: 0 class; returned sample not in dataset.\n");
	err=1;
      }
      if (err) break;
      ind2=bin_search(data[0]+n1, n2, sample2, -1);
      if (ind2<0 || ind2>=n2) {
        fprintf(stderr, "test_oppositesample: 1 class; returned sample not in dataset.\n");
	err=1;
      }
      if (err) break;

      k=n1*ind2+ind1;

      if (found[k]>=0) {
        fprintf(stderr, "test_oppositesample: generated duplicate sample pair, (%d, %d).\n", ind1, ind2);
	err=1;
	break;
      }
      found[k]=cnt;
      cnt++;
    }

    delete [] data[0];
    bordparam_clean(&param);

    if (err!=0) return err;

    if (cnt!=maxsample) {
      fprintf(stderr, "test_oppositesample: oppositesample didn't return all possible permutations (%d vs. %d)\n", cnt, maxsample);
      err=1;
    }

    sep=0;
    for (int i=0; i<maxsample; i++) {
      if (found[i]<0) {
        fprintf(stderr, "test_oppositesample: missing sample pair, (%d, %d).\n",i/n1, i%n1);
	err=1;
      }
      sep+=fabs(found[i]-i);
    }

    printf("test_oppositesample(%3d, %3d):  average separation out of %6d = %8.1f\n", n1, n2, maxsample, sep/n1/n2);

    return err;
  }

  template <class real>
  int test_oppositesample_small(nel_ta n) {
    int maxsample=n*(n-1)/2;
    bordparam<real> param;
    dim_ta D=1;
    real *data[n];
    int found[maxsample];
    int cnt;
    real sample1, sample2;
    int err=0;
    real sep;

    data[0]=new real [n];
    for (nel_ta i=0; i<n; i++) {
      data[i]=data[0]+i;
      data[i][0]=i;
    }
    for (int i=0; i<maxsample; i++) {
      found[i]=-1;
    }

    //fifth parameter shouldn't matter...
    bordparam_init(&param, data, D, n, n/2, 1);

    cnt=0;
    while (oppositesample_small((void *) &param, &sample1, &sample2)>=0) {
      long ind1, ind2, k;
      ind1=bin_search(data[0], n, sample1, -1);
      if (ind1<0 || ind1>=n) {
        fprintf(stderr, "test_oppositesample_small: 0 class; returned sample not in dataset.\n");
	err=1;
      }
      if (err) goto fail;
      ind2=bin_search(data[0], n, sample2, -1);
      if (ind2<0 || ind2>=n) {
        fprintf(stderr, "test_oppositesample_small: 1 class; returned sample not in dataset.\n");
	err=1;
      }
      if (err) goto fail;

      //assumes ind2 is the larger of the two:
      if (ind1>ind2) {
        //(but we can always swap them...)
        nel_ta swp=ind1;
	ind1=ind2;
	ind2=swp;
      } else if (ind1==ind2) {
        fprintf(stderr, "test_oppositesample_small: samples are the same (ind1==ind2)\n");
	err=1;
	goto fail;
      }
      k=(ind2)*(ind2-1)/2+ind1;

      if (k<0 || k>=maxsample) {
        fprintf(stderr, "test_oppositesample_small: something's wrong here (k=%d; maxsample=%d)\n", k, maxsample);
	err=1;
	goto fail;
      }

      //printf("k=%d\n", k);

      if (found[k]>=0) {
        fprintf(stderr, "test_oppositesample_small: generated duplicate sample pair, (%d, %d); k=%d.\n", ind1, ind2, k);
	err=1;
	break;
      }
      found[k]=cnt;
      cnt++;
    }

    if (cnt!=maxsample) {
      fprintf(stderr, "test_oppositesample_small: oppositesample_small didn't return all possible permutations (%d vs. %d)\n", cnt, maxsample);
      err=1;
    }

    for (int i=0; i<maxsample; i++) {
      if (found[i]<0) {
        fprintf(stderr, "test_oppositesample_small: missing sample pair, number %d", i);
	err=1;
	goto fail;
      }
      sep+=fabs(found[i]-i);
    }

    printf("test_oppositesample_small(%3d): average separation out of %6d = %8.1f\n", n, maxsample, 2*sep/n/(n-1));

    fail:
      delete [] data[0];
      bordparam_clean(&param);

    return err;
  }

  template int bordparam_init<float>(bordparam<float> *, float **, dim_ta, nel_ta, nel_ta, int); 
  template int bordparam_init<double>(bordparam<double> *, double **, dim_ta, nel_ta, nel_ta, int); 

  template void bordparam_clean<float>(bordparam<float> *);
  template void bordparam_clean<double>(bordparam<double> *);

  template void zero_bord_diag<float>(bord_diag<float> *);
  template void zero_bord_diag<double>(bord_diag<double> *);

  template void print_bord_diag<float>(FILE *, bord_diag<float> *diag, nel_ta);
  template void print_bord_diag<double>(FILE *, bord_diag<double> *diag, nel_ta);

  template int oppositesample<float>(void *, float *, float *);
  template int oppositesample<double>(void *, double *, double *);

  template int oppositesample_small<float>(void *, float *, float *);
  template int oppositesample_small<double>(void *, double *, double *);

  template int single_sample_cb<float>
		(float (*) (float *, void *, float *), 
		int (*) (void *, float *, float *), 
                void *, dim_ta, float, iter_ta,
                float *, float *, 
		bord_diag<float> *, float);

  template int single_sample_cb<double> 
		(double (*) (double *, void *, double *), 
                int (*) (void *, double *, double *), 
                void *,  dim_ta, double, iter_ta,
                double *, double *, 
		bord_diag<double> *, double);

  template nel_ta sample_class_borders<float>
		(float (*) (float *, void *, float *), 
		int (*) (void *, float *, float *), 
                void *, nel_ta, dim_ta, float, iter_ta, 
                float **, float **, float);

  template nel_ta sample_class_borders<double> 
		(double (*) (double *, void *, double *), 
                int (*) (void *, double *, double *), 
                void *,  nel_ta, dim_ta, double, iter_ta,
                double **, double **, double);

  template float border_classify0<float>(float **, float **, 
			dim_ta, nel_ta, float *, nel_ta &, float &);

  template double border_classify0<double>(double **, double **, 
			dim_ta, nel_ta, double *, nel_ta &, double &);

  //unit tests:
  template int test_oppositesample<float>(nel_ta, nel_ta);
  template int test_oppositesample<double>(nel_ta, nel_ta);

  template int test_oppositesample_small<float>(nel_ta);
  template int test_oppositesample_small<double>(nel_ta);

} //end namespace libagf

