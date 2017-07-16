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
// Central "engine" of adaptive Gaussian filtering: set filter width (variance)
// to keep the total of the weights constant.
//
// See:
//
// Peter Mills (2011). "Efficient statistical classification of satellite 
// measurements." International Journal of Remote Sensing 32 (21): 6109-6132.
//

#include <stdio.h>
#include <math.h>

#include <gsl/gsl_errno.h>
#include "supernewton.h"

#include "agf_lib.h"

using namespace std;
using namespace libpetey;

namespace libagf {

  long agf_global_weights_maxiter=WEIGHTS_MAXITER;
  real_a agf_global_weights_tol=WEIGHTS_TOL;

  //given a set of distances,
  //calculates the weights of an "adaptive Gaussian filter"

  template <class real>
  iter_ta agf_calc_w(real *d2,     //distances squared
              nel_ta k,           //number of distances
              real Wc,         //objective total weight
	      real var[2],	//initial filter width
              real *weight,	//returned weights
              real &var_f)         //returned final filter variance
  {

    iter_ta nd;                      //number of squarings of the weights
    real tw_old;                 //previous value of the total weight
    real tw;                     //current value of the total weight

    real wtint;                  //for interpolation of final weights

    //calculate the weights:
    for (nel_ta i=0; i<k; i++) weight[i]=exp(-d2[i]/var[1]/2);

    //repeatedly square the weights until the total is less than the threshold:
    tw=0;
    for (nel_ta i=0; i<k; i++) tw+=weight[i];
    nd=0;

    do {
      tw_old=tw;
      tw=0;
      for (nel_ta i=0; i<k; i++) {
        weight[i]*=weight[i];
        tw+=weight[i];
      }
      nd++;
    } while (tw > Wc && tw < tw_old);

    //interpolate between the current and previous weight:
    wtint=(log(tw_old)-log(Wc))/(log(tw_old)-log(tw))/2+0.5;

    //calculate the final filter width:
    var_f=var[1]/(1 << nd)/wtint;

    tw=0;
    for (nel_ta i=0; i<k; i++) {
      //weight[i]=pow(weight[i], wtint);
      //instead of interpolatting, why don't we just calculate the weights directly?
      //(the differences are insignificant, but still...)
      weight[i]=exp(-d2[i]/var_f/2);
    }

    //return the number of iterations as a diagnostic parameter:
    return nd;

  }

  template <class real>
  void agf_werr(real var, void *param, real *wdiff, real *dWdv) {
    void **param2=(void **) param;
    real W0;
    real W;
    nel_ta k;
    real *d2;
    real *w;

    W0=*(real *) param2[0];
    k=*(nel_ta *) param2[1];
    d2=(real *) param2[2];
    w=(real *) param2[3];

    W=0;
    *dWdv=0;
    for (nel_ta i=0; i<k; i++) {
      w[i]=exp(-d2[i]/var/2);
      (*dWdv)+=d2[i]*w[i];
      W+=w[i];
    }
    (*dWdv)/=var*var*2;

    *wdiff=W-W0;
  }

  template <class real>
  iter_ta agf_calc_w2(real *d2,		//distances squared
              nel_ta k,			//number of distances
              real Wc, 			//objective total weight
	      real var[2],		//initial filter width
              real *weight,		//returned weights
              real &var_f)		//returned final filter variance
  {
    supernewton_stat err;
    long nmov;		//number of times bracket has been moved
			//(count as extra iteration...)
    real werr1, werr2, dWdv1, dWdv2;
    gsl_error_handler_t * old_handler;
    real vtest;
    //real var1=var[0];
    //real var2=var[1];
    void *param[4];
    param[0]=&Wc;
    param[1]=&k;
    param[2]=d2;
    param[3]=weight;

    for (nmov=0; nmov < WEIGHTS_MAXITER; nmov++) {
      agf_werr<real>(var[0], param, &werr1, &dWdv1);
      if (werr1 <= 0) break;
      //make it "sticky":
      vtest=var[0]/2;
      //check for underflow:
      if (isnormal(vtest)) {
        var[0]=vtest;
      } else {
        fprintf(stderr, "agf_calc_w2: lower bracket underflow\n");
	var_f=var[0];
        break;
      }
      fprintf(stderr, "agf_calc_w2: lower bracket decreased to %g\n", var[0]);
    }
    for ( ; nmov < WEIGHTS_MAXITER; nmov++) {
      agf_werr<real>(var[1], param, &werr2, &dWdv2);
      if (werr2 >= 0) break;
      //make it "sticky":
      vtest=var[1]*2;
      //check for overflow:
      if (isnormal(vtest)) {
        var[1]=vtest;
      } else {
        fprintf(stderr, "agf_calc_w2: upper bracket overflow\n");
	var_f=var[1];
        break;
      }
      fprintf(stderr, "agf_calc_w2: upper bracket increased to %g\n", var[1]);
    }

    if (werr1 > 0 || werr2 < 0) {
      fprintf(stderr, "agf_calc_w2: failed to bracket root\n");
      fprintf(stderr, "             W_1(var_1=%g)=%f; W_2(var_2=%g)=%f\n", werr1+Wc, var[0], werr2+Wc, var[1]);
      return -1;
    }

    old_handler=gsl_set_error_handler(&agf_gsl_handler);

    var_f=supernewton(&agf_werr<real>, (void *) param, var[0], var[1], (real) 0., (real) agf_global_weights_tol, agf_global_weights_maxiter, &err, werr1, dWdv1, werr2, dWdv2);

    //set error handler back to previous one:
    gsl_set_error_handler(old_handler);

    return err.niter+nmov;

  }


  //calculate the gradient of the weights:
  template <class real>
  void agf_grad_w(real **x,		//matrix of samples
		dim_ta n,		//# of dimensions
		real *xvec,		//test point
		real *w,		//weights
		real *d2,		//squared distances
		nel_ta k,		//number of weights and indices
		real var_f,		//filter width
		real **dwdx)  		//returned gradient vectors
  {
    real dwcnum, dwcden;

    for (dim_ta j=0; j<n; j++) {
      dwcnum=0;
      dwcden=0;
      for (nel_ta i=0; i<k; i++) {
        dwcnum+=w[i]*(x[i][j]-xvec[j]);
        dwcden+=d2[i]*w[i];
      }
      for (nel_ta i=0; i<k; i++) {
        dwdx[i][j]=w[i]*(x[i][j]-xvec[j]-d2[i]*dwcnum/dwcden)/var_f;
      }
    }
  }

  template iter_ta agf_calc_w<float>(float *d2, nel_ta k, float Wc, 
		float var[2], float *weight, float &var_f);
  template iter_ta agf_calc_w2<float>(float *d2, nel_ta k, float Wc, 
		float var[2], float *weight, float &var_f);
  template void agf_grad_w<float>(float **x, dim_ta n, float *xvec, float *w,	
		float *d2, nel_ta k, float var_f, float **dwdx);

  template iter_ta agf_calc_w<double>(double *d2, nel_ta k, double Wc, 
		double var[2], double *weight, double &var_f);
  template iter_ta agf_calc_w2<double>(double *d2, nel_ta k, double Wc, 
		double var[2], double *weight, double &var_f);
  template void agf_grad_w<double>(double **x, dim_ta n, double *xvec, 
		double *w, double *d2, nel_ta k, double var_f, double **dwdx);

} //end namespace libagf
