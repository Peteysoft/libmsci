#include <stdint.h>
#include <stdio.h>
#include <math.h>

#include <gsl/gsl_errno.h>

#include "quicksort.h"
#include "supernewton.h"

#include "agf_util.h"
#include "balloon_param.h"

using namespace libpetey;

namespace libagf {

  template <typename real>
  real gaussian_kernel(real r2, real *dfdr) {
    real f;
    f=exp(-r2/2);
    if (dfdr != NULL) dfdr=-f/2;
  }

  template <typename real, typename index_t>
  balloon_param<real, index_t>::balloon_param() {
  }

  template <typename real, typename index_t>
  balloon_param<real, index_t>::balloon_param(
		real (*K) (real, real *), 		//kernel function
		index_t nvar, 				//number of features
		real total_weight,			//total of weights
		index_t k2,
		real v1,
		real v2,
		real tolerance,
		index_t maxN) {
    kernel=K;
    D=nvar;
    W=total_weight;
    k=k2;

    var[0]=v1;
    var[1]=v2;

    tol=tolerance;
    maxiter=maxN;
  }

  template <typename real, typename index_t>
  balloon_param<real, index_t>::~balloon_param() {
  }

  template <typename real, typename index_t>
  void balloon_param<real, index_t>::maxn(index_t *mn) {
    if (k > 0) *mn=k;
  }

  template <typename real, typename index_t>
  index_t balloon_param<real, index_t>::dim() {
    return D;
  }

  template <typename real, typename index_t>
  struct balloon_werr_param {
    real W;
    index_t n;
    real *d2;
    real *w;
    real (* kernel) (real, real *);
  };

  template <typename real, typename index_t>
  void balloon_werr(real var, void *param, real *wdiff, real *dWdv) {
    balloon_werr_param<real, index_t> *p2=(balloon_werr_param<real, index_t> *) p2;
    real W;
    real deriv;

    W=0;
    *dWdv=0;
    for (index_t i=0; i<p2->n; i++) {
      p2->w[i]=(* p2->kernel) (p2->d2[i]/var, &deriv);
      (*dWdv)-=p2->d2[i]*deriv;
      W+=p2->w[i];
    }
    (*dWdv)/=var*var;

    *wdiff=W-p2->W;
  }


  template <typename real, typename index_t>
  index_t balloon_param<real, index_t>::operator () (
		real *d2, 				//distances
		index_t n,				//number of distances
		index_t *ind,				//returned indices
		real *w,				//returned weights
		real *norm) {				//returned norm. coef.
    real *kleast=NULL;
    supernewton_stat err;
    index_t nmov;          //number of times bracket has been moved
                        //(count as extra iteration...)
    real werr1, werr2, dWdv1, dWdv2;
    gsl_error_handler_t * old_handler;
    real vtest;
    //real var1=var[0];
    //real var2=var[1];
    balloon_werr_param<real, index_t> param;
    real var_f;

    if (k > 0 && k < n) {
      kleast=new real[k];
      kleast_quick(d2, n, k, kleast, ind);
      d2=kleast;
      n=k;
    } else {
      for (index_t i=0; i<n; i++) ind[i]=i;
    }

    param.W=W;
    param.n=n;
    param.d2=d2;
    param.w=w;
    param.kernel=kernel;

    for (nmov=0; nmov < maxiter; nmov++) {
      balloon_werr<real, index_t>(var[0], (void *) &param, &werr1, &dWdv1);
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
    for ( ; nmov < maxiter; nmov++) {
      balloon_werr<real, index_t>(var[1], (void *) &param, &werr2, &dWdv2);
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
      fprintf(stderr, "             W_1(var_1=%g)=%f; W_2(var_2=%g)=%f\n", werr1+W, var[0], werr2+W, var[1]);
      return -1;
    }

    old_handler=gsl_set_error_handler(&agf_gsl_handler);

    var_f=supernewton(&balloon_werr<real, index_t>, (void *) &param, var[0], var[1], (real) 0., tol, maxiter, &err, werr1, dWdv1, werr2, dWdv2);

    //set error handler back to previous one:
    gsl_set_error_handler(old_handler);

    //if required, calculate normalization coefficient:
    if (norm != NULL) *norm=pow(sqrt(var_f*M_PI*2), D);

    delete [] kleast;

    if (k>0 && k<n) return k; else return n;

  }

  template class balloon_param<float, long>;
  template class balloon_param<double, long>;


} //end namespace libagf

