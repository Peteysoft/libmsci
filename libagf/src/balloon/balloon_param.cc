#include <stdint.h>
#include <math.h>

#include "quicksort.h"

#include "balloon_param.h"

using namespace libpetey;

  template <typename real>
  real gaussian_kernel(real r2, real *dfdr) {
    real f;
    f=exp(-r2/2);
    if (dfdr != NULL) dfdr=-f/2;
  }

  template <typename real, typename index_t>
  balloon_param<real, index_t>::balloon_param(
		real (*K) (real, real *), 		//kernel function
		index_t nvar, 				//number of features
		real total_weight,			//total of weights
		real k2,
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

    tol=1e-3;
    maxiter=1000;
  }

  template <typename real, typename index_t>
  balloon_param<real, index_t>::~balloon_param() {
  }

  template <typename real, typename index_t>
  void balloon_param<real, index_t>::maxn(index_t *mn) {
    if (k > 0) *mn=k;
  }

  template <typename real, typename index_t>
  void balloon_param<real, index_t>::dim() {
    return D;
  }

  template <typename real, typename index_t>
  void agf_werr(real var, void *param, real *wdiff, real *dWdv) {
    void **param2=(void **) param;
    real W0;
    real W;
    real deriv;
    index_t k;
    real *d2;
    real *w;
    real (*kernel) (real, real *);

    W0=*(real *) param2[0];
    k=*(nel_ta *) param2[1];
    d2=(real *) param2[2];
    w=(real *) param2[3];
    kernel=(* real) (real, real *) param2[4];

    W=0;
    *dWdv=0;
    for (index_t i=0; i<k; i++) {
      w[i]=(* kernel) (d2[i]/var, &deriv);
      (*dWdv)-=d2[i]*deriv;
      W+=w[i];
    }
    (*dWdv)/=var*var;

    *wdiff=W-W0;
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
    void *param[5];

    if (k > 0 && k < n) {
      kleast=new real[k];
      kleast_quick(d2, n, k, kleast, ind);
      d2=kleast;
      n=k;
    } else {
      for (index_t i=0; i<n; i++) ind[i]=i;
    }

    param[0]=&W;
    param[1]=&n;
    param[2]=d2;
    param[3]=w;
    param[4]=kernel;

    for (nmov=0; nmov < maxiter; nmov++) {
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
    for ( ; nmov < maxiter; nmov++) {
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

    var_f=supernewton(&agf_werr<real>, (void *) param, var[0], var[1], (real) 0., tol, maxiter, &err, werr1, dWdv1, werr2, dWdv2);

    //set error handler back to previous one:
    gsl_set_error_handler(old_handler);

    //if required, calculate normalization coefficient:
    if (norm != NULL) *norm=pow(sqrt(var_f*M_PI*2), D);

    delete [] kleast;

    if (k>0 && k<n) return k; else return n;

  }

  template class balloon_param<float, int32_t>;
  template class balloon_param<double, int32_t>;


} //end namespace libagf

