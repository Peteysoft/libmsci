#include <math.h>

#include "quicksort.h"

  template <typename real>
  real agf_kernel(real r2, real *dfdr) {
    real f;
    f=exp(-r2/2);
    if (dfdr != NULL) dfdr=-f/2;
  }

  template <typename real, typename index_t>
  class agf_param {
    private:
      //kernel function:
      real (* kernel) (real, real *);		//kernel function
      index_t k;				//number of weight to calc.
      real var[2];				//variance brackets
      index_t maxiter;
      real tol;
    protected:
      index_t D;				//number of features
      real W;					//total weight
    public:
      agf_param(real (*) (real), index_t, real, index_t);
      ~agf_param();

      virtual void n(index_t *);	//max size needed for weights
      virtual index_t dim();		//number of features

      //compute the weights, not including derivs:
      virtual index_t operator () (real *d2, index_t n, index_t *ind, real *w,
		      real *norm=NULL);
  };


  template <typename real, typename index_t>
  agf_param<real, index_t>::agf_param(
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
  agf_param<real, index_t>::~agf_param() {
  }

  template <typename real, typename index_t>
  void agf_param<real, index_t>::maxn(index_t *mn) {
    if (k > 0) *mn=k;
  }

  template <typename real, typename index_t>
  void agf_param<real, index_t>::dim() {
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
  index_t agf_param<real, index_t>::operator () (
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
    }

    param[0]=&W;
    param[1]=&n;
    param[2]=d2;
    param[3]=weight;
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

    //calculate final weights:
    for (index_t i=0; i<n; i++) w[i]=(* kernel) (d2[i]/var_f, NULL);

    //if required, calculate normalization coefficient:
    if (norm != NULL) *norm=pow(sqrt(var_f*M_PI*2), D);

    delete [] ind;
    delete [] w;
    delete [] kleast;

    if (k>0 && k<n) return k; else return n;

  }

  template <typename real, typename cls_t, typename index_t>
  cls_t agf_classify(agf_param *wt_calc,
		  index_t n,
		  real **x,
		  cls_t *cls,
		  cls_t ncls,
		  real *test,
		  real *p) {
    index_t D=wt_calc->dim();
    real d2[n];
    index_t maxn=n;

    index_t *ind;
    real *w;
    index_t k;
    real tw;

    for (index_t i=0; i<n; i++) {
      d2[i]=0;
      for (index_t j=0; j<D; j++) {
        real diff=x[i][j]-test[j];
	d2[i]+=diff*diff;
      }
    }

    wt_calc->maxn(&maxn);
    ind=new index_t[maxn];
    w=new index_t[maxn];

    k=wt_calc(d2, n, ind, w);

    for (index_t i=0; i<ncls; i++) p[i]=0;
    tw=0;
    for (index_t i=0; i<k; i++) {
      tw+=w[i];
      p[cls[ind[i]]]+=w[i];
    }

    for (index_t i=0; i<ncls; i++) p[i]/=tw;

    delete [] ind;
    delete [] wt;

    return choose_class(pdf, ncls);
  }

  template <typename real, typename index_t>
  class knn_param {
    private:
      index_t k;				//number of weight to calc.
    public:
      knn_param(real (*) (real), index_t, real, index_t);
      ~knn_param();

      virtual void n(index_t *);	//max size needed for weights

      //compute the weights, not including derivs:
      virtual index_t operator () (real *d2, index_t n, index_t *ind, real *w,
		      real *norm=NULL);
  };

  template <typename real, typename index_t>
  knn_param<real, index_t>::knn_param(index_t kay, index_t nvar) {
    k=kay;
    this->W=kay;
    this->D=nvar;
  }

  template <typename real, typename index_t>
  knn_param<real, index_t>::~knn_param() {
  }

  template <typename real, typename index_t>
  void knn_param<real, index_t>::maxn(index_t *mn) {
    *mn=k+1;
  }

  template <typename real, typename index_t>
  index_t knn_param<real, index_t>::operator () (
		real *d2, 				//distances
		index_t n,				//number of distances
		index_t *ind,				//returned indices
		real *w,				//returned weights
		real *norm) {				//returned norm. coef.
    real kleast[k+1];
    index_t D=this->D;
    real V;
    index_t extra=0;

    if (norm != NULL) extra=1;

    kleast_quick(d2, n, k+extra, kleast, ind);
    for (index_t i=0; i<k+extra; i++) w[i]=1;
    if (norm != NULL) {
      real sqrtpir=sqrt(M_PI)*(kleast[k-1]+kleast[k])/2;
      V=1.
      for (index_t i=0; i<D; i++) V*=sqrtpir;
      *norm=V/gsl_sf_gamma(D/2.+1);
    }
    return k;
  }

