#include <math.h>

#include "agf_balloon.h"

namespace libagf {
  //classify a test point using a variable bandwidth kernel density
  //estimator with constant weight:
  template <typename real, typename cls_t>
  cls_t balloon_classify(
		  balloon_param<real, long> *wt_calc,//classifier parameters
		  nel_ta n,			//number of training samples
		  real **x,			//coordinates
		  cls_t ncls,			//number of classes
		  cls_t *cls,			//discrete ordinates [0-ncls]
		  real *test,			//test point
		  real *p) {			//returned cond. prob.

    dim_ta D=wt_calc->dim();		//number of dim.
    real d2[n];				//squared distances
    long maxn=n;			//maximum number of weights

    long *ind;				//indices of samples in calc.
    real *w;				//weights
    nel_ta k;				//number of samples in calc.
    real tw;				//total of weights

    //calculate distances:
    for (nel_ta i=0; i<n; i++) {
      d2[i]=0;
      for (dim_ta j=0; j<D; j++) {
        real diff=x[i][j]-test[j];
	d2[i]+=diff*diff;
      }
    }

    //allocated space for indices and weights:
    wt_calc->maxn(&maxn);
    ind=new long[maxn];
    w=new real[maxn];

    //use parameter object to calculate weights:
    k=(* wt_calc)(d2, n, ind, w);

    //sum weights to find probabilities:
    for (cls_t i=0; i<ncls; i++) p[i]=0;
    tw=0;
    for (nel_ta i=0; i<k; i++) {
      tw+=w[i];
      p[cls[ind[i]]]+=w[i];
    }
    //normalize probabilities:
    for (cls_t i=0; i<ncls; i++) p[i]/=tw;

    //clean up:
    delete [] ind;
    delete [] w;

    //choose best class:
    return choose_class(p, ncls);
  }

  template <typename real, typename cls_t>
  cls_t balloon_classify_joint(balloon_param<real, long> *wt_calc,
		  nel_ta n,
		  real **x,
		  cls_t ncls,
		  cls_t *cls,
		  real *test,
		  real *p) {		//returned joint probabilities

    dim_ta D=wt_calc->dim();		//number of dimensions
    real d2[n];				//squared distances
    long maxn=n;			//maximum number of weights
    real norm;				//normalization coef.

    long *ind;			//index of points corr. to wts.
    real *w;				//weights
    nel_ta k;				//number of weights

    for (nel_ta i=0; i<n; i++) {
      d2[i]=0;
      for (dim_ta j=0; j<D; j++) {
        real diff=x[i][j]-test[j];
	d2[i]+=diff*diff;
      }
    }

    wt_calc->maxn(&maxn);
    ind=new long[maxn];
    w=new real[maxn];

    k=(*wt_calc)(d2, n, ind, w, &norm);

    for (cls_t i=0; i<ncls; i++) p[i]=0;
    for (nel_ta i=0; i<k; i++) {
      p[cls[ind[i]]]+=w[i];
    }

    for (cls_t i=0; i<ncls; i++) p[i]/=norm/n;

    delete [] ind;
    delete [] w;

    return choose_class(p, ncls);
  }

  //variable bandwidth kernel density estimator with constant weight:
  template <typename real>
  real balloon_density(balloon_param<real, long> *wt_calc,
		  nel_ta n,
		  real **x,
		  real *test) {
    dim_ta D=wt_calc->dim();
    real d2[n];
    real norm;
    long maxn=n;

    long *ind;
    real *w;
    nel_ta k;
    real tw;

    for (nel_ta i=0; i<n; i++) {
      d2[i]=0;
      for (dim_ta j=0; j<D; j++) {
        real diff=x[i][j]-test[j];
	d2[i]+=diff*diff;
      }
    }

    wt_calc->maxn(&maxn);
    ind=new long[maxn];
    w=new real[maxn];

    k=(*wt_calc)(d2, n, ind, w, &norm);

    tw=0;
    for (nel_ta i=0; i<k; i++) tw+=w[i];

    delete [] ind;
    delete [] w;

    return tw/norm/n;
  }

  template <typename real>
  real balloon_interpol(balloon_param<real, long> *wt_calc,
		  nel_ta n,
		  real **x,
		  real *y,
		  real *test,
		  real *e2) {
    dim_ta D=wt_calc->dim();
    real d2[n];
    long maxn=n;

    long *ind;
    real *w;
    nel_ta k;
    real tw;
    real result;

    for (nel_ta i=0; i<n; i++) {
      d2[i]=0;
      for (dim_ta j=0; j<D; j++) {
        real diff=x[i][j]-test[j];
	d2[i]+=diff*diff;
      }
    }

    wt_calc->maxn(&maxn);
    ind=new long[maxn];
    w=new real[maxn];

    k=(*wt_calc)(d2, n, ind, w);

    result=0;
    tw=0;
    for (nel_ta i=0; i<k; i++) {
      result+=w[i]*y[i];
      tw+=w[i];
    }
    result/=tw;

    if (e2 != NULL) {
      *e2=0;
      for (nel_ta i=0; i<k; i++) {
        real diff=w[i]*(y[i]-result);
        *e2+=diff*diff;
      }
    }

    delete [] ind;
    delete [] w;

    return result;
  }

  template cls_ta balloon_classify<float, cls_ta>(
		  balloon_param<float, long> *,
		  nel_ta, float **, cls_ta, cls_ta *, float *, float *);
  template cls_ta balloon_classify<double, cls_ta>(
		  balloon_param<double, long> *,
		  nel_ta, double **, cls_ta, cls_ta *, double *, double *);

  template cls_ta balloon_classify_joint<float, cls_ta>(
		  balloon_param<float, long> *,
		  nel_ta, float **, cls_ta, cls_ta *, float *, float *);
  template cls_ta balloon_classify_joint<double, cls_ta>(
		  balloon_param<double, long> *,
		  nel_ta, double **, cls_ta, cls_ta *, double *, double *);

  template float balloon_density<float>(
		  balloon_param<float, long> *,
		  nel_ta, float **, float *);
  template double balloon_density<double>(
		  balloon_param<double, long> *,
		  nel_ta, double **, double *);

  template float balloon_interpol<float>(
		  balloon_param<float, long> *,
		  nel_ta, float **, float *, float *, float *);
  template double balloon_interpol<double>(
		  balloon_param<double, long> *,
		  nel_ta, double **, double *, double *, double *);

} //end namespace libagf

