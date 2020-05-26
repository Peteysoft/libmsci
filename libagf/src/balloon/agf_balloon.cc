#include <math.h>

#include "balloon_param.h"

namespace libagf {
  //classify a test point using a variable bandwidth kernel density
  //estimator with constant weight:
  template <typename real, typename cls_t, typename index_t>
  cls_t balloon_classify(
		  balloon_param *wt_calc,	//classifier parameters
		  index_t n,			//number of training samples
		  real **x,			//coordinates
		  cls_t ncls,			//number of classes
		  cls_t *cls,			//discrete ordinates [0-ncls]
		  real *test,			//test point
		  real *p) {			//returned cond. prob.

    index_t D=wt_calc->dim();		//number of dim.
    real d2[n];				//squared distances
    index_t maxn=n;			//maximum number of weights

    index_t *ind;			//indices of samples in calc.
    real *w;				//weights
    index_t k;				//number of samples in calc.
    real tw;				//total of weights

    //calculate distances:
    for (index_t i=0; i<n; i++) {
      d2[i]=0;
      for (index_t j=0; j<D; j++) {
        real diff=x[i][j]-test[j];
	d2[i]+=diff*diff;
      }
    }

    //allocated space for indices and weights:
    wt_calc->maxn(&maxn);
    ind=new index_t[maxn];
    w=new index_t[maxn];

    //use parameter object to calculate weights:
    k=wt_calc(d2, n, ind, w);

    //sum weights to find probabilities:
    for (index_t i=0; i<ncls; i++) p[i]=0;
    tw=0;
    for (index_t i=0; i<k; i++) {
      tw+=w[i];
      p[cls[ind[i]]]+=w[i];
    }
    //normalize probabilities:
    for (index_t i=0; i<ncls; i++) p[i]/=tw;

    //clean up:
    delete [] ind;
    delete [] w;

    //choose best class:
    return choose_class(pdf, ncls);
  }

  template <typename real, typename cls_t, typename index_t>
  cls_t balloon_classify_joint(balloon_param *wt_calc,
		  index_t n,
		  real **x,
		  cls_t *cls,
		  cls_t ncls,
		  real *test,
		  real *p) {		//returned joint probabilities

    index_t D=wt_calc->dim();		//number of dimensions
    real d2[n];				//squared distances
    index_t maxn=n;			//maximum number of weights
    real norm;				//normalization coef.

    index_t *ind;			//index of points corr. to wts.
    real *w;				//weights
    index_t k;				//number of weights

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

    k=wt_calc(d2, n, ind, w, &norm);

    for (index_t i=0; i<ncls; i++) p[i]=0;
    tw=0;
    for (index_t i=0; i<k; i++) {
      p[cls[ind[i]]]+=w[i];
    }

    for (index_t i=0; i<ncls; i++) p[i]/=norm/n;

    delete [] ind;
    delete [] w;

    return choose_class(pdf, ncls);
  }

  //variable bandwidth kernel density estimator with constant weight:
  template <typename real, typename index_t>
  real balloon_density(balloon_param *wt_calc,
		  index_t n,
		  real **x,
		  real *test) {
    index_t D=wt_calc->dim();
    real d2[n];
    real norm;
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

    k=wt_calc(d2, n, ind, w, &norm);

    tw=0;
    for (index_t i=0; i<k; i++) tw+=w[i];

    delete [] ind;
    delete [] w;

    return tw/norm/n;
  }

  template <typename real, typename index_t>
  real balloon_interpol(balloon_param *wt_calc,
		  index_t n,
		  real **x,
		  real *y,
		  real *test,
		  real *e2) {
    index_t D=wt_calc->dim();
    real d2[n];
    index_t maxn=n;

    index_t *ind;
    real *w;
    index_t k;
    real tw;
    real result;

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

    result=0;
    tw=0;
    for (index_t i=0; i<k; i++) {
      result+=w[i]*y[i];
      tw+=w[i];
    }
    result/=tw;

    if (e2 != NULL) {
      *e2=0;
      for (index_t i=0; i<k; i++) {
        real diff=w[i]*(y[i]-result);
        *e2+=diff*diff;
      }
    }

    delete [] ind;
    delete [] w;

    return result;
  }

} //end namespace libagf

