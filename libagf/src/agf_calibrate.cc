#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_multifit.h>

#include "peteys_tmpl_lib.h"
#include "roots_mins.h"

#include "agf_lib.h"

using namespace std;
using namespace libagf;
using namespace libpetey;

namespace libagf {
  template <typename real>
  real uc_bin(real ntn,			//number of true negatives
		real nfp,		// " false positives
		real nfn,		// " false negatives
		real ntp) {		// " true positives
    double uc, ucr, ucs, uc1;
    real **ctab=new real *[3];
    ctab[0]=new real[9];
    ctab[1]=ctab[0]+3;
    ctab[2]=ctab[0]+6;
  
    ctab[0][0]=ntn;
    ctab[0][1]=nfp;
    ctab[1][0]=nfn;
    ctab[1][1]=ntp;
    ctab[0][2]=ntn+nfp;
    ctab[1][2]=nfn+ntp;
    ctab[2][0]=ntn+nfn;
    ctab[2][1]=nfp+ntp;
    ctab[2][2]=ntn+nfp+nfn+ntp;
    //print_contingency_table(ctab, 2, 2);
    //printf("\n");
    uc=uncertainty_coefficient(ctab, 2, 2, ucr, ucs);

    delete [] ctab[0];
    delete [] ctab;

    return uc;
  }

  template <typename real>
  real acc_bin(real ntn,		//number of true negatives
		real nfp,		// " false positives
		real nfn,		// " false negatives
		real ntp) {		// " true positives
    return (ntp+ntn)/(ntn+nfp+nfn+ntp);
  }

  template <typename real>
  real corr_bin(real ntn,		//number of true negatives
		real nfp,		// " false positives
		real nfn,		// " false negatives
		real ntp) {		// " true positives
    return (ntn*ntp-nfn*nfp)/sqrt((ntp+nfp)*(ntn+nfn)*(ntp+nfn)*(nfp+ntn));
  }

  //skill as a function of the discrimination border:
  template <typename real>
  void skill_vs_r0(real r0, void *param, real *uc) {
    nel_ta ntp, ntn, nfp, nfn;		//confusion matrix
    double ind;				//location of r0
    void **p2=(void **) param;
    nel_ta n=*(nel_ta *) p2[0];		//number of classes
    real *r=(real *) p2[1];		//sorted retrieved conditional prob. (not sorted)
    nel_ta *none=(nel_ta *) p2[2];	//cumulative number of ones
    real (*sfunc) (real, real, real, real)=(real (*) (real, real, real, real)) p2[3];

    printf("interpolating r stuff\n");
    ind=interpolate(r, (long) n, r0);
    printf("index=%lg\n", ind);

    nfn=none[(int) ind];
    ntn=ind-nfn;
    ntp=none[n-1]-nfn;
    nfp=n-ind-ntp;

    *uc=-(*sfunc) (ntn, nfp, nfn, ntp);
  }

  //optimizes the skill of a binary classifier using a golden-section search
  template <typename real, typename cls_t>
  real optimize_binary_skill(cls_t *nonecum, real *rsort, nel_ta n, real (*binary_skill) (real, real, real, real), real tol, long maxiter) {
    long niter;
    real max;			//optimized value
    real r0;			//border at optimal skill
    void **param[4];

    param[0]=&n;
    param[1]=rsort;
    param[2]=nonecum;
    param[3]=(void *) binary_skill;

    printf("maximizing skill\n");
    r0=min_golden(&skill_vs_r0<real>, (void *) param, (real) -1.0, (real) 0.0, 
		(real) 1.0, tol, maxiter, niter, max);

    return r0;
  }

  //optimizes the skill of a binary classifier rigourously, by going through the probabilities 
  //one-by-one using each as the threshold value and calculating the skill score
  template <typename real>
  real optimize_binary_skill_rig(nel_ta *nonecum, real *rsort, nel_ta n, real (*binary_skill) (real, real, real, real)) {
    real nfn, ntn, ntp, nfp;
    real skill;		//calculated skill
    real max=-2;	//max. skill
    real r0;		//optimal threshold

    for (nel_ta i=0; i<n; i++) {
      nfn=nonecum[i];
      ntn=i-nfn;
      ntp=nonecum[n-1]-nfn;
      nfp=n-i-ntp;
      skill=(* binary_skill)(ntn, nfp, nfn, ntp);
      if (skill>max) {
        max=skill;
	r0=rsort[i];
      }
    }

    return r0;
  }

  //sorts the conditional probabilities, cumulates the class labels:
  //optimizes the skill of a binary classifier rigourously, by going through the probabilities 
  //one-by-one using each as the threshold value and calculating the skill score
  template <typename real>
  real * calc_skill_vs_r0(nel_ta *nonecum, real *rsort, nel_ta n, real (*binary_skill) (real, real, real, real)) {
    real nfn, ntn, ntp, nfp;
    real *skill=new real[n];	//calculated skill

    for (nel_ta i=0; i<n; i++) {
      nfn=nonecum[i];
      ntn=i-nfn;
      ntp=nonecum[n-1]-nfn;
      nfp=n-i-ntp;
      skill[i]=(* binary_skill)(ntn, nfp, nfn, ntp);
    }

    return skill;
  }

  //sorts the conditional probabilities, cumulates the class labels:
  template <typename real, typename cls_t>
  void sortr_cumulate_ones(cls_t *truth, real *r, nel_ta n, nel_ta *nonecum, real *rsort) {
    long *sind=heapsort(r, n);
    nonecum[0]=truth[sind[0]];
    for (nel_ta i=1; i<n; i++) {
      nonecum[i]=nonecum[i-1]+truth[sind[i]];
      rsort[i]=r[sind[i]];
    }
    delete [] sind;
  }

  //integrate the ROC curve:
  template <typename real>
  real integrate_roc(nel_ta *nonecum, real *rsort, nel_ta n) {
    real hitrate[n];
    real farate[n];
    real rocarea=0;

    for (nel_ta i=0; i<n; i++) {
      hitrate[i]=(real) (nonecum[n-1]-nonecum[i])/(real) nonecum[n-1];
      farate[i]=(real) (n-i+nonecum[i]-nonecum[n-1])/(real) (n-nonecum[n-1]);
    }

    for (nel_ta i=1; i<n-1; i++) {
      rocarea+=hitrate[i]*(farate[i]-farate[i-1]);
    }
    rocarea+=hitrate[0]*(farate[1]-farate[0])/2;
    rocarea+=hitrate[n-1]*(farate[n-1]-farate[n-2])/2;

    return -rocarea;

    //integrate ROC curve:
    //(stupid way of doing it dammit...
    real thetaold=M_PI/2;
    real Rold=sqrt(hitrate[0]*hitrate[0]+(1-farate[0])*(1-farate[0]));
    for (nel_ta i=1; i<n; i++) {
      real theta=atan(hitrate[i]/(1-farate[i]));
      real Rcur=sqrt(hitrate[i]*hitrate[i]+(1-farate[i])*(1-farate[i]));
      rocarea+=(thetaold-theta)*(Rold+Rcur)*(Rold+Rcur)/4;
      thetaold=theta;
      Rold=Rcur;
    }
    return rocarea;
  }

  template <typename real>
  real ** calc_roc(nel_ta *nonecum, real *rsort, nel_ta n) {
    real *hitrate;
    real *farate;
    real **result=new real *[2];

    result[0]=new real[2*n];
    result[1]=result[0]+n;
    farate=result[0];
    hitrate=result[1];

    for (nel_ta i=0; i<n; i++) {
      farate[i]=(real) (n-i+nonecum[i]-nonecum[n-1])/(real) (n-nonecum[n-1]);
      hitrate[i]=(real) (nonecum[n-1]-nonecum[i])/(real) nonecum[n-1];
    }

    return result;
  }

  template <typename real, typename cls_t>
  real * calibrate_r(cls_t *cls, real *r, nel_ta n, int order, nel_ta k) {
    real d[k];		//distances
    real acc[n-k];	//actual accuracies as a function of r
    real rave[n-k];	//sorted r which we transform with inverse tangent
    nel_ta none=0;	//number of ones in the zone of interest
    long *sind=heapsort(r, n);
    //for fitting:
    gsl_matrix *A;		//A^T A x = A^T b		
    gsl_vector *x;
    gsl_vector *b;
    //for the fitting:
    gsl_matrix *cov;
    gsl_multifit_linear_workspace *work;
    double chisq;

    //constant term is not supplied:
    int cflag=0;
    real cterm=0;

    //results
    real * coef;
 
    rave[0]=0;
    for (nel_ta i=0; i<k; i++) none+=cls[sind[i]];

    //interesting that this itself is a density estimation problem:
    for (nel_ta i=0; i<n-k; i++) {
      rave[i]=0;
      for (nel_ta j=0; j<k; j++) rave+=atanh(r[sind[i+j-k/2]]);
      rave[i]=rave[i]/k;
      acc[i]=atanh(2*none/k-1);
      none-=cls[sind[i]];
      none+=cls[sind[i+k]];
    }

    //allocate vectors and matrices that define the fitting problem:
    //best fit for A x = b
    //where a_ij = atanh^j(r_i), tanh(b_i) is actual accuracy and x are the fitting coefficients
    printf("border_calibrated: preparing for fitting\n");
    A=gsl_matrix_alloc(n, order+1-cflag);
    x=gsl_vector_alloc(order+1-cflag);
    b=gsl_vector_alloc(n);

    //fill vectors and matrices:
    for (int i=0; i<n; i++) {
      //printf("%g %g\n", tab[0][i], tab[1][i]);
      for (int j=cflag; j<=order; j++) {
        //here we don't care so much about efficiency because we only do the
	//fitting once and it's a small problem...
        gsl_matrix_set(A, i, j-cflag, pow(rave[i]+cterm, j));
      }
      gsl_vector_set(b, i, acc[i]);
    }

    //allocate extra stuff we need to perform the fitting:
    work=gsl_multifit_linear_alloc(n, order+1-cflag);
    cov=gsl_matrix_alloc(order+1-cflag, order+1-cflag);

    //perform fitting:
    printf("border_calibrated: performing fitting\n");
    gsl_multifit_linear(A, b, x, cov, &chisq, work);

    //pull coefficients out from GSL vector type:
    printf("border_calibrated: storing coefficients\n");
    coef=new real[order+1];
    for (int i=cflag; i<=order; i++) coef[i]=gsl_vector_get(x, i-cflag);
    if (cflag) coef[0]=cterm;

    //clean up:
    printf("calibrate_r: cleaning up\n");
    gsl_multifit_linear_free(work);
    gsl_matrix_free(A);
    gsl_matrix_free(cov);
    gsl_vector_free(x);
    gsl_vector_free(b);

    delete [] sind;

    //return the coefficients:
    return coef;
  }


  template void sortr_cumulate_ones<float, cls_ta>(cls_ta *truth, float *r, nel_ta n, nel_ta *nonecum, float *rsort);
  template void sortr_cumulate_ones<double, cls_ta>(cls_ta *truth, double *r, nel_ta n, nel_ta *nonecum, double *rsort);

  template float optimize_binary_skill_rig<float>(nel_ta *nonecum, float *rsort, nel_ta n, float (*binary_skill) (float, float, float, float));
  template double optimize_binary_skill_rig<double>(nel_ta *nonecum, double *rsort, nel_ta n, double (*binary_skill) (double, double, double, double));

  template float integrate_roc<float>(nel_ta *, float *, nel_ta);
  template double integrate_roc<double>(nel_ta *, double *, nel_ta);

  template float * calc_skill_vs_r0(nel_ta *, float *, nel_ta, float (*) (float, float, float, float));
  template double * calc_skill_vs_r0(nel_ta *, double *, nel_ta, double (*) (double, double, double, double));

  template float **calc_roc<float>(nel_ta *, float *, nel_ta);
  template double **calc_roc<double>(nel_ta *, double *, nel_ta);

  template float uc_bin<float>(float, float, float, float);
  template double uc_bin<double>(double, double, double, double);

  template float acc_bin<float>(float, float, float, float);
  template double acc_bin<double>(double, double, double, double);

  template float corr_bin<float>(float, float, float, float);
  template double corr_bin<double>(double, double, double, double);

} //end namespace libagf

