
#include <math.h>
#include <stdint.h>

#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_fit.h>

#include "peteys_tmpl_lib.h"

#include "agf_defs.h"
#include "agf_eval.h"

using namespace std;
using namespace libpetey;

namespace libagf {

//contingency table includes row and column subtotals
template <class cls_t>
nel_ta **build_contingency_table(cls_t *truth, cls_t *ret, nel_ta n, cls_t &nclt, cls_t &nclr) {
  nel_ta **acc_mat;		//histogram (joint probability)

  nclt=0;
  nclr=0;
  for (nel_ta i=0; i<n; i++) {
    if (truth[i]>=nclt) nclt=truth[i]+1;
    if (ret[i]>=nclr) nclr=ret[i]+1;
  }

  acc_mat=new nel_ta *[nclt+1];
  acc_mat[0]=new nel_ta[(nclt+1)*(nclr+1)];
  for (cls_t i=1; i<=nclt; i++) acc_mat[i]=&acc_mat[0][i*(nclr+1)];
  for (long i=0; i<(nclr+1)*(nclt+1); i++) acc_mat[0][i]=0;
  
  for (nel_ta i=0; i<n; i++) {
    acc_mat[truth[i]][ret[i]]++;
  }

  for (nel_ta i=0; i<nclt; i++) {
    for (nel_ta j=0; j<nclr; j++) {
      acc_mat[nclt][j]+=acc_mat[i][j];		//column totals
      acc_mat[i][nclr]+=acc_mat[i][j];		//row totals
    }
    acc_mat[nclt][nclr]+=acc_mat[i][nclr];	//number of elements
  }

  return acc_mat;
}

template <typename cls_t, typename count_t>
double uncertainty_coefficient(count_t **acc_mat, cls_t nclt, cls_t nclr, double &ucr, double &uct) {
  double hrt, hr, ht;		//entropy measures
  //we've simplified the computation somewhat:
  count_t nt=acc_mat[nclt][nclr];
  double ntlognt=nt*log(nt);
  double uc;
  //calculate the uncertainty coefficient:
  
  hrt=0;
  hr=0;
  ht=0;

  for (cls_t i=0; i<nclt; i++) {
    for (cls_t j=0; j<nclr; j++) {
      if (acc_mat[i][j] != 0) {
        hrt-=acc_mat[i][j]*log(acc_mat[i][j]);
      }
    }
    if (acc_mat[i][nclr] != 0) {
      ht-=acc_mat[i][nclr]*log(acc_mat[i][nclr]);
    }
  }
  for (cls_t i=0; i<nclr; i++) {
    if (acc_mat[nclt][i] != 0) {
      hr-=acc_mat[nclt][i]*log(acc_mat[nclt][i]);
    }
  }

  uc=ht-hrt+hr+ntlognt;
  ucr=uc/(hr+ntlognt);
  uct=2*uc/(hr+ht+2*ntlognt);

  return uc/(ht+ntlognt);
}
  
template <class cls_t>
void print_contingency_table(nel_ta **acc_mat, cls_t nclt, cls_t nclr, FILE *fs) {
  fprintf(fs, "        retrieval\n");
  fprintf(fs, "truth ");
  for (cls_t i=0; i<nclr; i++) fprintf(fs, "%5d ", i);
  fprintf(fs, "    tot\n");
  for (cls_t i=0; i<nclt; i++) {
    fprintf(fs, "%3d   ", i);
    for (cls_t j=0; j<nclr; j++) {
      fprintf(fs, "%5d ", acc_mat[i][j]);
    }
    fprintf(fs, " %6d\n", acc_mat[i][nclr]);
  }
  fprintf(fs, "\ntotal ");
  for (cls_t i=0; i<nclr; i++) fprintf(fs, "%5d ", acc_mat[nclt][i]);
  fprintf(fs, "\n");
}
      
//calculate the uncertainty coefficient:
template <class cls_t>
double class_eval(cls_t *truth, cls_t *ret, nel_ta n, FILE *fs) {
  nel_ta **acc_mat;		//histogram (joint probability)
  double uc, ucr, uct;		//uncertainty coefficients
  cls_t nclt, nclr;		//number of classes
  nel_ta nt;			//number of true classes

  acc_mat=build_contingency_table(truth, ret, n, nclt, nclr);
  uc=uncertainty_coefficient(acc_mat, nclt, nclr, ucr, uct);

  nt=0;
  for (cls_t i=0; i<nclt && i<nclr; i++) nt+=acc_mat[i][i];

  if (fs != NULL) {
    print_contingency_table(acc_mat, nclt, nclr, fs);
    fprintf(fs, "\n");
    fprintf(fs, "U C (reverse): %f\n", (float) ucr);
    fprintf(fs, "U C (total): %f\n", (float) uct);
    fprintf(fs, "Uncertainty coefficient: %f\n", (float) uc);
    printf("Accuracy: %f\n", (float) nt/n);
  }

  //clean up:
  delete [] acc_mat[0];
  delete [] acc_mat;
  
  return uc;
  
}

//calculate the uncertainty coefficient:
template <class cls_t>
double class_eval_basic(cls_t *truth, cls_t *ret, nel_ta n, FILE *fs, flag_a Hflag) {
  nel_ta **acc_mat;		//histogram (joint probability)
  double uc, ucr, uct;		//uncertainty coefficients
  cls_t nclt, nclr;		//number of classes
  nel_ta nt;			//number of true classes

  acc_mat=build_contingency_table(truth, ret, n, nclt, nclr);
  uc=uncertainty_coefficient(acc_mat, nclt, nclr, ucr, uct);

  nt=0;
  for (cls_t i=0; i<nclt && i<nclr; i++) nt+=acc_mat[i][i];

  if (fs != NULL) {
    if (Hflag!=1) fprintf(fs, "Acc.   U.C.\n");
    //make it easier to read off in automated scripts:
    fprintf(fs, "%6.4f %6.4f\n", (float) nt/n, uc);
  }

  //clean up:
  delete [] acc_mat[0];
  delete [] acc_mat;
  
  return uc;
  
}

//test the accuracy of the confidence ratings:
template <class real, class cls_t>
real ** con_acc_table(cls_t *truth, cls_t *cls, real *con, nel_ta n, int nhist) {
  nel_ta total[nhist];
  nel_ta ind;
  real acc;
  real **result;

  result=new real*[2];
  result[0]=new real[2*nhist];
  result[1]=result[0]+nhist;

  for (nel_ta i=0; i<nhist; i++) {
    total[i]=0;
    result[0][i]=0;
    result[1][i]=0;
  }

  for (nel_ta i=0; i<n; i++) {
    ind=(nel_ta) (con[i]*nhist);
    if (ind<0) ind=0; else if (ind >= nhist) ind=nhist-1;
    total[ind]++;
    if (cls[i] == truth[i]) result[1][ind]++;
    result[0][ind] += con[i];
  }

  for (nel_ta i=0; i<nhist; i++) {
    result[0][i]=result[0][i]/total[i];
    result[1][i]=(real) result[1][i]/(real) total[i];
  }

  return result;

}

//test the accuracy of the confidence ratings
//for binary classifiers
template <class real, class cls_t>
real ** con_acc_table2(cls_t *truth, real *r, nel_ta n, int nhist) {
  nel_ta total[nhist*2];
  nel_ta ind;
  real acc;
  real **result;

  result=new real*[2];
  result[0]=new real[4*nhist];
  result[1]=result[0]+2*nhist;

  for (nel_ta i=0; i<nhist*2; i++) {
    total[i]=0;
    result[0][i]=0;
    result[1][i]=0;
  }

  for (nel_ta i=0; i<n; i++) {
    cls_t cls=convertR(r[i]);
    ind=(int) ((r[i]+1)*nhist);
    if (ind<0) ind=0; else if (ind >= 2*nhist) ind=2*nhist-1;
    total[ind]++;
    if (cls == truth[i]) result[1][ind]++;
    result[0][ind] += r[i];
  }

  for (nel_ta i=0; i<2*nhist; i++) {
    result[0][i]=result[0][i]/total[i];
    result[1][i]=2*(real) result[1][i]/(real) total[i]-1;
    if (i<nhist) result[1][i]=-result[1][i];
  }

  return result;

}

//for data read from the output files generated by libagf classifiers:
template <class real, class cls_t>
real ** con_acc_table2(cls_t *truth, cls_t *cls, real *con, nel_ta n, int nhist) {
  real r[n];
  for (int i=0; i<n; i++) r[i]=(2*cls[i]-1)*con[i];
  return con_acc_table2(truth, r, n, nhist);
}

//print the table accuracy vs. confidence rating:
template <class real>
void print_con_acc(real **table, int ncls, int nhist, FILE *fs) {

  fprintf(fs, "\nAccuracy vs. confidence:\n\n");
  for (nel_ta i=0; i<nhist; i++) {
    fprintf(fs, "%10.3f %10.3f\n", table[0][i], (ncls*table[1][i]-1.)/(ncls-1.));
  }
  fprintf(fs, "\n");

}

template <class real, class cls_t>
void check_confidence(cls_t *truth, cls_t *cls, real *con, nel_ta n, int nhist, FILE *fs) {
  real **table;
  cls_t ncls=1;
  //hell with it, just count 'em over:
  for (nel_ta i=0; i<n; i++) if (truth[i]>=ncls) ncls=truth[i]+1;

  if (ncls==2) {
    table=con_acc_table2(truth, cls, con, n, nhist);
    fprintf(fs, "\nAccuracy vs. probability difference:\n\n");
    for (nel_ta i=0; i<nhist*2; i++) {
      fprintf(fs, "%10.3f %10.3f\n", table[0][i], table[1][i]);
    }
    fprintf(fs, "\n");
  } else {
    table=con_acc_table(truth, cls, con, n, nhist);
    print_con_acc(table, ncls, nhist, fs);
  }

  delete [] table[0];
  delete [] table;
}

template <typename real, typename cls_t>
struct prob_struct {
  flag_a t;			//true or false?
  real p;			//probability
};

template <typename real, typename cls_t>
int prob_struct_comp(const void * v1, const void * v2) {
   prob_struct<real, cls_t> *p1=(prob_struct<real, cls_t> *) v1;
   prob_struct<real, cls_t> *p2=(prob_struct<real, cls_t> *) v2;
   if (p1->p > p2->p) {
     return 1;
   } else if (p1->p < p2->p) {
     return -1;
   } else if (p1->t > p2->t) {
     return 1;
   } else if (p1->t < p1->t) {
     return -1;
   } else {
     return 0;
   }
}
    
template <typename real, typename cls_t>
int validate_cond_prob(prob_struct<real, cls_t> *data, nel_ta nsamp, real thresh,
		double &corr, double &slope, real &brier, double &sumsqr, double &norm, FILE *fs) {
  real *ps;			//sorted probalities
  double *sum=NULL;		//for earlier, equal-step version
  double *sump;			//sum of probabilities
  double *nacc;			//number of correct guesses
  int midind;
  double r, m;			//correlation, slope
  double cov;			//covariance
  int exit_code=0;

  qsort(data, nsamp, sizeof(prob_struct<real, cls_t>), &prob_struct_comp<real, cls_t>);

  ps=new real[nsamp];
  for (nel_ta i=0; i<nsamp; i++) {
    ps[i]=data[i].p;
    //printf("%g\n", ps[i]);
  }

  midind=bin_search(ps, nsamp, thresh);
  if (midind<0) midind=0; else if (midind >= nsamp) midind=nsamp-1;

  delete [] ps;

  if (fs!=NULL) sum=new double[nsamp+1];
  sump=new double[nsamp+1];		//sum of probabilities
  nacc=new double[nsamp+1];		//number of correct guesses
  if (fs!=NULL) sum[midind]=0;
  sump[midind]=0;
  nacc[midind]=0;
  brier=0;
  norm=0;
  for (int i=midind-1; i>=0; i--) {
    real diff;
    sump[i]=sump[i+1]+data[i].p-1;
    nacc[i]=nacc[i+1]+data[i].t-1;
    diff=data[i].p-data[i].t;
    if (fs!=NULL) sum[i]=sum[i+1]-(data[i].t-1)/(1-data[i].p);
    brier+=diff*diff;
    norm+=(midind-i)*(midind-i);
  }

  for (int i=midind; i<nsamp; i++) {
    real diff;
    sump[i+1]=sump[i]+data[i].p;
    nacc[i+1]=nacc[i]+data[i].t;
    diff=data[i].p-data[i].t;
    if (fs!=NULL) sum[i+1]=sum[i]+data[i].t/data[i].p;
    brier+=diff*diff;
    norm+=(i-midind)*(i-midind);
  }

  if (fs!=NULL) {
    for (int i=0; i<nsamp; i++) {
      //rank[i]=i-midind;
      fprintf(fs, "%d %g %lg %lg %lg\n", i-midind, data[i].p, sum[i], nacc[i], sump[i]);
    }
  }

  //calculate correlation and fit the slope:
  corr=gsl_stats_correlation(nacc, 1, sump, 1, nsamp+1);
  exit_code=gsl_fit_mul(nacc, 1, sump, 1, nsamp+1, &slope, &cov, &sumsqr);

  //printf("sumsqr = %lg\n", sqrt(sumsqr/(nsamp-1)));
  //printf("norm. sumsqr = %lg\n", sqrt(sumsqr/norm));

  //clean up:
  delete [] sum;
  delete [] sump;
  delete [] nacc;
  delete [] data;

  return exit_code;

}

template <typename real, typename cls_t>
int validate_cond_prob(cls_t *class1, real **p, nel_ta n, cls_t ncls,
	real &corr, real &slope, real &brier, 
	real &rms, real &norm, FILE *fs) {
  double r, m;
  double norm1, sumsqr;
  int exit_code=0;
  prob_struct<real, cls_t> *data=new prob_struct<real, cls_t>[n*ncls];

  for (nel_ta i=0; i<n; i++) {
    for (cls_t j=0; j<ncls; j++) {
      data[i*ncls+j].p=p[i][j];
      data[i*ncls+j].t=(class1[i]==j);
    }
  }

  exit_code=validate_cond_prob(data, n*ncls, (real) 1./ncls, r, m, brier, 
		  sumsqr, norm1, fs);

  corr=r;
  slope=m;
  rms=sqrt(sumsqr/n);
  norm=sqrt(norm1/n);

  brier=sqrt(brier/n/ncls);

  return exit_code;
}

template <typename real, typename cls_t>
int validate_cond_prob(cls_t *class1, real *p, cls_t *class2, nel_ta n,
	real &corr, real &slope, real &brier, 
	real &rms, real &norm, FILE *fs) {
  double r, m;
  cls_t ncls;
  double norm1, sumsqr;
  int exit_code=0;
  prob_struct<real, cls_t> *data=new prob_struct<real, cls_t>[n];

  ncls=1;
  for (int i=0; i<n; i++) {
    data[i].p=p[i];
    data[i].t=(class1[i]==class2[i]);
    if (class1[i]>=ncls) ncls=class1[i]+1;
  }

  exit_code=validate_cond_prob(data, n, (real) 1./ncls, r, m, brier, 
		  sumsqr, norm1, fs);

  corr=r;
  slope=m;
  rms=sqrt(sumsqr/n);
  norm=sqrt(norm1/n);

  brier=sqrt(brier/n);

  return exit_code;
}

template double uncertainty_coefficient<cls_ta, nel_ta>(nel_ta **acc_mat, cls_ta nclt, cls_ta nclr, double &ucr, double &uct);
template double uncertainty_coefficient<cls_ta, float>(float **acc_mat, cls_ta nclt, cls_ta nclr, double &ucr, double &uct);
template double uncertainty_coefficient<cls_ta, double>(double **acc_mat, cls_ta nclt, cls_ta nclr, double &ucr, double &uct);

template double class_eval<int32_t>(int32_t *, int32_t *, nel_ta, FILE *);
template double class_eval_basic<int32_t>(int32_t *, int32_t *, nel_ta, FILE *, flag_a);
template float ** con_acc_table<float, int32_t>(int32_t *, int32_t *, float *, nel_ta, int);
template double ** con_acc_table<double, int32_t>(int32_t *, int32_t *, double *, nel_ta, int);
template float ** con_acc_table2<float, int32_t>(int32_t *, int32_t *, float *, nel_ta, int);
template double ** con_acc_table2<double, int32_t>(int32_t *, int32_t *, double *, nel_ta, int);
template float ** con_acc_table2<float, int32_t>(int32_t *, float *, nel_ta, int);
template double ** con_acc_table2<double, int32_t>(int32_t *, double *, nel_ta, int);
template void print_con_acc<float>(float **, int, int, FILE *);
template void print_con_acc<double>(double **, int, int, FILE *);
template void check_confidence<float, int32_t>(int32_t *, int32_t *, float *,
		nel_ta n, int nhist, FILE *);
template void check_confidence<double, int32_t>(int32_t *, int32_t *, double *,
		nel_ta n, int nhist, FILE *);

template int validate_cond_prob<float, int32_t>(int32_t *, float **, nel_ta, int32_t, 
	float &, float &, float &, float &, float &, FILE *);
template int validate_cond_prob<double, int32_t>(int32_t *, double **, nel_ta, int32_t, 
	double &, double &, double &, double &, double &, FILE *);
template int validate_cond_prob<float, int32_t>(int32_t *, float *, int32_t *, nel_ta, 
	float &, float &, float &, float &, float &, FILE *);
template int validate_cond_prob<double, int32_t>(int32_t *, double *, int32_t *, nel_ta,
	double &, double &, double &, double &, double &, FILE *);

} //end namespace libagf
