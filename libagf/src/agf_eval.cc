
//Go ahead and copy this piece of shit all you like.  And don't bother to give me a dime or 
//or even a thank you...

#include <math.h>
#include <stdint.h>

#include "agf_defs.h"
#include "agf_eval.h"

using namespace std;

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

template <class cls_t>
double uncertainty_coefficient(nel_ta **acc_mat, cls_t nclt, cls_t nclr, double &ucr, double &uct) {
  double hrt, hr, ht;		//entropy measures
  //we've simplified the computation somewhat:
  nel_ta nt=acc_mat[nclt][nclr];
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

} //end namespace libagf
