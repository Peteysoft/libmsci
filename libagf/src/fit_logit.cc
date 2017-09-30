#include <stdio.h>
#include <math.h>
#include <gsl/gsl_linalg.h>

#include "agf_lib.h"

#define ARBITRARY_LARGE 1000

void cost_function(*k, *skill, void *param) {
  void *p2[4]=(void **) param;
  int n=*(int *) p2[0];
  int D=*(int *) p2[1];
  real_a **x=(real_a**) p2[2];
  double *c=(double *) p2[3];
  cls_ta *truth=(cls_ta *) p2[4];
  real_a *result=new cls_ta[n];

  for (int i=0; i<n; i++) {
    real_a sum;
    sum=c[D];
    for (int j=0; j<D; j++) {
      sum+=c[j]*x[i][j];
    }
    result[i]=tanh(*k*sum);
  }




int main(int argc, char **argv) {
  FILE *fs=stdin;
  int n;
  real_a **train;
  cls_ta *cls;
  int D;

  n=read_lvq(fs, train, cls, D, 1);

  gsl_matrix *a=gsl_matrix_allocate(n, D+1);
  gsl_vector *b=gsl_vector_allocate(n);

  for (int i=0; i<n; i++) {
    for (int j=0; j<D; j++) {
      gsl_matrix_set(a, i, j, train[i][j]);
    }
    gsl_matrix_set(a, i, D, 1);
    gsl_vector_set(b, i, ARBITRARY_LARGE*(2*cls[i]-1));
  }

  gsl_vector *x=gsl_vector_allocate(D+1);

  gsl_solver(a, b, x);


