#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_statistics_double.h>

#include "full_util.h"
#include "peteys_tmpl_lib.h"

#include "agf_lib.h"

using namespace std;
using namespace libagf;
using namespace libpetey;

//compares two sets of class and outputs statistics on them:
int main(int argc, char ** argv) {
  FILE *fs;

  nel_ta n1, n2;
  dim_ta nc1;
  dim_ta nc2;

  real_a **p1;
  real_a **p2;

  real_a *ave1, *ave2;
  real_a *std1, *std2;

  real_a **corr;		//cross-correlation matrix

  gsl_matrix *a;
  gsl_vector *b;
  gsl_vector *x;

  gsl_matrix *vt;
  gsl_vector *s;
  gsl_vector *work;

  real_a **trans;		//transformation matrix
  real_a *cnst;

  int exit_code=0;
  char c;

  while ((c = getopt(argc, argv, "r:")) != -1) {
    switch (c) {
      case ('r'):
             fprintf(stderr, "optimal_r0: -r not implemented yet!\n");
	     break;
      case ('?'):
             fprintf(stderr, "Unknown option: -%c -- ignored\n", optopt);
	     exit_code=COMMAND_OPTION_PARSE_ERROR;
	     break;
      default:
	     fprintf(stderr, "Error parsing command line\n");
	     return FATAL_COMMAND_OPTION_PARSE_ERROR;
	     break;
    }
  }

  argc-=optind;
  argv+=optind;

  if (argc < 2) {
    printf("\nCompares two sets of classes.  Prints out a set of statistics\n");
    printf("Optionally checks the accuracy of estimated conditionaly probabilities\n");
    printf("by comparing with calculated accuracies\n\n");
    printf("usage:  cls_comp_stats [-b] file1 file2\n\n");
    printf("where:\n");
    printf("  file1  = binary file containing first set of classes ('truth')\n");
    printf("  file2  = binary file(s) containing second set of classes ('retrieved')\n");
    printf("             .cls for classification results\n");
    printf("             .con for confidence ratings\n\n");
    printf("  -b     = short output (uncertainty coefficient and accuracy)\n");
    printf("  -q     = number of histogram bins for evaluating confidence ratings [%d]\n", NCONHIST);
    printf("\n");
    return INSUFFICIENT_COMMAND_ARGS;
  }

  //read in the classes:
  p1=read_vecfile(argv[0], n1, nc1);
  if (n1 < 0) {
    fprintf(stderr, "Error reading input file: %s\n", argv[0]);
    return ALLOCATION_FAILURE;
  }
  if (p1 == NULL) {
    fprintf(stderr, "Unable to open file for reading: %s\n", argv[0]);
    return UNABLE_TO_OPEN_FILE_FOR_READING;
  }

  //read in confidence ratings:
  p2=read_vecfile(argv[1], n2, nc2);
  if (n2 < 0) {
    fprintf(stderr, "Error reading input file: %s\n", argv[1]);
    return ALLOCATION_FAILURE;
  }
  if (p2 == NULL) {
    fprintf(stderr, "Unable to open file for reading: %s\n", argv[1]);
    return UNABLE_TO_OPEN_FILE_FOR_READING;
  }
  if (n1 != n2) {
    fprintf(stderr, "Data elements in files, %s and %s, do not agree: %d vs. %d\n", 
		    argv[0], argv[1], n1, n2);
    return SAMPLE_COUNT_MISMATCH;
  }

  ave1=new real_a[nc1];
  ave2=new real_a[nc1];
  std1=new real_a[nc2];
  std2=new real_a[nc2];

  calc_norm(p1, nc1, n1, ave1, std1);
  calc_norm(p2, nc2, n1, ave2, std2);

  //norm_vec(p1, nc1, n1, ave1, std1);
  //norm_vec(p2, nc2, n1, ave2, std2);

  corr=zero_matrix<real_a, int32_t>(nc1, nc2);

  for (nel_ta i=0; i<n1; i++) {
    for (dim_ta j=0; j<nc1; j++) p1[i][j]-=ave1[j];
    for (dim_ta j=0; j<nc2; j++) p2[i][j]-=ave2[j];
  }

  for (dim_ta i=0; i<nc1; i++) {
    for (dim_ta j=0; j<nc2; j++) {
      for (nel_ta k=0; k<n1; k++) {
        corr[i][j]+=p1[k][i]*p2[k][j];
	//printf("%d: %g, %g\n", k, p1[k][i], p2[k][j]);
      }
      corr[i][j]/=std1[i]*std2[j]*(n1-1);
      printf("%10.6f ", corr[i][j]);
    }
    printf("\n");
  }

  gsl_matrix *u1, *u2;
  gsl_matrix *vt1, *vt2;
  gsl_vector *s1, *s2;
  gsl_vector *work1, *work2;
  real_a **corr2=zero_matrix<real_a, int32_t>(nc1, nc2);

  u1=gsl_matrix_alloc(n1, nc1);
  u2=gsl_matrix_alloc(n2, nc2);
  for (nel_ta i=0; i<n1; i++) {
    for (dim_ta j=0; j<nc1; j++) gsl_matrix_set(u1, i, j, p1[i][j]);
    for (dim_ta j=0; j<nc2; j++) gsl_matrix_set(u2, i, j, p2[i][j]);
  }

  vt1=gsl_matrix_alloc(nc1, nc1);
  vt2=gsl_matrix_alloc(nc2, nc2);
  s1=gsl_vector_alloc(nc1);
  s2=gsl_vector_alloc(nc2);
  work1=gsl_vector_alloc(nc1);
  work2=gsl_vector_alloc(nc2);

  gsl_linalg_SV_decomp(u1, vt1, s1, work1);
  gsl_linalg_SV_decomp(u2, vt2, s2, work2);

  printf("\n");
  for (dim_ta i=0; i<nc1; i++) {
    for (dim_ta j=0; j<nc2; j++) {
      for (nel_ta k=0; k<n1; k++) {
        corr2[i][j]+=gsl_matrix_get(u1, k, i)*gsl_matrix_get(u2, k, j);
	//printf("%d: %g, %g\n", k, p1[k][i], p2[k][j]);
      }
      //corr[i][j]/=(n1-1);
      printf("%10.6f ", corr2[i][j]);
    }
    printf("\n");
  }
  printf("\n");

  //do a simple partial correlation for off-diagonal elements only
  real_a m1, m2;
  real_a cov1, cov2, var1, var2;
  real_a **part=zero_matrix<real_a, int32_t>(nc1, nc2);

  for (dim_ta i=0; i<nc1; i++) {
    for (dim_ta j=0; j<nc2; j++) {
      if (i!=j) {
        cov1=0;
        cov2=0;
	var1=0;
	var2=0;
	for (nel_ta k=0; k<n1; k++) {
          cov1+=p1[k][i]*p1[k][j];
	  cov2+=p2[k][i]*p1[k][j];
	  var1+=p1[k][j]*p1[k][j];
	  //var2+=p2[k][i]*p2[k][i];
	}
	m1=cov1/var1;
	m2=cov2/var1;
	for (nel_ta k=0; k<n1; k++) {
          real_a diff1=p1[k][i]-p1[k][j]*m1;
          real_a diff2=p2[k][j]-p1[k][j]*m2;
	  part[i][j]+=diff1*diff2;
	}
	part[i][j]/=std1[i]*std2[j]*(n1-1);
      } else {
        part[i][j]=corr[i][j];
      }
      printf("%10.6f ", part[i][j]);
    }
    printf("\n");
  }

  long *sind[n1];			//sort each of the results and correlate
  				//by rank
  real_a corr3[nc2];

  for (nel_ta i=0; i<n1; i++) sind[i]=heapsort(p2[i], nc2);

  printf("\n");
  for (dim_ta i=0; i<nc2; i++) {
    double rms=0;
    double bias=0;
    //just pass this shit to the gsl routines:
    double vec1[n1];
    double vec2[n2];

    for (nel_ta j=0; j<n1; j++) {
      double diff;
      vec1[j]=p1[j][sind[j][i]];
      vec2[j]=p2[j][sind[j][i]];
      diff=vec2[j]-vec1[j];
      bias+=diff;
      rms+=diff*diff;
    }

    corr3[i]=gsl_stats_correlation(vec1, 1, vec2, 1, n1);
    //printf("%10.6g ", sqrt(rms/(n1-1)));
    printf("%10.6g ", bias/n1);
  }
  printf("\n");
  for (dim_ta i=0; i<nc2; i++) {
    printf("%10.6g ", corr3[i]);
  }
  
  //print_matrix(stdout, corr, nc1, nc2);

  exit(0);

/*  for (int i=0; i<n1; i++) {
    for (int j=0; j<ncls; j++) {
      if (p[i][j]!=0) {
        p[i][j]=log(p[i][j]);
      } else {
        p[i][j]=log(1e-5);
      }
    }
  }*/

  //fill up the matrix and result vector:
  a=gsl_matrix_alloc(nc2*n1, (nc1+1)*nc2);
  //a=gsl_matrix_alloc(2*n1, ncls*ncls);
  b=gsl_vector_alloc(nc2*n1);
  //b=gsl_vector_alloc(2*n1);
  for (int i=0; i<n1; i++) {
    for (int j=0; j<nc2; j++) {
      for (int k=0; k<nc1; k++) {
        gsl_matrix_set(a, i*nc2+j, j*nc1+k, p1[i][k]);
      }
      gsl_matrix_set(a, i*nc2+j, nc2*nc1+j, 1);		//const. term
      gsl_vector_set(b, i*nc2+j, p2[i][j]);
    }
  }

  vt=gsl_matrix_alloc((nc1+1)*nc2, (nc1+1)*nc2);
  s=gsl_vector_alloc((nc1+1)*nc2);
  work=gsl_vector_alloc((nc1+1)*nc2);
  x=gsl_vector_alloc((nc1+1)*nc2);

  gsl_linalg_SV_decomp(a, vt, s, work);
  gsl_linalg_SV_solve(a, vt, s, b, x);

  //solution vector is transformation matrix:
  trans=allocate_matrix<real_a, nel_ta>(nc2, nc1);
  cnst=new real_a[nc2];

  for (int i=0; i<nc2; i++) {
    for (int j=0; j<nc1; j++) {
      trans[i][j]=gsl_vector_get(x, i*nc1+j);
    }
    cnst[i]=gsl_vector_get(x, nc1*nc2+i);
  }

  print_matrix(stdout, trans, nc2, nc1);
  printf("\n");
  for (int j=0; j<nc2; j++) printf("%g ", cnst[j]);
  printf("\n");
  printf("\n");

  delete [] p1;
  delete [] p2;

  delete_matrix(trans);

  gsl_matrix_free(a);
  gsl_matrix_free(vt);
  gsl_vector_free(b);
  gsl_vector_free(s);
  gsl_vector_free(work);
  gsl_vector_free(x);

  return exit_code;

}

