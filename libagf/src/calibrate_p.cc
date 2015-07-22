#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>

#include <gsl/gsl_linalg.h>

#include "full_util.h"

#include "agf_lib.h"

using namespace std;
using namespace libagf;
using namespace libpetey;

//compares two sets of class and outputs statistics on them:
int main(int argc, char ** argv) {
  FILE *fs;

  char *clsfile;
  char *confile;
  size_t slen;

  cls_ta *class1;		//true classes
  cls_ta *class2;		//classes after transformation
  nel_ta n1, n2;
  dim_ta ncls;
  dim_ta *nlab;

  real_a **p;

  gsl_matrix *a;
  gsl_vector *b;
  gsl_vector *x;

  gsl_matrix *vt;
  gsl_vector *s;
  gsl_vector *work;

  real_a **trans;		//transformation matrix
  real_a *cnst;
  real_a **p2;			//transformed p

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

  if (argc < 1) {
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

  slen=strlen(argv[0]);
  clsfile=new char[slen+5];
  strcpy(clsfile, argv[0]);
  strcat(clsfile, ".cls");
  confile=new char[slen+5];
  strcpy(confile, argv[0]);
  strcat(confile, ".vec");

  //read in the classes:
  class1=read_clsfile(clsfile, n1);
  if (n2 < 0) {
    fprintf(stderr, "Error reading input file: %s\n", clsfile);
    return ALLOCATION_FAILURE;
  }
  if (class1 == NULL) {
    fprintf(stderr, "Unable to open file for reading: %s\n", clsfile);
    return UNABLE_TO_OPEN_FILE_FOR_READING;
  }

  //read in confidence ratings:
  p=read_vecfile(confile, n2, ncls);
  if (n2 < 0) {
    fprintf(stderr, "Error reading input file: %s\n", confile);
    return ALLOCATION_FAILURE;
  }
  if (p == NULL) {
    fprintf(stderr, "Unable to open file for reading: %s\n", confile);
    return UNABLE_TO_OPEN_FILE_FOR_READING;
  }
  if (n1 != n2) {
    fprintf(stderr, "Data elements in files, %s and %s, do not agree: %d vs. %d\n", 
		    clsfile, confile, n1, n2);
    return SAMPLE_COUNT_MISMATCH;
  }

/*  for (int i=0; i<n1; i++) {
    for (int j=0; j<ncls; j++) {
      if (p[i][j]!=0) {
        p[i][j]=log(p[i][j]);
      } else {
        p[i][j]=log(1e-5);
      }
    }
  }*/

  //count the number of each class label:
  nlab=new cls_ta[ncls];
  for (int i=0; i<n1; i++) nlab[i]=0;
  for (int i=0; i<n1; i++) nlab[class1[i]]++;

  //fill up the matrix and result vector:
  a=gsl_matrix_alloc((ncls+1)*n1+ncls, ncls*(ncls+1));
  //a=gsl_matrix_alloc(2*n1, ncls*ncls);
  b=gsl_vector_alloc((ncls+1)*n1+ncls);
  //b=gsl_vector_alloc(2*n1);
  for (int i=0; i<n1; i++) {
    for (int j=0; j<ncls; j++) {
      for (int k=0; k<ncls; k++) {
        gsl_matrix_set(a, i*(ncls+1)+j, j*ncls+k, n1*p[i][k]/nlab[j]/nlab[j]);
        gsl_matrix_set(a, i*(ncls+1)+ncls, j*ncls+k, 1*p[i][k]);
        //gsl_matrix_set(a, i*2+1, j*ncls+k, 1*p[i][k]);
      }
      gsl_matrix_set(a, i*(ncls+1)+j, ncls*ncls+j, 1.);
      gsl_matrix_set(a, i*(ncls+1)+ncls, ncls*ncls+j, 1);
      if (class1[i]==j) {
        gsl_vector_set(b, i*(ncls+1)+j, (1.*n1)/nlab[j]/nlab[j]);
        //for (int k=0; k<ncls; k++) gsl_matrix_set(a, i*2, j*ncls+k, 1*p[i][k]);
        //gsl_vector_set(b, i*2, 1);
      } else {
        gsl_vector_set(b, i*(ncls+1)+j, 0);
      }
      gsl_vector_set(b, i*(ncls+1)+ncls, 1);
      //gsl_vector_set(b, i*2+1, 1);
    }
  }

  for (int i=0; i<ncls; i++) {
    real_a aval=0;
    for (int j=0; j<n1; j++) aval+=p[j][i];
    for (int j=0; j<ncls; j++) {
      gsl_matrix_set(a, n1*(ncls+1)+j, j*ncls+i, aval/nlab[i]/nlab[i]);
    }
    gsl_vector_set(b, n1*(ncls+1)+i, 1./nlab[i]);
  }
    

  vt=gsl_matrix_alloc(ncls*(ncls+1), ncls*(ncls+1));
  s=gsl_vector_alloc(ncls*(ncls+1));
  work=gsl_vector_alloc(ncls*(ncls+1));
  x=gsl_vector_alloc(ncls*(ncls+1));

  gsl_linalg_SV_decomp(a, vt, s, work);
  gsl_linalg_SV_solve(a, vt, s, b, x);

  //solution vector is transformation matrix:
  trans=allocate_matrix<real_a, nel_ta>(ncls, ncls);
  cnst=new real_a[ncls];
  p2=allocate_matrix<real_a, nel_ta>(n1, ncls);
  class2=new cls_ta[n1];

  for (int i=0; i<ncls; i++) {
    for (int j=0; j<ncls; j++) {
      trans[i][j]=gsl_vector_get(x, i*ncls+j);
    }
    cnst[i]=gsl_vector_get(x, ncls*ncls+i);
  }

  matrix_mult_t(p, trans, p2, n1, ncls, ncls);

  print_matrix(stdout, trans, ncls, ncls);
  printf("\n");
  for (int j=0; j<ncls; j++) printf("%g ", cnst[j]);
  printf("\n");
  printf("\n");

  for (int i=0; i<n1; i++) {
    class2[i]=0;
    for (int j=0; j<ncls; j++) {
      //p2[i][j]=exp(p2[i][j]);
      p2[i][j]+=cnst[j];
      printf("%g ", p2[i][j]);
      if (p2[i][j]>p2[i][class2[i]]) class2[i]=j;
    }
    printf("%d\n", class2[i]);
  }

  class_eval(class1, class2, n1);

  delete [] class1;
  delete [] class2;
  delete [] p;
  delete [] p2;

  delete_matrix(trans);

  gsl_matrix_free(a);
  gsl_matrix_free(vt);
  gsl_vector_free(b);
  gsl_vector_free(s);
  gsl_vector_free(work);
  gsl_vector_free(x);

  delete [] clsfile;
  delete [] confile;

  return exit_code;

}

