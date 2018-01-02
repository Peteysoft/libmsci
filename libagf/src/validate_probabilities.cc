#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>

#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_fit.h>

#include "peteys_tmpl_lib.h"
#include "agf_lib.h"

using namespace std;
using namespace libagf;
using namespace libpetey;

//compares two sets of class and outputs statistics on them:
int main(int argc, char ** argv) {
  FILE *fs;
  char *tfile;
  char *pfile;
  char *ofile;
  size_t slen;

  cls_ta *class1;
  cls_ta *class2;
  nel_ta n1, n2;
  dim_ta ncls;

  real_a **p;			//probabilities
  real_a *ps;			//sorted probabilities
  long *sind;			//sort indices

  long midind;			//index of probability closest to 0.5

  double *sum;			//cumulator for even point intervals
  double *sump;			//sum of probabilities
  double *nacc;			//number of correct guesses
  double *rank;

  //slope and regression:
  double r, m;
  //stuff I don't need:
  double cov, sumsqr;

  flag_a Cflag=0;
  flag_a Hflag=0;
  int exit_code=0;

  char c;

  while ((c = getopt(argc, argv, "HC")) != -1) {
    switch (c) {
      case ('C'):
             Cflag=1;
	     break;
      case ('H'):
	     Hflag=1;
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
    printf("\nPoint-by-point validation of class conditional probabilities\n");
    printf("\n");
    printf("usage:  validate_probabilities class prob [output]\n\n");
    printf("where:\n");
    printf("  class  = binary file containing ture classes\n");
    printf("  prob   = ASCII file containing estimated probabilities\n");
    printf("  output = output file for plotting\n");
    printf("  -C     = no class data in ASCII file\n");
    printf("  -H     = no header in ASCII file\n");
    printf("\n");
    return INSUFFICIENT_COMMAND_ARGS;
  }

  tfile=argv[0];
  pfile=argv[1];
  if (argc > 2) ofile=argv[2]; else ofile=NULL;

  class1=read_clsfile<cls_ta>(tfile, n1);
  if (n1 < 0) {
    fprintf(stderr, "Error reading input file: %s\n", tfile);
    return ALLOCATION_FAILURE;
  }
  if (class1 == NULL) {
    fprintf(stderr, "Unable to open file for reading: %s\n", tfile);
    return UNABLE_TO_OPEN_FILE_FOR_READING;
  }
 
  n2=read_lvq(pfile, p, class2, ncls, Hflag+2*Cflag); 
  if (n2 < 0) {
    fprintf(stderr, "Error reading input file: %s\n", pfile);
    return ALLOCATION_FAILURE;
  }
  if (p == NULL) {
    fprintf(stderr, "Unable to open file for reading: %s\n", pfile);
    return UNABLE_TO_OPEN_FILE_FOR_READING;
  }
  if (n1 != n2) {
    fprintf(stderr, "Data elements in files, %s and %s, do not agree: %d vs. %d\n", 
		    tfile, pfile, n1, n2);
    return SAMPLE_COUNT_MISMATCH;
  }

  sind=heapsort(p[0], n1*ncls);
  ps=map_vector(p[0], sind, n1*ncls);

  //for (int i=0; i<n1*ncls; i++) printf("%g\n", ps[i]);
  
  //to save memory:
  delete [] p[0];
  delete [] p;

  midind=bin_search(ps, n1*ncls, (float) 1./ncls);
  //printf("%g\n", ps[n1*ncls/2]);
  sum=new double[n1*ncls+1];
  sump=new double[n1*ncls+1];		//sum of probabilities
  nacc=new double[n1*ncls+1];		//number of correct guesses
  sum[midind]=0;
  sump[midind]=0;
  nacc[midind]=0;
  for (int i=midind-1; i>=0; i--) {
    long k=sind[i]/ncls;
    cls_ta cls=sind[i]%ncls;
    sump[i]=sump[i+1]+ps[i]-1;
    if (cls!=class1[k]) {
      sum[i]=sum[i+1]-1/(1-ps[i]);
      nacc[i]=nacc[i+1]-1;
    } else {
      sum[i]=sum[i+1];
      nacc[i]=nacc[i+1];
    }
  }

  for (int i=midind+1; i<=n1*ncls; i++) {
    long k=sind[i-1]/ncls;
    cls_ta cls=sind[i]%ncls;
    sump[i]=sump[i-1]+ps[i-1];
    if (cls==class1[k]) {
      sum[i]=sum[i-1]+1/ps[i-1];
      nacc[i]=nacc[i-1]+1;
    } else {
      sum[i]=sum[i-1];
      nacc[i]=nacc[i-1];
    }
  }

  //save more memory:
  delete [] sind;

  rank=new double[n1*ncls];
  if (ofile!=NULL) {
    fs=fopen(ofile, "w");
  }
  for (int i=0; i<n1*ncls; i++) {
    rank[i]=i-midind;
    if (ofile!=NULL) {
      fprintf(fs, "%d %g %g %g %g\n", i-midind, ps[i], sum[i], nacc[i], sump[i]);
    }
  }
  if (ofile!=NULL) fclose(fs);

  //find correlation and slope:
  //lets waste some compute cycles by calculating them independently...
  //r=gsl_stats_correlation(rank, 1, sum, 1, n1*ncls);
  //exit_code=gsl_fit_mul(rank, 1, sum, 1, n1*ncls, &m, &cov, &sumsqr);
  r=gsl_stats_correlation(nacc, 1, sump, 1, n1*ncls);
  exit_code=gsl_fit_mul(nacc, 1, sump, 1, n1*ncls, &m, &cov, &sumsqr);


  printf("r   = %15.8lg\n", r);
  printf("m   = %15.8lg\n", m);
  printf("rms = %15.8lg\n", sqrt(sumsqr/(n1-1)));

  delete [] class1;
  delete [] class2;
  delete [] ps;
  delete [] sum;

  return exit_code;

}

