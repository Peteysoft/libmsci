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
  nel_ta nsamp;
  dim_ta ncls;

  real_a **p;			//probabilities
  real_a *ps;			//sorted probabilities
  long *sind;			//sort indices

  long midind;			//index of probability closest to 0.5

  double *sum;			//cumulator for even point intervals
  double *sump;			//sum of probabilities
  double *nacc;			//number of correct guesses
  double *rank;
  double brier=0;		//Brier score

  //slope and regression:
  double r, m;
  //stuff I don't need:
  double cov, sumsqr;

  flag_a Cflag=0;
  flag_a Hflag=0;
  flag_a bflag=0;
  int exit_code=0;

  char c;

  while ((c = getopt(argc, argv, "HCb")) != -1) {
    switch (c) {
      case ('C'):
             Cflag=1;
	     break;
      case ('H'):
	     Hflag=1;
	     break;
      case ('b'):
	     bflag=1;
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
    printf("  -b     = winning probabilities only\n");
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

  if (bflag) {
    char *yfile;
    char *confile;
    real_a *con;		//"confidence ratings"
    yfile=new char[strlen(pfile)+5];
    sprintf(yfile, "%s.cls", pfile);
    confile=new char[strlen(pfile)+5];
    sprintf(confile, "%s.con", pfile);
    class2=read_clsfile<cls_ta>(yfile, n2);
    if (n2 < 0 || class2==NULL) {
      fprintf(stderr, "Error reading input file, %s\n", yfile);
      exit(FILE_READ_ERROR);
    }
    if (n1 != n2) {
      fprintf(stderr, "Data elements in files, %s and %s, do not agree: %d vs. %d\n", 
		    tfile, yfile, n1, n2);
      return SAMPLE_COUNT_MISMATCH;
    }
    con=read_datfile<real_a>(confile, n2);
    if (n2 < 0 || con==NULL) {
      fprintf(stderr, "Error reading input file, %s\n", con);
      exit(FILE_READ_ERROR);
    }
    if (n1 != n2) {
      fprintf(stderr, "Data elements in files, %s and %s, do not agree: %d vs. %d\n", 
		    tfile, confile, n1, n2);
      return SAMPLE_COUNT_MISMATCH;
    }
    //convert confidence ratings back to conditional probabilities:
    p=new real_a *[1];
    p[0]=new real_a[n1];
    ncls=1;
    for (nel_ta i=0; i<n1; i++) if (class1[i]>=ncls) ncls=class1[i]+1;
    for (nel_ta i=0; i<n1; i++) p[0][i]=(con[i]*(ncls-1)+1)/ncls;
    printf("ncls=%d\n", ncls);
    nsamp=n1;
    delete [] con;
    delete [] yfile;
    delete [] confile;
  } else {
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
    nsamp=n1*ncls;
  }

  sind=heapsort(p[0], nsamp);
  ps=map_vector(p[0], sind, nsamp);

  //for (int i=0; i<nsamp; i++) printf("%g\n", ps[i]);
  
  //to save memory:
  delete [] p[0];
  delete [] p;

  midind=bin_search(ps, nsamp, (float) 1./ncls);
  //printf("%g\n", ps[n1*ncls/2]);
  sum=new double[nsamp+1];
  sump=new double[nsamp+1];		//sum of probabilities
  nacc=new double[nsamp+1];		//number of correct guesses
  sum[midind]=0;
  sump[midind]=0;
  nacc[midind]=0;
  for (int i=midind-1; i>=0; i--) {
    double diff;
    long k;
    cls_ta cls;
    if (bflag) {
      k=sind[i];
      cls=class2[k];
    } else {
      k=sind[i]/ncls;
      cls=sind[i]%ncls;
    }
    sump[i]=sump[i+1]+ps[i]-1;
    if (cls!=class1[k]) {
      sum[i]=sum[i+1]-1/(1-ps[i]);
      nacc[i]=nacc[i+1]-1;
      diff=ps[i];
    } else {
      sum[i]=sum[i+1];
      nacc[i]=nacc[i+1];
      diff=ps[i]-1;
    }
    brier+=diff*diff;
  }

  for (int i=midind; i<nsamp; i++) {
    double diff;
    long k;
    cls_ta cls;
    if (bflag) {
      k=sind[i];
      cls=class2[k];
    } else {
      k=sind[i]/ncls;
      cls=sind[i]%ncls;
    }
    sump[i+1]=sump[i]+ps[i];
    if (cls==class1[k]) {
      sum[i+1]=sum[i]+1/ps[i];
      nacc[i+1]=nacc[i]+1;
      diff=ps[i]-1;
    } else {
      sum[i+1]=sum[i];
      nacc[i+1]=nacc[i];
      diff=ps[i];
    }
    brier+=diff*diff;
  }

  //save more memory:
  delete [] sind;

  rank=new double[nsamp];
  if (ofile!=NULL) {
    fs=fopen(ofile, "w");
  }
  for (int i=0; i<nsamp; i++) {
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
  r=gsl_stats_correlation(nacc, 1, sump, 1, nsamp+1);
  exit_code=gsl_fit_mul(nacc, 1, sump, 1, nsamp+1, &m, &cov, &sumsqr);


  printf("correlation = %15.8lg\n", r);
  printf("slope       = %15.8lg\n", m);
  printf("Brier score = %15.8lg\n", sqrt(brier/(nsamp-1)));

  exit(0);

  delete [] class1;
  //delete [] class2;
  delete [] ps;
  delete [] sum;

  return exit_code;

}

