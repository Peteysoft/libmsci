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
  FILE *fs=NULL;
  char *tfile;
  char *pfile;
  char *ofile;
  size_t slen;

  cls_ta *class1;
  cls_ta *class2=NULL;
  real_a **p;
  nel_ta n1, n2;
  dim_ta ncls;

  real_a brier;			//Brier score

  //slope and regression:
  real_a r, m;
  //other dumb measures that don't really work as well as the Brier score:
  real_a rms, norm;

  flag_a Cflag=0;
  flag_a Hflag=0;
  flag_a bflag=0;
  flag_a Bflag=0;
  int exit_code=0;

  char c;

  while ((c = getopt(argc, argv, "HCbB")) != -1) {
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
      case ('B'):
	     Bflag=1;
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
    printf("  -B     = Brier score only\n");
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
    //printf("ncls=%d\n", ncls);
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
  }

  if (ofile!=NULL) {
    fs=fopen(ofile, "w");
  }

  if (bflag) {
    validate_cond_prob(class1, p[0], class2, n1, r, m, brier, rms, norm, fs);
  } else {
    validate_cond_prob(class1, p, n1, ncls, r, m, brier, rms, norm, fs);
  }

  if (ofile!=NULL) fclose(fs);

  if (Bflag) {
    printf("%g\n", brier);
  } else {
    printf("1-r         = %15.8g\n", 1.-r);
    printf("m-1         = %15.8g\n", m-1.);
    printf("Brier score = %15.8g\n", brier);
    printf("norm. rmse  = %15.8g\n", rms/norm);
  }

  //clean up:
  delete [] p[0];
  delete [] p;
  delete [] class1;
  delete [] class2;

  return exit_code;

}

