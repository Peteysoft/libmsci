#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>

#include "peteys_tmpl_lib.h"
#include "roots_mins.h"

#include "agf_lib.h"

using namespace std;
using namespace libagf;
using namespace libpetey;

//uncertainty coefficient as a function of the discrimination border:
void uc_vs_r0(real_a r0, void *param, real_a *uc) {
  void **p2=(void **) param;
  nel_ta n=*(nel_ta *) p2[0];		//number of classes
  cls_ta *truth=(cls_ta *) p2[1];	//"true" classes
  real_a *r=(real_a *) p2[2];		//retrieved conditional prob.
  cls_ta *ret=(cls_ta *) p2[3];		//classes derived from cond. prob. & r0
  cls_ta **ctab;			//contingency table
  cls_ta nclt, nclr;			//numbers of class labels
  double uc1, ucr, ucs;			//uncertainty coefficient

  for (nel_ta i=0; i<n; i++) {
    if (r[i]>r0) ret[i]=1; else ret[i]=0;
  }
  
  ctab=build_contingency_table(truth, ret, n, nclt, nclr);
  uc1=uncertainty_coefficient(ctab, nclt, nclr, ucr, ucs);

  delete ctab[0];
  delete ctab;

  *uc=-uc1;

}
  

//compares two sets of class and outputs statistics on them:
int main(int argc, char ** argv) {
  FILE *fs;

  char *tfile;
  char *clsfile;
  char *confile;
  size_t slen;

  cls_ta *class1;
  cls_ta *class2;
  nel_ta n1, n2;

  real_a *con;

  real_a uc_max;		//maximum uncertainty coeff.
  real_a r0;			//optimal threshold

  void *param[4];		//parameters to pass to function

  long maxiter=1000;
  long niter;

  cls_ta nc1, nc2;		//numbers of each class

  int exit_code=0;

  char c;

  while ((c = getopt(argc, argv, "r:")) != -1) {
    switch (c) {
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
    printf("\nOptimizes the decision border.\n");
    printf("\n\n");
    printf("usage:  optimize_r0 [-b] file1 file2\n\n");
    printf("where:\n");
    printf("  file1  = binary file containing true class values\n");
    printf("  file2  = pair of binary files containing corresponding decision functions:\n");
    printf("             .cls for classification results\n");
    printf("             .con for confidence ratings\n\n");
    printf("  -b     = optimize based on uncertainty coefficient\n");
    printf("\n");
    return INSUFFICIENT_COMMAND_ARGS;
  }

  tfile=argv[0];
  slen=strlen(argv[1]);
  clsfile=new char[slen+5];
  strcpy(clsfile, argv[1]);
  strcat(clsfile, ".cls");
  confile=new char[slen+5];
  strcpy(confile, argv[1]);
  strcat(confile, ".con");

  //read in "true" classes:
  class1=read_clsfile(tfile, n1);
  if (n1 < 0) {
    fprintf(stderr, "Error reading input file: %s\n", tfile);
    return ALLOCATION_FAILURE;
  }
  if (class1 == NULL) {
    fprintf(stderr, "Unable to open file for reading: %s\n", tfile);
    return UNABLE_TO_OPEN_FILE_FOR_READING;
  }
  
  //read in retrieved classes:
  class2=read_clsfile(clsfile, n2);
  if (n2 < 0) {
    fprintf(stderr, "Error reading input file: %s\n", clsfile);
    return ALLOCATION_FAILURE;
  }
  if (class2 == NULL) {
    fprintf(stderr, "Unable to open file for reading: %s\n", clsfile);
    return UNABLE_TO_OPEN_FILE_FOR_READING;
  }
  if (n1 != n2) {
    fprintf(stderr, "Data elements in files, %s and %s, do not agree: %d vs. %d\n", 
		    tfile, clsfile, n1, n2);
    return SAMPLE_COUNT_MISMATCH;
  }

  //read in confidence ratings:
  con=read_datfile(confile, n2);
  if (n2 < 0) {
    fprintf(stderr, "Error reading input file: %s\n", confile);
    return ALLOCATION_FAILURE;
  }
  if (con == NULL) {
    fprintf(stderr, "Unable to open file for reading: %s\n", confile);
    return UNABLE_TO_OPEN_FILE_FOR_READING;
  }
  if (n1 != n2) {
    fprintf(stderr, "Data elements in files, %s and %s, do not agree: %d vs. %d\n", 
		    tfile, confile, n1, n2);
    return SAMPLE_COUNT_MISMATCH;
  }

  //convert class and confidence rating to difference in cond. prob.:
  nc1=0;
  nc2=0;
  for (int i=0; i<n1; i++) {
    if (class2[i]<1) {
      con[i]=-con[i];
      nc1++;
    } else {
      nc2++;
    }
  }

  //fill up parameters to pass to minimization routine:
  param[0]=&n1;
  param[1]=class1;
  param[2]=con;
  param[3]=class2;

  r0=min_golden(&uc_vs_r0, (void *) param, (real_a) -1.0, (real_a) 0.5, 
		(real_a) 1.0, (real_a) 1e-7, maxiter, niter, uc_max);

  uc_max=-uc_max;

  printf("r0 = %g; uc = %g\n", r0, uc_max);

  long *sind=heapsort(con, n1);
  fs=fopen("rhist.txt", "w");
  for (int i=0; i<n1; i++) {
    fprintf(fs, "%d %g\n", i, con[sind[i]]);
  }
  fclose(fs);

  double ucm2=0;
  nel_ta ind=0;

  fs=fopen("accvsr0.txt", "w");
  for (int i=0; i<n1; i++) {
    double ucr, ucs, uc1;
    cls_ta **ctab;
    cls_ta nclt, nclr;

    r0=con[sind[i]];
    for (int j=0; j<n1; j++) if (con[j]<r0) class2[j]=0; else class2[j]=1;
    ctab=build_contingency_table(class1, class2, n1, nclt, nclr);
    uc1=uncertainty_coefficient(ctab, nclt, nclr, ucr, ucs);


    //fprintf(fs, "%g", r0);
    //class_eval_basic(class1, class2, n1, fs);
    fprintf(fs, "%g %g %lg\n", r0, real_a(ctab[0][0]+ctab[1][1])/n1, uc1);

    if (uc1>ucm2) {
      ucm2=uc1;
      ind=i;
    }

    delete ctab[0];
    delete ctab;
  }
  fclose(fs);

  printf("r0 = %g; uc = %lg\n", con[sind[ind]], uc_max);

  real_a relent=0;		//total relative entropy
  real_a p1, p2;

  for (nel_ta i=0; i<n1; i++) {
    real_a ent=0;
    p1=(1-con[i])/2;
    p2=(con[i]+1)/2;
    if (p1!=0) ent=-p1*log(n1*p1/nc1);
    if (p2!=0) ent-=p2*log(n2*p2/nc2);
    //printf("p1=%g; p2=%g; e=%g\n", p1, p2, ent);
    relent+=ent;
  }
  printf("relative entropy=%g\n", relent/log(2));
  printf("uc * n =         %g\n", uc_max*n1);

  delete [] class1;
  delete [] class2;
  delete [] con;

  delete [] clsfile;
  delete [] confile;

  return exit_code;

}

