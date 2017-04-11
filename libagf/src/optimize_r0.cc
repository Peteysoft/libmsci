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

//compares two sets of class and outputs statistics on them:
int main(int argc, char ** argv) {
  real_a (*skill_func) (real_a, real_a, real_a, real_a)=&acc_bin;
  FILE *fs;

  char *tfile;
  char *clsfile;
  char *confile;
  size_t slen;
  char *outfile=NULL;

  cls_ta *class1;
  cls_ta *class2;
  nel_ta n1, n2;

  real_a *con;

  real_a max_skill;		//maximum uncertainty coeff.
  real_a r0;			//optimal threshold
  real_a area;			//area under ROC curve

  long maxiter=1000;
  long niter;
  int nhist=50;
  int type=0;			//what to do: see help screen
  int skill_type=0;		//type of skill score: 2 accuracy; 3 U.C.; 4 Pearson corr.
  char skill_name[10]="accuracy";

  cls_ta nc1, nc2;		//numbers of each class

  int exit_code=0;

  char c;

  while ((c = getopt(argc, argv, "bq:Q:c:")) != -1) {
    switch (c) {
      case ('c'):
	     skill_type=atoi(optarg);
	     break;
      case ('q'):
	     nhist=atoi(optarg);
	     break;
      case ('Q'):
	     type=atoi(optarg);
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
    printf("\nAnd calculate the area under the ROC curve.\n");
    printf("\n\n");
    printf("usage:  optimize_r0 [-q samples] [-b] file1 file2 [output]\n\n");
    printf("where:\n");
    printf("  file1    = binary file containing true class values\n");
    printf("  file2    = pair of binary files containing corresponding decision functions:\n");
    printf("               .cls for classification results\n");
    printf("               .con for confidence ratings\n\n");
    printf("  output   = where to output graphs (if applicable)\n");
    printf("  -q samp  = number of samples for integrating ROC curve\n");
    printf("  -Q type  = type of graph to print out:\n");
    printf("               0 print no graph\n");
    printf("               1 ROC curve\n");
    printf("               2 r cumulative distribution\n");
    printf("               3 skill vs. r_0\n");
    printf("  -c skill = type of skill to optimize:\n");
    printf("               2 accuracy\n");
    printf("               3 uncertainty coefficient\n");
    printf("               4 correlation\n");
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

  if (argc > 2) outfile=argv[2];

  //read in "true" classes:
  class1=read_clsfile<cls_ta>(tfile, n1);
  if (n1 < 0) {
    fprintf(stderr, "Error reading input file: %s\n", tfile);
    return ALLOCATION_FAILURE;
  }
  if (class1 == NULL) {
    fprintf(stderr, "Unable to open file for reading: %s\n", tfile);
    return UNABLE_TO_OPEN_FILE_FOR_READING;
  }
  
  //read in retrieved classes:
  class2=read_clsfile<cls_ta>(clsfile, n2);
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
  con=read_datfile<real_a>(confile, n2);
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

  long *sind=heapsort(con, n1);
  real_a rsort[n1];
  nel_ta nonecum[n1];		//cumulative number of ones
  real_a *skill;

  //prepare for rest of work:
  sortr_cumulate_ones(class1, con, n1, nonecum, rsort);
  area=integrate_roc(nonecum, rsort, n1);
  printf("ROC area = %g\n", area);

  switch(skill_type) {
    case(2):
      skill_func=&acc_bin;
      break;
    case(3):
      skill_func=&uc_bin;
      strcpy(skill_name, "U.C.");
      break;
    case(4):
      skill_func=&corr_bin;
      strcpy(skill_name, "corr.");
      break;
  }
  skill=calc_skill_vs_r0(nonecum, rsort, n1, skill_func);
  nel_ta i0=0;
  while (isfinite(skill[i0])!=1) {i0++;};
  max_skill=skill[i0];
  for (nel_ta i=i0+1; i<n1; i++) {
    if (skill[i]>max_skill) {
      max_skill=skill[i];
      r0=rsort[i];
    }
  }
  printf("r0 = %g\n", r0);
  printf("max %s = %g\n", skill_name, max_skill);

  if (outfile!=NULL) fs=fopen(outfile, "w"); else fs=stdout;
  if (type==1) {
    //ROC curve:
    real_a **roc;
    roc=calc_roc(nonecum, rsort, n1);
    for (nel_ta i=0; i<n1; i++) {
      fprintf(fs, "%g %g\n", roc[0][i], roc[1][i]);
    }
  } else if (type==2) {
    for (int i=0; i<n1; i++) {
      fprintf(fs, "%d %g\n", i, rsort[i]);
    }
    //fclose(fs);
  } else if (type==3) {
    for (int i=0; i<n1; i++) {
      fprintf(fs, "%g %g\n", rsort[i], skill[i]);
    }
  }

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
  printf("%s * n =         %g\n", skill_name, max_skill*n1);

  fclose(fs);
  delete [] class1;
  delete [] class2;
  delete [] con;

  delete [] clsfile;
  delete [] confile;

  delete [] skill;

  return exit_code;

}

