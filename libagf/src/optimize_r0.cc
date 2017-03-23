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

//accuracy as a function of the discrimination border:
void acc_vs_r0(real_a r0, void *param, real_a *ac) {
  void **p2=(void **) param;
  nel_ta n=*(nel_ta *) p2[0];		//number of classes
  cls_ta *truth=(cls_ta *) p2[1];	//"true" classes
  real_a *r=(real_a *) p2[2];		//retrieved conditional prob.
  cls_ta *ret=(cls_ta *) p2[3];		//classes derived from cond. prob. & r0
  nel_ta ntrue=0;			//number of true retrievals

  for (nel_ta i=0; i<n; i++) {
    if (r[i]>r0) ret[i]=1; else ret[i]=0;
    if (ret[i]==truth[i]) ntrue++;
  }
  
  *ac=(real_a) ntrue/(real_a) n;

}
  
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

real_a uc_bin(real_a ntn,		//number of true negatives
		real_a nfp,		// " false positives
		real_a nfn,		// " false negatives
		real_a ntp) {		// " true positives
  double uc, ucr, ucs;
  cls_ta **ctab=new cls_ta*[3];
  ctab[0]=new cls_ta[9];
  ctab[1]=ctab[0]+3;
  ctab[2]=ctab[0]+6;
  
  ctab[0][0]=ntn;
  ctab[0][1]=nfp;
  ctab[1][0]=nfn;
  ctab[1][1]=ntp;
  ctab[0][2]=ntn+nfp;
  ctab[1][2]=nfn+ntp;
  ctab[2][0]=ntn+nfn;
  ctab[2][1]=nfp+ntp;
  ctab[2][2]=ntn+nfp+nfn+ntp;
  //print_contingency_table(ctab, 2, 2);
  //printf("\n");
  uc=uncertainty_coefficient(ctab, 2, 2, ucr, ucs);

  delete [] ctab[0];
  delete [] ctab;

  return uc;

}
//skill as a function of the discrimination border:
void skill_vs_r0(real_a r0, void *param, real_a *uc) {
  nel_ta ntp, ntn, nfp, nfn;		//confusion matrix
  double ind;				//location of r0
  void **p2=(void **) param;
  nel_ta n=*(nel_ta *) p2[0];		//number of classes
  real_a *r=(real_a *) p2[1];		//sorted retrieved conditional prob. (not sorted)
  cls_ta *none=(nel_ta *) p2[2];	//cumulative number of ones
  real_a (*sfunc) (real_a, real_a, real_a, real_a)=(real_a (*) (real_a, real_a, real_a, real_a)) p2[3];

  printf("interpolating r stuff\n");
  ind=interpolate(r, (long) n, r0);
  printf("index=%lg\n", ind);

  nfn=none[(int) ind];
  ntn=ind-nfn;
  ntp=none[n-1]-nfn;
  nfp=n-ind-ntp;
  //correction based on comparison with more conventional method:
  ntp--;
  ntn++;

  //*uc=-(*sfunc) (ntn, nfp, nfn, ntp);
  *uc=-uc_bin(ntn, nfp, nfn, ntp);
}

//compares two sets of class and outputs statistics on them:
int main(int argc, char ** argv) {
  //function to minimize:
  void (*func) (real_a, void *, real_a *)=&acc_vs_r0;
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
  int nhist=50;
  int type=0;

  cls_ta nc1, nc2;		//numbers of each class

  int exit_code=0;

  char c;

  while ((c = getopt(argc, argv, "bq:Q:")) != -1) {
    switch (c) {
      case ('b'):
             func=&uc_vs_r0;
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
    printf("usage:  optimize_r0 [-q samples] [-b] file1 file2\n\n");
    printf("where:\n");
    printf("  file1   = binary file containing true class values\n");
    printf("  file2   = pair of binary files containing corresponding decision functions:\n");
    printf("              .cls for classification results\n");
    printf("              .con for confidence ratings\n\n");
    printf("  -b      = optimize based on uncertainty coefficient\n");
    printf("  -q samp = number of samples for integrating ROC curve\n");
    printf("  -Q type = type of graph to print out:\n");
    printf("              0 print no graph\n");
    printf("              1 ROC curve\n");
    printf("              2 r cumulative distribution\n");
    printf("              3 accuracy vs. r_0\n");
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

  //fill up parameters to pass to minimization routine:
  /*
  param[0]=&n1;
  param[1]=class1;
  param[2]=con;
  param[3]=class2;

  r0=min_golden(&uc_vs_r0, (void *) param, (real_a) -1.0, (real_a) 0.0, 
		(real_a) 1.0, (real_a) 1e-7, maxiter, niter, uc_max);

  uc_max=-uc_max;

  printf("r0 = %g; uc = %g\n", r0, uc_max);
  */

  nel_ta ind=0;			//index of maximum skill
  long *sind=heapsort(con, n1);
  real_a rsort[n1];
  //real_a *rsort=map_vector(con, sind, n1);
  //ROC curve:
  real_a hitrate[n1];
  real_a farate[n1];
  nel_ta nonecum[n1];		//cumulative number of ones
  fs=stdout;
  //nonecum[0]=class1[0];
  nonecum[0]=class1[sind[0]];
  for (nel_ta i=1; i<n1; i++) nonecum[i]=nonecum[i-1]+class1[sind[i]];
  for (nel_ta i=0; i<n1; i++) {
    hitrate[i]=(real_a) (nonecum[n1-1]-nonecum[i])/(real_a) nonecum[n1-1];
    farate[i]=(real_a) (n1-i+nonecum[i]-nonecum[n1])/(real_a) (n1-nonecum[n1-1]);
  }
  printf("sorting prob.\n");
  for (nel_ta i=0; i<n1; i++) rsort[i]=con[sind[i]];

  param[0]=&n1;
  param[1]=rsort;
  param[2]=nonecum;
  param[3]=(void *) &uc_bin;

  printf("maximizing skill\n");
  r0=min_golden(&skill_vs_r0, (void *) param, (real_a) -1.0, (real_a) 0.0, 
		(real_a) 1.0, (real_a) 1e-7, maxiter, niter, uc_max);

  uc_max=-uc_max;

  printf("r0 = %g; uc = %g\n", r0, uc_max);

  //integrate ROC curve:
  real_a thetaold=M_PI/2;
  real_a Rold=sqrt(hitrate[0]*hitrate[0]+(1-farate[0])*(1-farate[0]));
  real_a rocarea=0;
  for (nel_ta i=1; i<n1; i++) {
    real_a theta=atan(hitrate[i]/(1-farate[i]));
    real_a Rcur=sqrt(hitrate[i]*hitrate[i]+(1-farate[i])*(1-farate[i]));
    rocarea+=(thetaold-theta)*(Rold+Rcur)*(Rold+Rcur)/4;
    thetaold=theta;
    Rold=Rcur;
  }

  printf("area = %g\n", rocarea/2);

  //maximize uncertainty coefficient just based on sorting confidence ratings:
  real_a uc1;
  double ucm2=0;
  cls_ta **ctab=new cls_ta*[3];
  ctab[0]=new cls_ta[9];
  ctab[1]=ctab[0]+3;
  ctab[2]=ctab[0]+6;
  for (int i=0; i<n1; i++) {
    double ucr, ucs, uc1;
    cls_ta nclt, nclr;
    //number of true positives, true negatives, false positives, false negatives:
    nel_ta ntp, ntn, nfp, nfn;
    real_a hi, hj, hij;
    real_a uc2;

    ctab=new cls_ta*[3];
    ctab[0]=new cls_ta[9];
    ctab[1]=ctab[0]+3;
    ctab[2]=ctab[0]+6;
    nclt=2;
    nclr=2;

    nfn=nonecum[i];
    ntn=i-nfn;
    ntp=nonecum[n1-1]-nfn;
    nfp=n1-i-ntp;

    ctab[0][0]=ntn;
    ctab[0][1]=nfp;
    ctab[1][0]=nfn;
    ctab[1][1]=ntp;
    ctab[0][2]=ntn+nfp;
    ctab[1][2]=nfn+ntp;
    ctab[2][0]=ntn+nfn;
    ctab[2][1]=nfp+ntp;
    ctab[2][2]=ntn+nfp+nfn+ntp;
    //print_contingency_table(ctab, nclt, nclr);
    //printf("\n");
    uc1=uncertainty_coefficient(ctab, nclt, nclr, ucr, ucs);

    if (uc1>ucm2) {
      ucm2=uc1;
      ind=i;
    }
  }

  delete ctab[0];
  delete ctab;

  if (type==1) {
    for (nel_ta i=0; i<n1; i++) {
      fprintf(fs, "%g %g\n", farate[i], hitrate[i]);
    }
  } else if (type==2) {
    fs=stdout;
    for (int i=0; i<n1; i++) {
      fprintf(fs, "%d %g\n", i, con[sind[i]]);
    }
    //fclose(fs);
  } else if (type==3) {
    fs=fopen("accvsr0.txt", "w");
    for (int i=0; i<n1; i++) {
      double ucr, ucs, uc1;
      cls_ta **ctab;
      cls_ta nclt, nclr;
      //number of true positives, true negatives, false positives, false negatives:
      nel_ta ntp, ntn, nfp, nfn;
      real_a uc2;

      r0=con[sind[i]];
      for (int j=0; j<n1; j++) if (con[j]<=r0) class2[j]=0; else class2[j]=1;
      //very wasteful since we can use the stuff before to eliminate much of
      //the computation:
      /*
      ctab=new cls_ta*[3];
      ctab[0]=new cls_ta[9];
      ctab[1]=ctab[0]+3;
      ctab[2]=ctab[0]+6;
      nclt=2;
      nclr=2;
      */
      ctab=build_contingency_table(class1, class2, n1, nclt, nclr);
      uc1=uncertainty_coefficient(ctab, nclt, nclr, ucr, ucs);
      print_contingency_table(ctab, nclt, nclr);

      nfn=nonecum[i];
      ntn=i-nfn;
      ntp=nonecum[n1-1]-nfn;
      nfp=n1-i-ntp;

      ctab[0][0]=ntn;
      ctab[0][1]=nfp;
      ctab[1][0]=nfn;
      ctab[1][1]=ntp;
      ctab[0][2]=ntn+nfp;
      ctab[1][2]=nfn+ntp;
      ctab[2][0]=ntn+nfn;
      ctab[2][1]=nfp+ntp;
      ctab[2][2]=ntn+nfp+nfn+ntp;
      print_contingency_table(ctab, nclt, nclr);
      printf("\n");
      uc2=uncertainty_coefficient(ctab, nclt, nclr, ucr, ucs);

      //fprintf(fs, "%g", r0);
      //class_eval_basic(class1, class2, n1, fs);
      //fprintf(fs, "%g %g %lg\n", r0, real_a(ctab[0][0]+ctab[1][1])/n1, uc1);
      fprintf(fs, "%g %g %lg\n", r0, real_a(ntp+ntn)/n1, uc2);

      delete ctab[0];
      delete ctab;
    }
    fclose(fs);
    printf("r0 = %g; uc = %lg\n", con[sind[ind]], ucm2);
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
  printf("uc * n =         %g\n", uc_max*n1);

  delete [] class1;
  delete [] class2;
  delete [] con;

  delete [] clsfile;
  delete [] confile;

  return exit_code;

}

