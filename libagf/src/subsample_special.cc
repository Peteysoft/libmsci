
//
// For subsampling keeping relative class ratios constant or making sure 
// smallest classes have some representation for very uneven class 
// distributions.
//
// Should be merged with agf_preprocess.
//

#include <math.h>
#include <string.h>
#include <stdio.h>

#include <gsl/gsl_linalg.h>

#include "randomize.h"
#include "peteys_tmpl_lib.h"
#include "full_util.h"
#include "roots_mins.h"
#include "agf_lib.h"

using namespace std;
using namespace libagf;
using namespace libpetey;

void ffromk(real_a k, void *param, real_a *f) {
  void **p2=(void **) param;
  cls_ta ncls=*(cls_ta *) p2[0];
  nel_ta *cind=(nel_ta *) p2[1];
  real_a C=pow(*(nel_ta *) p2[2], k);
  real_a f0=*(real_a *) p2[3];
  real_a sum=0;


  for (cls_ta i=0; i<ncls; i++) {
    nel_ta ni=cind[i+1]-cind[i];
    sum+=C*pow(ni, -k)*ni;
  }
  //printf("C=%g; sum=%g; nt=%d\n", C, sum, cind[ncls]);

  *f=sum/cind[ncls]-f0;
}

int main(int argc, char *argv[]) {
  char *vecfile=NULL;		//training data
  char *clsfile=NULL;		//class data
  char *outbase;
  char *outvec=NULL;			//output file
  char *outcls=NULL;
  char *testvec=NULL;
  char *testcls=NULL;
  FILE *fs;
  FILE *diagfs;			//print diagnostics to this file stream

  nel_ta ntrain;		//number of training data points
  nel_ta ntrain2;
  dim_ta nvar;		//number of variables
  cls_ta ncls;		//number of classes

  real_a **train;	//training data vectors
  real_a *all;		//train[0] to this for deletion
  cls_ta *ord=NULL;
  size_t ordsize;	//size of ordinates

  int exit_value;

  agf_command_opts opt_args;

  exit_value=0;

  //parse the command line arguments:
  //exit_value=agf_parse_command_opts(argc, argv, "c:d:f:AMRz", &opt_args);
  exit_value=agf_parse_command_opts(argc, argv, "d:f:CRz", &opt_args);
  if (exit_value==FATAL_COMMAND_OPTION_PARSE_ERROR) return exit_value;

  //print help page if there are not enough mandatory arguments:
  if (argc<3) {
    printf("\n");
    printf("Syntax:	subsample_special [-d ndiv] [-f frac] \n");
    printf("                          input output test\n");
    printf("\n");
    printf("arguments:\n");
    printf("  input        base name for binary input files\n");
    printf("  output       base name for binary output files:\n");
    printf("                 .vec for coordinate data (features)\n");
    printf("                 .cls for class data\n");
    printf("  test         base name for test data\n");
    printf("\n");
    //printf("file options:\n");
    //printf("  -A           operate on ASCII files\n");
    //printf("  -M           specifies LIBSVM format\n");
    //printf("\n");
    printf("  -C           keep relative class numbers constant\n");
    printf("  -z           randomly permute data\n");
    printf("  -R           data separation works by random selection rather than\n");
    printf("                 permutation (output files are not of definite size)\n");
    printf("  -d ndiv      number of separate output files\n");
    printf("  -f frac      separate into test and training (over-rides -d)\n");
    printf("\n");
    printf("\n");
    return INSUFFICIENT_COMMAND_ARGS;
  }

  diagfs=stderr;
  ran_init();

  if (opt_args.asciiflag) {
    vecfile=new char[strlen(argv[0])+1];
    sprintf(vecfile, "%s", argv[0]);
  } else {
    vecfile=new char[strlen(argv[0])+5];
    sprintf(vecfile, "%s.vec", argv[0]);
  }

  fs=fopen(vecfile, "r");
  if (fs==NULL) {
    fprintf(stderr, "agf_preprocess: Unable to open coordinate data for reading\n");
    exit(UNABLE_TO_OPEN_FILE_FOR_READING);
  }

  //read in training data and set the output file names:
  //fprintf(stderr, "agf_preprocess: Reading training data\n");
  outbase=argv[1];
  if (opt_args.asciiflag) {
    if (opt_args.Mflag) {
      ntrain=read_svm(fs, train, ord, nvar);
    } else {
      ntrain=read_lvq(fs, train, ord, nvar);
    }
    fclose(fs);
  } else {
    int32_t nt1, nv1;

    train=read_matrix<real_a, int32_t>(fs, nt1, nv1);
    if (nt1 == -1) {
      fprintf(stderr, "agf_preprocess: Error reading coordinate data\n");
      exit(FILE_READ_ERROR);
    }
    ntrain=nt1;
    nvar=nv1;
    fclose(fs);

    outvec=new char[strlen(outbase)+5];
    sprintf(outvec, "%s.vec", outbase);

    clsfile=new char[strlen(argv[0])+5];
    if (opt_args.Lflag) {
      sprintf(clsfile, "%s.dat", argv[0]);
    } else {
      sprintf(clsfile, "%s.cls", argv[0]);
    }

    outcls=new char[strlen(outbase)+5];
    if (opt_args.Lflag) {
      sprintf(outcls, "%s.dat", outbase);
    } else {
      sprintf(outcls, "%s.cls", outbase);
    }

    ord=read_clsfile<cls_ta>(clsfile, ntrain2);
    if (ntrain2 == -1) {
      fprintf(stderr, "agf_preprocess: Error reading file: %s\n", clsfile);
      exit(FILE_READ_ERROR);
    }
    if (ord==NULL) {
        fprintf(stderr, "agf_preprocess: Unable to open file for reading: %s\n", clsfile);
        exit(UNABLE_TO_OPEN_FILE_FOR_READING);
    }
    fprintf(diagfs, "%d class labels found in file: %s\n", ntrain2, clsfile);
    if (ntrain2!=ntrain) {
      fprintf(stderr, "agf_preprocess: Sample count mismatch\n");
      exit(SAMPLE_COUNT_MISMATCH);
    }
  }

  testvec=new char[strlen(argv[2])+5];
  sprintf(testvec, "%s.vec", argv[2]);
  testcls=new char[strlen(argv[2])+5];
  if (opt_args.Lflag) {
    sprintf(testcls, "%s.dat", argv[2]);
  } else {
    sprintf(testcls, "%s.cls", argv[2]);
  }

  all=train[0];		//need this when it's time to clean up

  //randomly permute the data:
  if (opt_args.zflag) {
    //fprintf(stderr, "agf_preprocess: randomizing data\n");
    real_a **tnew=new real_a *[ntrain];
    long *rind=randomize(ntrain);
    for (nel_ta i=0; i<ntrain; i++) {
      tnew[i]=train[rind[i]];
    }
    delete [] train;
    train=tnew;
    cls_ta *cnew=new cls_ta[ntrain];
    for (nel_ta i=0; i<ntrain; i++) {
      cnew[i]=ord[rind[i]];
    }
    ord=cnew;
    delete [] rind;
    //randomize_vec(train, nvar, ntrain, cls);
  }

  //fprintf(stderr, "agf_preprocess: sorting ordinates\n");
  ncls=1;
  for (int i=0; i<ntrain; i++) if (ord[i]>=ncls) ncls=ord[i]+1;
  nel_ta *cind=sort_classes(train, ntrain, ord, ncls);
  fprintf(diagfs, "Class locations:\n");
  nel_ta n0=ntrain;
  nel_ta sumnc2=0;		//sum of class numbers squared
  nel_ta nmax=0;
  for (cls_ta i=0; i<ncls; i++) {
    fprintf(diagfs, "%d: %d %d\n", i, cind[i], cind[i+1]-cind[i]);
    if (cind[i+1]-cind[i]<n0) {
      n0=cind[i+1]-cind[i];
    }
    if (cind[i+1]-cind[i]>nmax) nmax=cind[i+1]-cind[i];
    sumnc2+=(cind[i+1]-cind[i])*(cind[i+1]-cind[i]);
  }

  FILE *trainvecfs;
  FILE *trainclsfs;
  FILE *testvecfs;
  FILE *testclsfs;
  trainvecfs=fopen(outvec, "w");
  testvecfs=fopen(testvec, "w");
  fwrite(&nvar, sizeof(nvar), 1, trainvecfs);
  fwrite(&nvar, sizeof(nvar), 1, testvecfs);
  trainclsfs=fopen(outcls, "w");
  testclsfs=fopen(testcls, "w");

  //find the parameter to get the desired fraction:
  real_a C;
  real_a fmax;
  if (opt_args.Cflag!=1) {
    void *param[4];
    long niter;
    param[0]=&ncls;
    param[1]=cind;
    param[2]=&n0;
    param[3]=&opt_args.ftest;
    fmax=root_false_position(&ffromk, (void *) param, (real_a) 0., (real_a) 1., (real_a) 1e-6, (real_a) 1e-6, (long) 1000, niter);
    C=pow(n0, fmax);
  }
  for (cls_ta i=0; i<ncls; i++) {
    nel_ta ni=cind[i+1]-cind[i];
    size_t l1, l2;
    //l1=size_t(ntrain*(1-opt_args.ftest));	//don't do this...
    if (opt_args.Cflag) {
      l2=size_t(ni*opt_args.ftest);
    } else {
      l2=size_t(C*ni*pow(ni, -fmax));
      printf("%g: %d\n", C*pow(ni, -fmax), l2);
    }
    l1=cind[i+1]-cind[i]-l2;
    for (nel_ta j=cind[i]; j<cind[i]+l1; j++) {
      fwrite(train[j], sizeof(real_a), nvar, trainvecfs);
    }
    for (nel_ta j=cind[i]+l1; j<cind[i+1]; j++) {
      fwrite(train[j], sizeof(real_a), nvar, testvecfs);
    }
    fwrite(ord+cind[i], sizeof(cls_ta), l1, trainclsfs);
    fwrite(ord+cind[i]+l1, sizeof(cls_ta), l2, testclsfs);
  }

  fclose(trainvecfs);
  fclose(testvecfs);
  fclose(trainclsfs);
  fclose(testclsfs);

  delete [] cind;

  //clean up:
  delete [] all;
  delete [] train;

  delete [] ord;

  if (clsfile!=NULL) delete [] clsfile;
  if (outcls!=NULL) delete [] outcls;

  if (vecfile!=NULL) delete [] vecfile;
  if (outvec!=NULL) delete [] outvec;

  ran_end();

  return exit_value;

}


