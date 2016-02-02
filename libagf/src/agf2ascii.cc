#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#include "getopt.h"

#include "full_util.h"
#include "agf_lib.h"

using namespace std;
using namespace libagf;
using namespace libpetey;

//#include "parse_command_opts.h"

//#include "/home/Petey/software_projects/include/class_lib.h"

int main(int argc, char ** argv) {

  char *vecfile;
  char *clsfile;

  FILE *fs;

  nel_ta n;
  dim_ta nvar;

  real_a **vec;
  cls_ta *cls;
  real_a *ord;

  real_a *ave, *std;
  agf_command_opts opts;

  int errcode=0;

/*
  int normflag=0;
  int svmflag=0;
  int hflag=0;
  int cflag=0;
*/

  //keep it stand-alone:
/*
  void *optarg[10];
  int flag[10];

  parse_command_opts(argc, argv, "nMhc", "%%%%", optarg, flag);
  normflag=flag[0];
  svmflag=flag[1];
  hflag=flag[2];
  cflag=flag[3];
*/

  agf_parse_command_opts(argc, argv, "CHLMn", &opts);

//  long *clind;
//  long ncl;
//
  
  if (argc < 1 && opts.Cflag==0) {
    printf("Converts vector and class data in libAGF-compatible binary format\n");
    printf("to LVQPAK or LIBSVM compatible ASCII format.\n\n");
    printf("usage: agf2ascii [-n] [-M] [-H] [-C] basename\n\n");
    printf("    where:\n");
    printf("-n       = option to normalise the data\n");
    printf("-M       = specifies LIBSVM format (LVQPAK is the default)\n");
    printf("-H       = do not print header\n");
    printf("-C       = do not print class data (**note: does not append file extension)\n");
    printf("-L       = floating point ordinates (.dat extension)\n");
    printf("basename = base name of binary input files:\n");
    printf("           .vec for vector data, .cls for class data\n");
    exit(1);
  }

  //generate the input file names:
  if (opts.Cflag) {
    vecfile=argv[0];
  } else {
    vecfile=new char[strlen(argv[0])+5];
    sprintf(vecfile, "%s.vec", argv[0]);

    if (opts.Lflag) {
      clsfile=new char[strlen(argv[0])+5];
      sprintf(clsfile, "%s.dat", argv[0]);
    } else {
      clsfile=new char[strlen(argv[0])+5];
      sprintf(clsfile, "%s.cls", argv[0]);
    }
  }

  //read in the training data:
  if (argc >= 1) {
    vec=read_vecfile<real_a>(vecfile, n, nvar);
  } else {
    int32_t nt1, nv1;
    vec=read_matrix<real_a, int32_t>(stdin, nt1, nv1);
    n=nt1;
    nvar=nv1;
  }
  if (nvar==-1 || n==-1) {
    fprintf(stderr, "agf2ascii: error reading file: %s\n", vecfile);
    exit(FILE_READ_ERROR);
  }
  if (vec==NULL) {
    fprintf(stderr, "agf2ascii: could not open file, %s, for reading\n", vecfile);
    exit(UNABLE_TO_OPEN_FILE_FOR_READING);
  }

  //not sure why we need this...?
  if (opts.normflag) {
    real_a diff;

    ave=new real_a[nvar];
    std=new real_a[nvar];

    for (nel_ta i=0; i<nvar; i++) {
      ave[i]=0;
      std[i]=0;
      for (nel_ta j=0; j<n; j++) ave[i]+=vec[j][i];
      ave[i]/=n;
      for (nel_ta j=0; j<n; j++) {
        diff=vec[j][i]-ave[i];
        vec[j][i]=diff;
        std[i]+=diff*diff;
      }
      std[i]=sqrt(std[i]/(n-1));
      for (nel_ta j=0; j<n; j++) vec[j][i]/=std[i];
    }

    delete [] ave;
    delete [] std;
  }

  if (opts.Cflag==0) {
    nel_ta n1;
    if (opts.Lflag) {
      ord=read_datfile<real_a>(clsfile, n1);
      if (ord==NULL) {
        fprintf(stderr, "agf2ascii: unable to open file, %s, for reading\n", clsfile);
        exit(UNABLE_TO_OPEN_FILE_FOR_READING);
      }
    } else {
      cls=read_clsfile<cls_ta>(clsfile, n1);
      if (cls==NULL) {
        fprintf(stderr, "agf2ascii: unable to open file, %s, for reading\n", clsfile);
        exit(UNABLE_TO_OPEN_FILE_FOR_READING);
      }
    }
    if (n1<0) {
      fprintf(stderr, "agf2ascii: error reading file: %s\n", clsfile);
      exit(FILE_READ_ERROR);
    }
    if (n1!=n) {
      fprintf(stderr, "agf2ascii: number of samples (=%d) in ordinate file (%s)\n", n1, clsfile);
      fprintf(stderr, "  does not agree with samples (=%d) in coordinate file (%s)\n", n, vecfile);
      exit(SAMPLE_COUNT_MISMATCH);
    }
  }

//  clind=sort_classes(vec, nvar, cls, n, ncl);

  if (argc>1) {
    fs=fopen(argv[1], "w");
  } else {
    fs=stdout;
  }

  if (opts.Mflag) {
    for (nel_ta i=0; i<n; i++) {
      if (opts.Cflag==0) {
        if (opts.Lflag) {
          fprintf(fs, "%g ", ord[i]);
        } else {
          fprintf(fs, "%d ", cls[i]);
        }
      }
      for (dim_ta j=0; j<nvar; j++) fprintf(fs, "%d:%g ", j+1, vec[i][j]);
      fprintf(fs, "\n");
    }
  } else {
    if (opts.Hflag==0) fprintf(fs, "%d\n", nvar);
    for (nel_ta i=0; i<n; i++) {
      for (dim_ta j=0; j<nvar; j++) fprintf(fs, "%g ", vec[i][j]);
      if (opts.Cflag==0) {
        if (opts.Lflag) {
          fprintf(fs, "%g", ord[i]);
        } else {
          fprintf(fs, "%d", cls[i]);
        }
      }
      fprintf(fs, "\n");
    }
  }

  if (argc>1) fclose(fs);

  delete [] vec[0];
  delete [] vec;
  if (opts.Cflag==0) {
    delete [] vecfile;
    delete [] clsfile;
    if (opts.Lflag) delete [] ord; else delete [] cls;
  }

  return errcode;

}

