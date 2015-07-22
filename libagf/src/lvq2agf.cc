#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <getopt.h>

#include "error_codes.h"

#include "agf_lib.h"

using namespace std;
using namespace libagf;

int main(int argc, char **argv) {
  char *basename;
  char *vecfile;
  char *clsfile;

  FILE *ifs;
  FILE *ofs;

  dim_ta nvar;
  nel_ta ntrain;
  real_a **train;
  cls_ta *cls;

  int hflag, cflag, oflag;

  char c;

  int errcode=0;

  hflag=0;
  cflag=0;
  oflag=0;

  while ((c = getopt(argc, argv, "oHC")) != -1) {
    switch (c) {
      case ('H'):
             hflag=1;
	     break;
      case ('C'):
             cflag=1;
	     break;
      case ('o'):
             oflag=1;
	     break;
      case ('?'):
             fprintf(stderr, "Unknown option: -%c -- ignored\n", optopt);
	     errcode=10;
	     break;
      default:
	     fprintf(stderr, "Error parsing command line\n");
	     exit(21);
    }
  }

  argc-=optind;
  argv+=optind;

  if (argc < 1 && (cflag==0)) {
    printf("Usage: lvq2agf [-C] [-o] [-H] [lvqfile [agfbase]]\n");
    printf("\nConverts LVQPAK compatible ASCII files\n");
    printf("to libAGF compatible binary files.\n\n");
    printf("where:\n\n");
    printf("lvqfile is the input file name\n");
    printf("agfbase is the base name of the output files:\n");
    printf("           .vec for vectors, .cls for classes\n");
    printf("\n");
    printf("options:\n");
    printf("  -o     omit class data\n");
    printf("  -C     no class data\n");
    printf("  -H     no header\n");
    exit(1);
  }

  ifs=stdin;
  if (argc==1) {
    basename=argv[0];
  } else if (argc>=2) {
    ifs=fopen(argv[0], "r");
    basename=argv[1];
    if (ifs==NULL) {
      fprintf(stderr, "lvq2agf: unable to open file, %s, for input\n", argv[0]);
      exit(UNABLE_TO_OPEN_FILE_FOR_READING);
    }
  }

  ntrain=read_lvq(ifs, train, cls, nvar, hflag+2*cflag+4*oflag);
  if (ntrain<0) {
    fprintf(stderr, "lvq2agf: error reading input at line %d.\n", -ntrain);
    exit(FILE_READ_ERROR);
  }

  if (argc>0) {
    //generate the output file names:
    vecfile=new char[strlen(basename)+5];
    strcpy(vecfile, basename);
    strcat(vecfile, ".vec");

    if (cflag==0 && oflag==0) {
      clsfile=new char[strlen(basename)+5];
      strcpy(clsfile, basename);
      strcat(clsfile, ".cls");
      ofs=fopen(clsfile, "w");
      if (ofs==NULL) {
        fprintf(stderr, "lvq2agf: error opening output file, %s\n", clsfile);
        exit(UNABLE_TO_OPEN_FILE_FOR_WRITING);
      }
      if (fwrite(cls, sizeof(cls_ta), ntrain, ofs)!=ntrain) {
        fprintf(stderr, "lvq2agf: error writing vector data; %s\n", clsfile);
        exit(FILE_WRITE_ERROR);
      }
      fclose(ofs);
      delete [] clsfile;
    }
    ofs=fopen(vecfile, "w");
    if (ofs==NULL) {
      fprintf(stderr, "lvq2agf: error opening output file, %s\n", vecfile);
      exit(UNABLE_TO_OPEN_FILE_FOR_WRITING);
    }
    delete [] vecfile;
  } else {
    ofs=stdout;
  }
  if (fwrite(&nvar, sizeof(dim_ta), 1, ofs)!=1) {
    fprintf(stderr, "lvq2agf: error writing vector data; %s\n", vecfile);
    exit(FILE_WRITE_ERROR);
  }
  if (fwrite(train[0], sizeof(real_a), ntrain*nvar, ofs)!=ntrain*nvar) {
    fprintf(stderr, "lvq2agf: error writing vector data; %s\n", vecfile);
    exit(FILE_WRITE_ERROR);
  }
  if (argc>0) fclose(ofs);

  delete [] train[0];
  delete [] train;
  if (cls!=NULL) delete [] cls;

  return errcode;

}

