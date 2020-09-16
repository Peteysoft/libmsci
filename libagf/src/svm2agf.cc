#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <assert.h>

#include "error_codes.h"

#include "agf_lib.h"

using namespace std;
using namespace libpetey;
using namespace libagf;

int main(int argc, char **argv) {
  char *outfile;
  char *vecfile;
  char *clsfile;

  FILE *ifs;
  FILE *ofs1, *ofs2;

  dim_ta nvar;
  long n;
  int dum;
  real_a **vec;			//feature data
  cls_ta *cls;			//class data

  real_a missing=0.;		//need to be able to set this...

  char c;
  int hflag=0;
  int Uflag=0;
  int Lflag=0;

  int errcode=0;

  while ((c = getopt(argc, argv, "LUHE:")) != -1) {
    switch (c) {
      case ('L'):
             Lflag=1;
	     break;
      case ('U'):
             Uflag=1;
	     break;
      case ('H'):
             hflag=1;
	     break;
      case ('E'):
             if (sscanf(optarg, "%g", &missing)!=1) {
               fprintf(stderr, "svm2agf: error parsing argument -E\n");
               exit(FATAL_COMMAND_OPTION_PARSE_ERROR);
             }
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

  if (argc < 1) {
    printf("Converts class training data in LIBSVM compatible ascii format\n");
    printf("to libAGF compatible binary format.\n");
    printf("* Note: will not convert multi-label formats.\n\n");
    printf("usage: svm2agf2 [-H] [-R] [-U] [-E missing] [infile] basename\n\n");
    printf("    where:\n");
    printf("infile   = name of ASCII input file\n");
    printf("basename = base name of binary output files:\n");
    printf("           .vec for vector data, .cls for class data\n");
    printf("missing  = value for missing feature data\n");
    printf("-H       = omit header in output vector file\n");
    printf("-L       = input file contains floating point ordinates.\n");
    printf("-U       = re-label classes so that they go from [0-nc).\n");
    exit(0);
  }

  if (argc==1) {
    ifs=stdin;
    outfile=argv[0];
  } else {
    ifs=fopen(argv[0], "r");
    outfile=argv[1];
  }

  n=fscanf(ifs, "%d", &dum);
  rewind(ifs);
  if (Lflag) {
    real_a *ord;
    n=read_svm(ifs, vec, ord, nvar, missing);
    clsfile=new char[strlen(outfile)+5];
    sprintf(clsfile, "%s.dat", outfile);
    ofs2=fopen(clsfile, "w");
    fwrite(ord, sizeof(real_a), n, ofs2);
    fclose(ofs2);
    delete [] ord;
  } else {
    n=read_svm(ifs, vec, cls, nvar, missing, Uflag);
    clsfile=new char[strlen(outfile)+5];
    sprintf(clsfile, "%s.cls", outfile);
    ofs2=fopen(clsfile, "w");
    fwrite(cls, sizeof(cls_ta), n, ofs2);
    fclose(ofs2);
    delete [] cls;
  }

  if (n<1) {
    fprintf(stderr, "svm2agf: an error occurred\n");
    exit(FILE_READ_ERROR);
  }

  vecfile=new char[strlen(outfile)+5];
  sprintf(vecfile, "%s.vec", outfile);

  ofs1=fopen(vecfile, "w");
  if (hflag==0) fwrite(&nvar, sizeof(nvar), 1, ofs1);
  fwrite(vec[0], sizeof(real_a), nvar*n, ofs1);
  fclose(ofs1);

  delete [] vecfile;
  delete [] clsfile;

  delete [] vec[0];
  delete [] vec;

  return errcode;

}

