#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <getopt.h>

#include "agf_lib.h"

#include "sample_class1_obj.h"
#include "sample_class2_obj.h"

using namespace std;
using namespace libagf;

int main(int argc, char **argv) {
  sample_class1_obj<real_a> sc1;
  sample_class2_obj<real_a> sc2;

  FILE *fs;
  char *clsfile;
  char *confile;

  real_a **test;		//test data
  dim_ta d;		//must be two...
  nel_ta n;		//number of points
  real_a ratio;		//as in doc
  real_a p1, p2;		//pdfs of first and second class resp.
  real_a pt;
  cls_ta *cls;		//classes
  real_a *con;		//confidences

  int jflag=0;		//print joint instead of conditional probabilities

  int exit_code=0;

  char c;

  while ((c=getopt(argc, argv, "j")) != -1) {
    switch (c) {
      case ('j'):
        jflag=1;
        break;
      case ('?'):
        fprintf(stderr, "Unknown option: -%c -- ignored\n", optopt);
        //kind of stupid since it comes out as 0 anyway:
        exit_code=COMMAND_OPTION_PARSE_ERROR;
        break;
      default:
        fprintf(stderr, "Error parsing command line\n");
        exit(FATAL_COMMAND_OPTION_PARSE_ERROR);
    }
  }

  argc-=optind;
  argv+=optind;


  if (argc < 2) {
    printf("\nDetermines the most likely class of a set of test points based on\n");
    printf("analytically and semi-analytically computed pdfs of the synthetic\n");
    printf("test classes.  In theory this should provide the best possible\n");
    printf("result for this pair of classes\n\n");
    printf("usage:  classify_sc testfile outfile [ratio]\n\n");
    printf("where:\n");
    printf("testfile = binary file containing test points\n");
    printf("outfile  = base name of the binary output files\n");
    printf("           (.cls for classes .con for confidence ratings)\n");
    printf("ratio    = ratio of the size of second class to the first\n");
    printf("           (default is two (2))\n\n");
    exit(3-argc);
  }

  test=read_vecfile(argv[0], n, d);
  if (d != 2) {
    fprintf(stderr, "Error: classes are only defined in two dimensions!\n");
    fprintf(stderr, "Test data has %d.\n", d);
    exit(2-d);
  }

  if (argc > 2) sscanf(argv[2], "%f", &ratio); else ratio=2.;

  cls=new cls_ta[n];
  con=new real_a[n];

  for (nel_ta i=0; i<n; i++) {
    p1=sc1.pdf(test[i][0], test[i][1])/(1.+ratio);
    p2=sc2.pdf(test[i][0], test[i][1])*ratio/(1.+ratio);
    if (p1 > p2) cls[i]=0; else cls[i]=1;
    pt=p1+p2;
    con[i]=fabs(p1-p2)/pt;
    //printf("%8d %4d   %8.6f %8.6f\n", i, cls[i], p1/pt, p2/pt);
    if (jflag) {
      printf("%8.6f %8.6f %4d\n", p1, p2, cls[i]);
    } else {
      printf("%8.6f %8.6f %4d\n", p1/pt, p2/pt, cls[i]);
    }
  }

  clsfile=new char[strlen(argv[1])+5];
  confile=new char[strlen(argv[1])+5];

  sprintf(clsfile, "%s.cls", argv[1]);
  sprintf(confile, "%s.con", argv[1]);

  fs=fopen(clsfile, "w");
  fwrite(cls, sizeof(cls_ta), n, fs);
  fclose(fs);

  fs=fopen(confile, "w");
  fwrite(con, sizeof(real_a), n, fs);
  fclose(fs);

  delete [] test[0];
  delete [] test;
  delete [] cls;
  delete [] con;

  delete [] clsfile;
  delete [] confile;

  return exit_code;

}

