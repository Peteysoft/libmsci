#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <getopt.h>

//takes a stream of data from stdin and generates a GMT-compatible
//file with contour levels in it

int main(int argc, char **argv) {
  char c;
  int ncon;
  int gflag=0;
  int bflag=0;
  int tflag=0;
  int hflag=0;
  float bottom;
  float top;
  int n;

  FILE *fs;

  FILE *help;

  n=21;

  while ((c = getopt(argc, argv, "z:I:F:gh")) != -1) {
    switch (c) {
      case ('z'):
        ncon=sscanf(optarg, "%d", &n);
        if (ncon != 1) {
          fprintf(stderr, "Warning: garbled command option argument: -z %s", optarg);
          exit(2);
        }
        break;
      case ('I'):
        ncon=sscanf(optarg, "%f", &bottom);
        if (ncon != 1) {
          fprintf(stderr, "Warning: garbled command option argument: -I %s", optarg);
          exit(2);
        }
	bflag=1;
        break;
      case ('F'):
        ncon=sscanf(optarg, "%f", &top);
        if (ncon != 1) {
          fprintf(stderr, "Warning: garbled command option argument: -F %s", optarg);
          exit(2);
        }
	tflag=1;
        break;
      case ('g'):
	     gflag=1;
	     break;
      case ('h'):
	     hflag=1;
	     break;
      case ('?'):
             fprintf(stderr, "Unknown option: -%c -- ignored\n", optopt);
             break;
      default:
             fprintf(stderr, "Error parsing command line\n");
             exit(2);
             break;
    }
  }
  argc-=optind;
  argv+=optind;

  help=stdout;

  if (hflag == 1) {
    fprintf(help, "\n");
    fprintf(help, "Syntax:  gen_zgrid [-h] [-g] [-I bottom] [-F top] [-z nz]\n");
    fprintf(help, "\n");
    fprintf(help, "Based on data from stdin,\n"); 
    fprintf(help, "generates a series of contour lines.\n");
    fprintf(help, "Contour intervals are printed to stdout\n");
    fprintf(help, "in a Generic Mapping Tools (GMT)-compatible format.\n");
    fprintf(help, "\n");
    fprintf(help, "where:\n");
    fprintf(help, "  bottom   is the bottom-most contour\n");
    fprintf(help, "  top      is the top-most contour\n");
    fprintf(help, "  nz       is the number of contour intervals [%d]\n", n);
    fprintf(help, "  -g       signifies a geometric progresson\n");
    fprintf(help, "  -h       print this help screen\n");
    fprintf(help, "\n");
    exit(1);
  }

  fs=stdin;

  long k=0;
  float x, y;

  if (bflag == 0 || tflag == 0) {
    while(feof(fs) == 0) {
      if (fscanf(fs, "%f", &x)==0) break;
      if (k==0) {
        if (bflag==0) bottom=x;
        if (tflag==0) top=x;
      }
      if (bflag==0 && x<bottom) bottom=x;
      if (tflag==0 && x>top) top=x;
      k++;
    }
  }

  if (gflag==1) {
    bottom=log(bottom);
    top=log(top);
  }

  for (int i=0; i<n; i++) {
    float z=i*(top-bottom)/(n-1)+bottom;
    if (gflag==1) z=exp(z);
    fprintf(stdout, "%g\n", z);
  }
}

