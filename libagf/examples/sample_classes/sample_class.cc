#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/timeb.h>

//#include "nr.h"
#include <gsl/gsl_rng.h>
#include "randomize.h"

#include "agf_lib.h"

#include "sample_class1_obj.h"
#include "sample_class2_obj.h"

using namespace std;
using namespace libagf;
using namespace libpetey;

int main(int argc, char ** argv) {
  gsl_rng *rann;

  char *vecfile;
  char *clsfile;

  FILE *vecfs;
  FILE *clsfs;

  sample_class1_obj<real_a> sc1;
  sample_class2_obj<real_a> sc2;

  //number of samples:
  nel_ta n1, n2;

  real_a x, y;
  dim_ta ndim=2;
  cls_ta cls;

  int rflag;

  if (argc < 4 || argc > 6) {
    printf("\nGenerates a pair of synthetic classes for validation\n");
    printf("of classification algorithms\n\n");
    printf("usage:  sample_class [-R] n1 n2 output\n\n"); 
	// [ c1file [sfile]]\n\n");
    printf("      where:\n");
    printf("-R      = randomly pick class based on ratio between the two\n");
    printf("n1      = number of samples in first class\n");
    printf("n2      = number of samples in second class\n");
    printf("output  = base name of output files\n");
    //printf("c1file  = ascii file location and size of first class\n");
    //printf("sfile   = ascii file defining 'spine' of second class\n");
    printf("\n");
    exit(1);
  }

  rflag=0;
  if (argc > 4) {
    if (argv[1][0]=='-') {
      if (argv[1][1]='R') {
        rflag=1;
      } else {
        fprintf(stderr, "Unrecognized option: %s\n", argv[1]);
      }
      argv++;
    }
  }

  sscanf(argv[1], "%d", &n1);
  sscanf(argv[2], "%d", &n2);

  vecfile=new char[strlen(argv[3])+5];
  sprintf(vecfile, "%s.vec", argv[3]);
  clsfile=new char[strlen(argv[3])+5];
  sprintf(clsfile, "%s.cls", argv[3]);

  vecfs=fopen(vecfile, "w");
  clsfs=fopen(clsfile, "w");

  fwrite(&ndim, sizeof(ndim), 1, vecfs);

  printf("%d\n", ndim);

  if (rflag) {
    nel_ta nt=n1+n2;

    rann=gsl_rng_alloc(agf_gsl_rng_type);
    gsl_rng_set(rann, seed_from_clock());

    for (nel_ta i=0; i<nt; i++) {
      if (gsl_rng_uniform(rann)*nt/n1 < 1) {
        sc1.sample(x, y);
        cls=0;
      } else {
        sc2.sample(x, y);
        cls=1;
      }
      printf("%f %f %d\n", x, y, cls);
      fwrite(&x, sizeof(x), 1, vecfs);
      fwrite(&y, sizeof(y), 1, vecfs);
      fwrite(&cls, sizeof(cls), 1, clsfs);
    }
    gsl_rng_free(rann);
  } else {
    cls=0;
    for (nel_ta i=0; i<n1; i++) {
      sc1.sample(x, y);
      printf("%f %f %d\n", x, y, cls);
      fwrite(&x, sizeof(x), 1, vecfs);
      fwrite(&y, sizeof(y), 1, vecfs);
      fwrite(&cls, sizeof(cls), 1, clsfs);
    }

    cls=1;
    for (nel_ta i=0; i<n2; i++) {
      sc2.sample(x, y);
      printf("%f %f %d\n", x, y, 1);
      fwrite(&x, sizeof(x), 1, vecfs);
      fwrite(&y, sizeof(y), 1, vecfs);
      fwrite(&cls, sizeof(cls), 1, clsfs);
    }
  }

  fclose(vecfs);
  fclose(clsfs);

  delete [] vecfile;
  delete [] clsfile;

}

