#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <assert.h>

#include "parse_command_opts.h"

#include "ctraj_defaults.h"
#include "contour_anal.h"
#include "global_metric.h"

using namespace libpetey;
using namespace ctraj;

int main(int argc, char **argv) {
  FILE *fs;
  float *lon, *lat;
  float *v;
  int npt;		//total number of points
  int nf;		//number of files
  int k;

  float s;

  void *optarg[20];
  int optflag[20];

  int32_t nvar;
  int wrapflag;

  argc=parse_command_opts(argc, argv, "o?", "%%", optarg, optflag, OPT_WHITESPACE);

  if (optflag[1]) {
    fprintf(stdout, "syntax: contour_length [-o] [file1 [file 2...]]\n");
    ctraj_optargs(stdout, "o?");
    return 0;
  }
  wrapflag=optflag[0];

  if (argc==1) {
    npt=read_ascii_contour(stdin, lon, lat);
    s=0;
    for (int i=1; i<npt; i++) {
      s+=sqrt(sdist(lon[i-1], lat[i-1], lon[i], lat[i]));
    }
    if (wrapflag==0) s+=sqrt(sdist(lon[0], lat[0], lon[npt-1], lat[npt-1]));
    delete [] lon;
    delete [] lat;
  } else {
    nf=argc-1;
    s=0;
    for (int i=0; i<nf; i++) {
      fs=fopen(argv[i+1], "r");
      fread(&nvar, sizeof(int32_t), 1, fs);
      assert(nvar==2);
      fseek(fs, 0, SEEK_END);
      npt=(ftell(fs)-sizeof(nvar))/nvar/sizeof(float);
      v=new float[npt*2];
      fseek(fs, sizeof(nvar), SEEK_SET);
      fread(v, sizeof(float), npt*2, fs);
      fclose(fs);
      for (int i=1; i<npt; i++) {
        s+=sqrt(sdist(v[(i-1)*2], v[2*i-1], v[2*i], v[2*i+1]));
      }
      delete [] v;
    }
  }

  printf("%f\n", s);

}

