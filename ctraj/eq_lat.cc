#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include "error_codes.h"
#include "parse_command_opts.h"

#include "time_class.h"
#include "peteys_tmpl_lib.h"

#include "coordtran.h"
#include "tracer_anal.h"
#include "ctraj_defaults.h"
#include "av.h"

#define MAXN 10000
#define MAXLL 200
#define TLEN 30

using namespace libpetey;
using namespace libsparse;
using namespace ctraj;

int main(int argc, char **argv) {

  FILE *fs;
  FILE *fs2;
  char *tfile;
  char *outfile;
  char c;
  int32_t np;

  long nall;
  long fsize;

  float *q;
  float *eqlat;

  void *optarg[20];
  int flag[20];

  //parse the command line arguments:
  argc=parse_command_opts(argc, argv, "?", "%", optarg, flag, OPT_WHITESPACE);

  if (argc<3 || flag[0]) {
    printf("\n");
    printf("usage: eq_lat [-?] tracer eqlat \n");
    printf("\n");
    printf("where:\n");
    printf("  tracer       = binary file containing tracer field\n");
    printf("  eqlat        = binary file containing equivalent latitudes\n");
    printf("\n");
    if (flag[0]) return 0; else return -1;
  }

  tfile=argv[1];
  outfile=argv[2];

  //read in the array of tracers:
  fs=fopen(tfile, "r");
  fread(&np, sizeof(np), 1, fs);
  fseek(fs, 0, SEEK_END);
  fsize=ftell(fs);

  if ((fsize-sizeof(np)) % (sizeof(float)*np) != 0) {
    fprintf(stderr, "Warning: not an even number of tracer fields in file, %s\n", tfile);
  }

  nall=(fsize-sizeof(np))/sizeof(float)/np;

  fseek(fs, sizeof(np), SEEK_SET);

  fprintf(stderr, "%d tracer fields found in %s\n", tfile);

  fclose(fs);

  q=new float[np];
  eqlat=new float[np];
  fs2=fopen(outfile, "w");
  fwrite(&np, sizeof(np), 1, fs2);
  for (long i=0; i<nall; i++) {
    fread(q, sizeof(float), nall, fs);
    printf("%g\n", eq_lat(q, np, eqlat));
    fwrite(eqlat, sizeof(float), nall, fs2);
  }

  fclose(fs);
  fclose(fs2);

  delete [] q;
  delete [] eqlat;

  return 0;
}

