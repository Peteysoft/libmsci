#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <sys/timeb.h>

#include "error_codes.h"
#include "randomize.h"
#include "parse_command_opts.h"
//#include "nr.h"
#include "meas_data.h"
#include "ctraj_defaults.h"

using namespace libpetey;
using namespace ctraj;

int main(int argc, char **argv) {
  meas_data *result;
  float lon;
  float lat;

  timeb now;

  time_class t0, tf;
  int32_t n=0;		//initializing here removes spurious errors returned by valgrind

  char tstring[30];
  int32_t i, j;

  int hemi=0;
  int32_t nt;
  float tstep;

  void *optarg[10];
  int flag[10];

  optarg[4]=&tstep;
  optarg[5]=&nt;

  t0.read_string("1948/1/1");

  argc=parse_command_opts(argc, argv, "if+-hN?", "%s%s%%%g%d%", optarg, flag, OPT_WHITESPACE);
  if (argc<0) {
    fprintf(stderr, "random_global: error parsing command line\n");
    return FATAL_COMMAND_OPTION_PARSE_ERROR;
  }
  if (flag[0]) t0.read_string((char *) optarg[0]);
  if (flag[1]) tf.read_string((char *) optarg[1]);
  if (flag[2]) hemi=1;
  if (flag[3]) hemi=-1;
  
  if (argc != 2 || flag[6]) {
    FILE *docfs;
    int err;
    if (flag[6]) {
      docfs=stdout;
      err=0;
    } else {
      docfs=stderr;
      err=INSUFFICIENT_COMMAND_ARGS;
    }
    fprintf(docfs, "Generates a set of random initial conditions.\n");
    fprintf(docfs, "\n");
    fprintf(docfs, "Syntax: random_global [-i t0] [-f tf] [-+|--] [-h dt -N nt] n\n");
    fprintf(docfs, "	where:\n");
    fprintf(docfs, "n:		number of samples\n");
    fprintf(docfs, "\n");
    fprintf(docfs, "options:\n");
    ctraj_optargs(docfs, "if+-hN?", 1);
    return err;
  }

  if (flag[1]==0) {
    tf="1970/1/1";
    ftime(&now);
    tf.add(now.time/86400.+now.millitm/86400000.);
  }

  //n is properly initialized here (don't listen to valgrind...):
  sscanf(argv[1], "%d", &n);

  ran_init();
  result=generate_random_global(t0, tf, n, hemi);
  ran_end();

  for (i=0; i<n; i++) {
    result[i].t.write_string(tstring);
    if (flag[4] && flag[5]) {
      printf("%23s%10.3f%10.3f%9.3f% d\n", tstring, result[i].lon, result[i].lat, tstep, nt);
    } else {
      printf("%23s%10.3f%10.3f\n", tstring, result[i].lon, result[i].lat);
    }
  }

}

