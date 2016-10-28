#include <getopt.h>
#include <stdlib.h>

#include "parse_command_opts.h"
#include "ctraj_defaults.h"

#include "meas_data.h"
#include "error_codes.h"

using namespace libpetey;
using namespace ctraj;

int main(int argc, char **argv) {
  FILE *infs;
  FILE *outfs;

  time_class t1;
  time_class t2;

  meas_data *data;
  meas_data *sdata;
  long n;
  meas_data *selected;
  long n2;

  int32_t twid=TFIELD_WIDTH;
  int hemi=0;
  char c;
  long ncon;

  //latitude range:
  float min=-90;
  float max=90;

  void *optarg[10];
  int flag[10];
  int sflag;

  optarg[2]=&twid;
  optarg[7]=&min;
  optarg[8]=&max;
  argc=parse_command_opts(argc, argv, "ifd+-S?IF", "%s%s%d%%%%%g%g", 
		  optarg, flag, OPT_WHITESPACE);
  if (argc<0) {
    fprintf(stderr, "select_meas: command option parse error\n");
    return FATAL_COMMAND_OPTION_PARSE_ERROR;
  }
  if (flag[3]) hemi=1; else if (flag[4]) hemi=-1;
  if (flag[0]) t1.read_string((char *) optarg[0]);
  if (flag[1]) t2.read_string((char *) optarg[1]);
  sflag=flag[5];

  if (flag[6] || argc>3) {
    FILE *docfs;
    int err;
    if (flag[6]) {
      docfs=stdout;
      err=0;
    } else {
      docfs=stderr;
      err=INSUFFICIENT_COMMAND_ARGS;
    }
    fprintf(docfs, "Usage: select_meas [-i t0] [-f tf] [-+|--] [-d dwid] [-S] infile\n");
    fprintf(docfs, "options:\n");
    fprintf(docfs, "  -I   select measurements at this latitude or above\n");
    fprintf(docfs, "  -F   select measurements at this latitude or below\n"); 
    ctraj_optargs(docfs, "ifd+-S?");
    return err;
  }

  if (argc > 1) {
    infs=fopen(argv[1], "r");
    if (argc > 2) {
      outfs=fopen(argv[2], "w");
    } else {
      outfs=stdout;
    }
  } else {
    infs=stdin;
  }

  data=read_meas(infs, &n, twid);

  if (sflag) {
    sdata=sort_meas(data, n);
    delete [] data;
    data=sdata;
  }

  if (flag[0]==0) t1=data[0].t;
  if (flag[1]==0) t2=data[n-1].t;

  selected=select_meas(t1, t2, data, n, &n2, hemi);

  if (flag[7] || flag[8]) {
    meas_data *sel2;
    sel2=select_lat_range(selected, n2, min, max, &n2);
    delete [] selected;
    selected=sel2;
  }

  write_meas(selected, n2, stdout);

  delete [] data;
  delete [] selected;

}

