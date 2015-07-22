#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include <stdint.h>

#include "error_codes.h"
#include "parse_command_opts.h"
#include "time_class.h"
#include "linked.h"
#include "contour_anal.h"

#include "ctraj_defaults.h"

using namespace libpetey;
using namespace ctraj;

int main(int argc, char **argv) {
  FILE *fs;

  int32_t n=1;
  int32_t npt;
  long npt2;
  int32_t magic;
  int16_t yr, mon, dy, hr, min;
  float sec;

  time_class t;

  float *x, *y;

  void *optarg[10];
  int optflag[10];

  argc=parse_command_opts(argc, argv, "?", "%", optarg, optflag, OPT_WHITESPACE);

  if (argc < 2 || optflag[0]) {
    FILE *docfs;
    int err;
    if (optflag[0]) {
      docfs=stdout;
      err=0;
    } else {
      docfs=stderr;
      err=INSUFFICIENT_COMMAND_ARGS;
    }
    fprintf(docfs, "Usage: bev2xy outfile [date]\n");
    fprintf(docfs, "\nConverts an ascii file containing point pairs to a single-record\n");
    fprintf(docfs, "binary file.  Date is optional since it is ignored for initial contours.\n");
    fprintf(docfs, "Takes input from stdin\n");
    fprintf(docfs, "\nWhere:\n\n");
    //printf("- infile is the input ascii file of point pairs\n");
    fprintf(docfs, "- outfile is the binary output file\n");
    fprintf(docfs, "- date is the date as: yyyy/mm/dd-hh:mm:ss\n");
    return err;
  }

  if (argc > 2) {
    t.read_string(argv[2]);
  } else {
    //arbitrary date:
    t.init(2000, 1, 1, 0, 0, 0);
  }

  //fs=fopen(argv[1], "r");
  fs=stdin;

  npt2=read_ascii_contour(fs, x, y);

  npt=npt2;

  //fclose(fs);

  fs=fopen(argv[1], "w");

  magic=MAGIC;

  fwrite(&magic, sizeof(magic), 1, fs);
  fwrite(&n, sizeof(n), 1, fs);

  t.get_fields(yr, mon, dy, hr, min, sec);

  fwrite(&yr, sizeof(yr), 1, fs);
  fwrite(&mon, sizeof(mon), 1, fs);
  fwrite(&dy, sizeof(dy), 1, fs);
  fwrite(&hr, sizeof(hr), 1, fs);
  fwrite(&min, sizeof(min), 1, fs);
  fwrite(&sec, sizeof(sec), 1, fs);
  fwrite(&npt, sizeof(npt), 1, fs);

  fwrite(x, sizeof(float), npt, fs);
  fwrite(y, sizeof(float), npt, fs);

  fclose(fs);

  delete [] x;
  delete [] y;

}


