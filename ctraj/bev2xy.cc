#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>

#include <stdint.h>

#include "parse_command_opts.h"
#include "error_codes.h"
#include "time_class.h"

#include "ctraj_defaults.h"

using namespace libpetey;
using namespace ctraj;

//#include "boundary_list_g.h"

int main(int argc, char **argv) {
  FILE *fs;

  int32_t n;
  int32_t npt;
  int32_t magic;
  int32_t ind;
  int16_t yr, mon, dy, hr, min;
  float sec;

  time_class t;
  char tstring[30];

  float *x, *y;
  int32_t i;
  char response[3];

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
    fprintf(docfs, "Usage: bev2xy file [index]\n");
    fprintf(docfs, "\nConverts a single record in a binary file containing contour\n");
    fprintf(docfs, "advection results to ascii: sends results to standard out\n");
    fprintf(docfs, "\nWhere:\n\n");
    fprintf(docfs, "- file is the input file\n");
    fprintf(docfs, "- index is the desired record\n");
    fprintf(docfs, "  (to list the contents of the file, omit index)\n");
    return err;
  }

  fs=fopen(argv[1], "r");
  if (fs==NULL) {
    fprintf(stderr, "bev2xy: cannot open file, %s, for reading\n", argv[1]);
    return UNABLE_TO_OPEN_FILE_FOR_READING;
  }
  if (argc > 2) sscanf(argv[2], "%d", &ind);

  fread(&magic, sizeof(magic), 1, fs);
  if (magic != MAGIC) {
    fprintf(stderr, "bev2xy: file, %s, is wrong type\n", argv[1]);
    fclose(fs);
    return FILE_READ_ERROR;
  }
  fread(&n, sizeof(n), 1, fs);

  if (argc > 2 && (ind < 0 || ind >= n)) {
    fprintf(stderr, "bev2xy: index, %d, out of range [0, %d]\n", ind, n-1);
    return PARAMETER_OUT_OF_RANGE;
  }

  if (argc == 2) ind=n;

  for (i=0; i<ind; i++) {
    if (feof(fs)!=0) break;
    fread(&yr, sizeof(yr), 1, fs);
    if (feof(fs)!=0) goto incomplete_record;
    fread(&mon, sizeof(mon), 1, fs);
    if (feof(fs)!=0) goto incomplete_record;
    fread(&dy, sizeof(dy), 1, fs);
    if (feof(fs)!=0) goto incomplete_record;
    fread(&hr, sizeof(hr), 1, fs);
    if (feof(fs)!=0) goto incomplete_record;
    fread(&min, sizeof(min), 1, fs);
    if (feof(fs)!=0) goto incomplete_record;
    fread(&sec, sizeof(sec), 1, fs);
    if (feof(fs)!=0) goto incomplete_record;
    fread(&npt, sizeof(npt), 1, fs);
    if (feof(fs)!=0) goto incomplete_record;
    if (argc == 2) {
      t.init(yr, mon, dy, hr, min, sec);
      t.write_string(tstring);
      printf("%d %s %d\n", i, tstring, npt);
    }
    if (fseek(fs, 8*npt, SEEK_CUR)!=0) goto incomplete_record;
  }

  if (argc==2) {
    if (n!=i) goto fix_header;
    return 0;
  }
  if (feof(fs)!=0) goto fix_header;

  fread(&yr, sizeof(yr), 1, fs);
  fread(&mon, sizeof(mon), 1, fs);
  fread(&dy, sizeof(dy), 1, fs);
  fread(&hr, sizeof(hr), 1, fs);
  fread(&min, sizeof(min), 1, fs);
  fread(&sec, sizeof(sec), 1, fs);
  fread(&npt, sizeof(npt), 1, fs);

  x=new float[npt];
  y=new float[npt];

  fread(x, sizeof(float), npt, fs);
  fread(y, sizeof(float), npt, fs);

  fclose(fs);

  for (i=0; i<npt; i++) {
    printf("%g %g\n", x[i], y[i]);
  }

  delete [] x;
  delete [] y;

  return 0;

  incomplete_record:
    fprintf(stderr, "bev2xy: record %d is incomplete\n", i);

  fix_header:
    fclose(fs);
    fprintf(stderr, "bev2xy: header data is incorrect (%d vs. %d)\n", i, n);
    fprintf(stderr, "        fix now (y/n)?\n");
    scanf("%1s", response);
    if (strcmp(response, "y")==0) {
      fs=fopen(argv[1], "r+");
      fseek(fs, sizeof(magic), SEEK_SET);
      n=i;
      fwrite(&n, sizeof(n), 1, fs);
      fclose(fs);
    }
  return INTERNAL_ERROR;

}


