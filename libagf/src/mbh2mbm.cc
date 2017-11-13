#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#include "getopt.h"

#include "full_util.h"
#include "agf_lib.h"

using namespace std;
using namespace libagf;
using namespace libpetey;

//#include "parse_command_opts.h"

//#include "/home/Petey/software_projects/include/class_lib.h"

int main(int argc, char ** argv) {
  agf_command_opts opts;
  multiclass_hier<real_a, cls_ta> *classifier;

  FILE *fs;

  int errcode=0;

  errcode=agf_parse_command_opts(argc, argv, "y:Z", &opts);

  if (argc < 2) {
    printf("Converts multi-borders hierarchical (recursive control file) to\n");
    printf("multi-borders modes (single ASCII file) using one of three modes:\n");
    printf("one-vs-one, one-vs-the-rest, partitioning of adjacent classes.\n");
    printf("\n");
    printf("usage: mbh2mbm control outfile\n\n");
    printf("    where:\n");
    printf("control  = name of control file.\n");
    printf("outfile  = output file in ASCII format.\n");
    exit(1);
  }

  classifier=new multiclass_hier<real_a, cls_ta>(argv[0], 0, opts.path, NULL, 0, 0, opts.Zflag);
  fs=fopen(argv[1], "w");
  classifier->save(fs);

  delete classifier;

  return errcode;

}

