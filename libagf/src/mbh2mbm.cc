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

  multiclass_hier<real_a, cls_ta> *classifier;

  FILE *fs;

  int errcode=0;

/*
  int normflag=0;
  int svmflag=0;
  int hflag=0;
  int cflag=0;
*/

  //keep it stand-alone:
/*
  void *optarg[10];
  int flag[10];

  parse_command_opts(argc, argv, "nMhc", "%%%%", optarg, flag);
  normflag=flag[0];
  svmflag=flag[1];
  hflag=flag[2];
  cflag=flag[3];
*/

  //agf_parse_command_opts(argc, argv, "CHLMn", &opts);

//  long *clind;
//  long ncl;
//
  
  if (argc < 3) {
    printf("Converts multi-borders hierarchical (recursive control file) to\n");
    printf("multi-borders modes (single ASCII file).\n");
    printf("\n");
    printf("usage: mbh2mbm control outfile\n\n");
    printf("    where:\n");
    printf("control  = base name of binary input files.\n");
    printf("outfile  = output file in ASCII format.\n");
    exit(1);
  }

  classifier=new multiclass_hier<real_a, cls_ta>(argv[1]);
  fs=fopen(argv[2], "w");
  classifier->save(fs);

  delete classifier;

  return errcode;

}

