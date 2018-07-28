#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "error_codes.h"
#include "parse_command_opts.h"
#include "agf_defs.h"

using namespace std;
using namespace libagf;
using namespace libpetey;

//dirt simple application to convert floating point data to a pair
//of classes based on a threshold:

int main(int argc, char ** argv) {
  FILE *fs;
  nel_ta n;			//number of ordinates
  real_a * data;		//floating point data
  cls_ta * cls;			//transformed discrete data

  void *optarg[20];
  int flag[20];

  argc=parse_command_opts(argc, argv, "", "", optarg, flag, 1);
  if (argc<0) {
    fprintf(stderr, "float_to_class: command option parse error\n");
    exit(FATAL_COMMAND_OPTION_PARSE_ERROR);
  }

  if (argc < 3) {
    printf("Usage:  class_to_float class_file float_file\n");
    printf("	where:\n");
    printf("class_file	input file containing class data\n");
    printf("float_file	output file containing floating-point data\n");
    exit(1);
  }

  fs=fopen(argv[1], "r");
  fseek(fs, 0, SEEK_END);
  n=ftell(fs)/sizeof(real_a);
  fseek(fs, 0, SEEK_SET);
  data=new real_a[n];
  cls=new cls_ta[n];
  fread(cls, n, sizeof(cls_ta), fs);
  fclose(fs);

  for (nel_ta i=0; i<n; i++) data[i]=cls[i];

  fs=fopen(argv[2], "w");
  fwrite(data, n, sizeof(real_a), fs);
  fclose(fs);

  delete [] data;
  delete [] cls;

  return 0;

}
    

