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

void find_minmax(real_a *data, nel_ta n, real_a &min, real_a &max) {
  min=data[0];
  max=data[0];
  for (nel_ta i=1; i<n; i++) {
    if (data[i]>max) max=data[i];
    else if (data[i]<min) min=data[i];
  }
}

void find_minmax(real_a *data, nel_ta n, real_a &min, real_a &max, real_a tmiss) {
  min=data[0];
  max=data[0];
  for (nel_ta i=1; i<n; i++) {
    if (data[i]>max && data[i] < tmiss) max=data[i];
    else if (data[i]<min) min=data[i];
  }
}


int main(int argc, char ** argv) {
  FILE *fs;
  nel_ta n;			//number of ordinates
  real_a * data;		//floating point data
  real_a *thresh;		//list of threshold values ("bins")
  cls_ta maxcls;		//maximum value for class label
  cls_ta * cls;			//transformed discrete data

  void *optarg[20];
  int flag[20];
  int Mflag;			//find min/max
  int gflag;			//geometric progression
  int Qflag;			//print class break-down
  real_a min, max;
  real_a min1, max1;
  real_a tmiss;			//"out-of-bounds"
  int Tflag;

  cls_ta *bin;			//histogram of classes

  int nread;

  optarg[0]=&maxcls;
  optarg[3]=&tmiss;
  argc=parse_command_opts(argc, argv, "qgMTQ", "%d%%%%", optarg, flag, 5);
  if (argc<0) {
    fprintf(stderr, "float_to_class: command option parse error\n");
    exit(FATAL_COMMAND_OPTION_PARSE_ERROR);
  }
  //-q option specifies number of classes, not divisions:
  if (flag[0]) maxcls--;

  gflag=flag[1];
  Mflag=flag[2];
  Tflag=flag[3];
  Qflag=flag[4];

  if (argc < 2) {
    printf("Usage:  float_to_class [-M] [-q nz] [-g] float_file [class_file \n");
    printf("              {[min max] | thresh1 [thresh2 [ thresh3... ]]]}\n");
    printf("	where:\n");
    printf("float_file	file containing floating-point data\n");
    printf("class_file	output file containing class data\n");
    printf("threshN	threshold value(s) defining classes\n");
    printf("min		minimum value in divison scheme\n");
    printf("max		maximum value in division scheme\n");
    printf("-M          print out min and max values\n");
    printf("-q		number of classes/(divisions-1) to use\n");
    printf("-g		use geometric progression\n");
    printf("-T		threshold for missing data or outliers\n");
    printf("-Q		print class break-down; otherwise, print bin thresholds\n");
    exit(1);
  }

  fs=fopen(argv[1], "r");
  fseek(fs, 0, SEEK_END);
  n=ftell(fs)/sizeof(real_a);
  fseek(fs, 0, SEEK_SET);
  data=new real_a[n];
  cls=new cls_ta[n];
  fread(data, n, sizeof(real_a), fs);
  fclose(fs);

  if (Mflag || (flag[0] && argc<5)) {
    if (Tflag) find_minmax(data, n, min, max, tmiss);
    		else find_minmax(data, n, min, max);
  }

  if (argc == 2) {
    if (Mflag) {
      printf("%g\n%g\n", min, max);
      return 0;
    } else {
      for (nel_ta i=0; i<n; i++) printf("%g\n", data[i]);
      return 0;
    }
  }

  if (Mflag) {
    min1=min;
    max1=max;
  }
  if (flag[0] && argc>3) {
    nread=sscanf(argv[3], "%g", &min);
    if (nread != 1) fprintf(stderr, "float_to_class: error parsing third command line argument.\n");
    if (argc>4) {
      sscanf(argv[4], "%g", &max);
      if (nread != 1) {
        real_a dum;
        fprintf(stderr, "float_to_class: error parsing fourth command line argument.\n");
        if (Tflag) find_minmax(data, n, min, max, tmiss);
    		else find_minmax(data, n, min, max);
      }
    }
  }

  if (flag[0]) {
    real_a dthresh;
    thresh=new real_a[maxcls];
    if (gflag) {
      dthresh=pow(max/min, 1./(maxcls-1));
      thresh[0]=min;
      for (cls_ta i=1; i<maxcls; i++) thresh[i]=thresh[i-1]*dthresh;
    } else {
      dthresh=(max-min)/(maxcls-1);
      for (cls_ta i=0; i<maxcls; i++) thresh[i]=min+dthresh*i;
    }
  } else {
    maxcls=argc-3;
    thresh=new real_a[maxcls];
    for (int i=3; i<argc; i++) {
      nread=sscanf(argv[i], "%g", thresh+i-3);
      if (nread!=1) {
        fprintf(stderr, "float_to_class: error parsing %dth command line argument", i);
        exit(FATAL_COMMAND_OPTION_PARSE_ERROR);
      }
    }
  }

  if (Mflag) printf("%g\n%g\n", min1, max1);

  bin=new nel_ta[maxcls+1];
  for (cls_ta i=0; i<=maxcls; i++) bin[i]=0;
  for (nel_ta i=0; i<n; i++) {
    if (data[i]<=thresh[0]) {
      cls[i]=0;
      bin[0]++;
    } else if (data[i]>thresh[maxcls-1]) {
      cls[i]=maxcls;
      bin[maxcls]++;
    } else {
      for (cls_ta j=1; j<maxcls; j++) {
        if (data[i]>thresh[j-1] && data[i]<=thresh[j]) {
          cls[i]=j;
          bin[j]++;
          break;
        }
      }
    }
    //printf("%g: %d\n", data[i], cls[i]);
  }

  if (Qflag) {
    printf("Class break-down:\n");
    for (cls_ta i=0; i<=maxcls; i++) {
      printf("%2d:  %d\n", i, bin[i]);
    }
  } else {
    for (cls_ta i=0; i<maxcls; i++) printf("%g\n", thresh[i]);
  }
    

  fs=fopen(argv[2], "w");
  fwrite(cls, n, sizeof(cls_ta), fs);
  fclose(fs);

  delete [] data;
  delete [] cls;

  return 0;

}
    

