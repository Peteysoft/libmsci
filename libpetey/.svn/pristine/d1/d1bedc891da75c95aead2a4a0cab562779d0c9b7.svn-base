#include <stdio.h>
#include <stdlib.h>

#include "quicksort.h"
#include "peteys_tmpl_lib.h"
#include "parse_command_opts.h"
#include "tree_lg.h"
#include "randomize.h"
#include "error_codes.h"
#include "read_ascii_all.h"
#include "kselect.h"

using namespace std;
using namespace libpetey;

//part of the test suite

//simple executable that takes a list of numbers and spits out 
//a sorted list or the k least 

int main(int argc, char **argv) {
  //tree_lg<float> sorter;
  float *data;
  long n;
  int32_t k;
  FILE *fs;
  char **line;

  //parse the arguments:
  void *opt[20];
  int flag[20];
  opt[0]=&k;

  argc=parse_command_opts(argc, argv, "k?qhtriCs", "%d%%%%%%%%", opt, flag, 1);
  if (argc<0) {
    fprintf(stderr, "sorter: error parsing command line\n");
    //no second chances:
    exit(FATAL_COMMAND_OPTION_PARSE_ERROR);
  }

  if (flag[1]) {
    printf("Usage: sorter [-k k] [-q] [-h] [-t] [file]\n");
    printf("  where:\n");
    printf("k      select k nearest neighbours\n");
    printf("file   optional input file (or take data from stdin)\n");
    printf("-q     quicksort\n");
    printf("-h     heapsort\n");
    printf("-t     treesort\n");
    printf("-s     selection sort\n");
    printf("-r     randomize the list instead of sorting it\n");
    printf("-i     output indices\n");
    printf("-C     check if input is sorted\n");
    exit(0);
  }

  if (argc > 1) fs=fopen(argv[1], "r"); else fs=stdin;

  line=read_ascii_all(fs, &n);
  data=new float[n];
  for (long i=0; i<n; i++) {
    if (sscanf(line[i], "%g", data+i)!=1) {
      n=i;
      break;
    }
    delete [] line[i];
  }
  fclose(fs);

  delete [] line;

  if (flag[7]) {
    for (long i=1; i<n; i++) {
      if (data[i]<data[i-1]) {
        fprintf(stderr, "sorter: input list is not sorted\n");
        exit(SYNTAX_ERROR);
      }
    }
    delete [] data;
    return 0;
  } else if (flag[6]) {  
    long *result;
    long nres;
    if (flag[0]) { 
      kiselect_base<float> *selector;
      float *dum;

      if (flag[2]) {
        selector=new kiselect_quick<float>(k);
      } else if (flag[3]) {
        selector=new kiselect_heap<float>(k);
      } else if (flag[4]) {
        selector=new kiselect_tree<float>(k);
      } else if (flag[8]) {
        selector=new kiselect_naive<float>(k);
      } else {
        selector=new kiselect_quick<float>(k);
      }

      for (long i=0; i<n; i++) selector->add(data[i]);
      dum=new float[k];
      result=new long[k];
      selector->get(dum, result);
      nres=k;

      delete selector;
      delete [] dum;
    } else if (flag[5]) {
      ran_init();
      result=randomize(n);
      nres=n;
      ran_end();
    } else {

      if (flag[2]) {
        result=new long[n];
        quicksort(data, result, n);
      } else if (flag[3]) {
        result=heapsort(data, n);
      } else if (flag[4]) {
        result=treesort(data, n);
      } else {
        result=new long[n];
        quicksort(data, result, n);
      }
      nres=n;
    }

    for (long i=0; i<nres; i++) printf("%d\n", result[i]);

    delete [] result;
  } else {
    float *result;
    long nres;
    if (flag[0]) { 
      kselect_base<float> *selector;

      if (flag[2]) {
        selector=new kselect_quick<float>(k);
      } else if (flag[3]) {
        selector=new kselect_heap<float>(k);
      } else if (flag[4]) {
        selector=new kselect_tree<float>(k);
      } else {
        selector=new kselect_quick<float>(k);
      }

      for (long i=0; i<n; i++) selector->add(data[i]);
      result=new float[k];
      selector->get(result);
      nres=k;

      delete selector;
    } else if (flag[5]) {
      ran_init();
      long *ind=randomize(n);
      result=map_vector(data, ind, n);
      delete [] ind;
      nres=n;
      ran_end();
    } else {

      if (flag[2]) {
        quicksort(data, n);
      } else if (flag[3]) {
        heapsort_inplace(data, n);
      } else if (flag[4]) {
        treesort(data, n);
      } else {
        quicksort(data, n);
      }
      result=data;
      nres=n;
    }

    for (long i=0; i<nres; i++) printf("%g\n", result[i]);

    if (flag[0] || flag[5]) delete [] result;
  }

  delete [] data;

}


