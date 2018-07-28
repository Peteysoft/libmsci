#include <stdio.h>
#include <stdlib.h>

#include "tree_lg.h"

using namespace std;
using namespace libpetey;

//simple executable that takes a list of numbers and spits out the
//k least

int main(int argc, char **argv) {
  tree_lg<float> sorter;
  float data;
  long k;
  float *result;
  FILE *fs;

  if (argc < 2) {
    printf("Usage: kleast k [file]\n");
    exit(1);
  }

  sscanf(argv[1], "%d", &k);

  if (argc > 2) fs=fopen(argv[2], "r"); else fs=stdin;

  result=new float[k];

  for (long i=0; i<k; i++) {
    fscanf(fs, "%g", &data);
    sorter.add(data);
  }

  while(feof(fs) == 0) {
    fscanf(fs, "%g", &data);
    sorter.add(data);
    sorter.delete_greatest();
  }

  sorter.decompose(result, k);

/*  for (long i=0; i<k; i++) {
    printf("%g\n", result[i]);
  }*/

  delete [] result;

}


