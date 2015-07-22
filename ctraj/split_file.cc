#include <stdio.h>
#include "randomize.h"

using namespace libpetey;

int main(int argc, char **argv) {
  char line[200];
  FILE *infs;
  FILE **outfs;
  FILE *docfs;

  int n, j;

  ran_init();

  if (argc < 3) {
    docfs=stdout;
    fprintf(docfs, "Usage: split_file infile outfile1 outfile2 ... outfileN\n");
    fprintf(docfs, "\n");
    fprintf(docfs, "Splits an ascii file into N different files\n");
    fprintf(docfs, "by randomly selecting lines\n");
    exit(0);
  }

  n=argc-2;

  infs=fopen(argv[1], "r");

  outfs=new FILE*[n];
  for (int i=0; i<n; i++) {
    outfs[i]=fopen(argv[i+2], "w");
  }

  while (feof(infs)==0) {
    j=ranu()*n;
    fgets(line, 200, infs);
    fprintf(outfs[j], "%s", line);
  }

  fclose(infs);
  for (int i=0; i<n; i++) fclose(outfs[i]);
  
  delete [] outfs;

}  

