#include "../simple_temp.h"
#include "../dependent_temp.h"
#include <stdio.h>

int main(int argc, char *argv[]) {
  simple<float> simple1;
  simple<float> simple2;
  dependent<float> *dep;
  FILE *infile;
  long nrow, ncol;
  long nlines;
  char line[200];
  long size;
  float value1, value2, value3;

  if (argc != 2) return 1;
  infile=fopen(argv[1], "r");

  fgets(line, size, infile);
  sscanf(line, "%d", &nlines);

  dep=new dependent<float>(&simple1, &simple2);

  for (long i=0; i<nlines; i++) {
    fgets(line, size, infile);
    sscanf(line, "%f%f%f", &value1, &value2, &value3);
    dep->cel(value3, simple1.add_el(value1), simple2.add_el(value2));

    ncol=simple1.nel();
    nrow=simple2.nel();
    
    printf("          ");
    for (long i=0; i<ncol; i++) {
      simple1.get(value1, i);
      printf("%10.2f", value1);
    }
    printf("\n");

    for (long i=0; i<nrow; i++) {
      simple2.get(value1, i);
      printf("%10.2f", value1);
      for (long j=0; j<ncol; j++) {
        dep->get(value1, j, i);
        printf("%10.2f", value1);
      }
      printf("\n");
    }
    printf("\n");
  }

  ncol=simple1.nel();
  nrow=simple2.nel();

  printf("          ");
  for (long i=0; i<ncol; i++) {
    simple1.get(value1, i);
    printf("%10.2f", value1);
  }
  printf("\n");

  for (long i=0; i<nrow; i++) {
    simple2.get(value1, i);
    printf("%10.2f", value1);
    for (long j=0; j<ncol; j++) {
      dep->get(value1, j, i);
      printf("%10.2f", value1);
    }
    printf("\n");
  }
  printf("\n");

  delete dep;
  fclose(infile);

}



