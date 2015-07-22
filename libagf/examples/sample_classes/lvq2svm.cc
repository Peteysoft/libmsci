#include <stdio.h>

int main(int argc, char **argv) {
  FILE *ifs;
  FILE *ofs;

  long nvar;
  float *vec;
  long cls;

  ifs=fopen(argv[1], "r");
  fscanf(ifs, "%d", &nvar);

  vec=new float [nvar];

  ofs=fopen(argv[2], "w");

  while (!feof(ifs)) {
    for (long i=0; i<nvar; i++) fscanf(ifs, "%g", vec+i);
    fscanf(ifs, "%d", &cls);
    fprintf(ofs, "%d ", cls);
    for (long i=0; i<nvar; i++) fprintf(ofs, "%d:%f ", i+1, vec[i]);
    fprintf(ofs, "\n");
  }

  fclose(ifs);
  fclose(ofs);

  delete[] vec;

}

