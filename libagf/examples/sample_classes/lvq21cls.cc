#include <stdio.h>
#include <string.h>

int main(int argc, char **argv) {
  char *vecfile;
  char *clsfile;

  char line [200];
  long lptr1, lptr;

  FILE *ifs;
  FILE *ofs;

  long nvar;
  float *vec;
  long cls;

  ifs=fopen(argv[1], "r");
  fgets(line, 200, ifs);
  sscanf(line, "%d", &nvar);

  ofs=fopen(argv[2], "w");

  fprintf(ofs, "%d\n", nvar);

  vec=new float [nvar];

  while (!feof(ifs)) {
    fgets(line, 200, ifs);
    lptr=0;
    for (long i=0; i<nvar; i++) {
      sscanf(line+lptr, "%g%n", vec+i, &lptr1);
      lptr+=lptr1;
    }
    sscanf(line+lptr, "%d", &cls);
    for (long i=0; i<nvar; i++) fprintf(ofs, "%g ", vec[i]);
    cls=1;
    fprintf(ofs, "%d\n", cls);
  }

  fclose(ifs);
  fclose(ofs);

  delete[] vec;

}

