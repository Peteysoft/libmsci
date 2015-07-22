#include <stdio.h>
#include <string.h>

#include "../class_lib.h"

#define MAXLL 200

int main(int argc, char **argv) {

  FILE *fs;
  char line[MAXLL];
  char *err;
  long nvar;
  long xdim, ydim;
  long n;
  char topl[10];
  float **cod;
  float *vec;
  long cls;
  long strptr1, strptr2;
  long *nhit;		//number of hits for a given code vector
  long *cp;
  float d2, mind2;
  long minind;

  //read in the code vectors:
  fs=fopen(argv[1], "r");

  fgets(line, MAXLL, fs);
  sscanf(line, "%d %s %d %d", &nvar, topl, &xdim, &ydim);
  n=xdim*ydim;

  cod=new float*[n];
  for (long i=0; i<n; i++) {
    cod[i]=new float[nvar];
    fgets(line, MAXLL, fs);
//    printf("%s", line);
    strptr1=0;
    for (long j=0; j<nvar; j++) {
      sscanf(line+strptr1, "%g%n", cod[i]+j, &strptr2);
//      printf("%g ", cod[i][j]);
      strptr1+=strptr2;
    }
//    printf("\n");
  }

  fclose(fs);


  nhit=new long[n];
  cp=new long[n];

  for (long i=0; i<n; i++) {
    nhit[i]=0;
    cp[i]=0;
  }

  //read in the training data one by one:
  fs=fopen(argv[2], "r");

  vec=new float[nvar];
  while (!feof(fs)) {
    err=fgets(line, MAXLL, fs);
    if (err==NULL || strlen(line) == 0) break;
    strptr1=0;
    for (long j=0; j<nvar; j++) {
      sscanf(line+strptr1, "%g%n", vec+j, &strptr2);
      strptr1+=strptr2;
    }
    sscanf(line+strptr1, "%d", &cls);

    //search for nearest code vector:
    minind=0;
    mind2=metric2(cod[0], vec, nvar);
    for (long i=1; i<n; i++) {
      d2=metric2(cod[i], vec, nvar);
      if (d2 < mind2) {
        mind2=d2;
        minind=i;
      }
    }
    nhit[minind]++;
    if (cls == 0) cp[minind]--; else cp[minind]++;
  }

  fclose(fs);

  printf("%d", nvar);
  for (long i=0; i<n; i++) {
    for (long j=0; j<nvar; j++) printf("%g ", cod[i][j]);
    if (nhit[i]==0) printf("%g\n", nhit[i]); 
		else printf("%g\n", (float) cp[i]/(float) nhit[i]);
  }

  delete [] nhit;
  delete [] cp;
  delete [] vec;
  for (long i=0; i<n; i++) delete [] cod[i];
  delete [] cod;

}

