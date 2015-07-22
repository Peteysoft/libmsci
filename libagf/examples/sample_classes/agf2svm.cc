#include <stdio.h>
#include <string.h>
#include <assert.h>

//#include "/home/Petey/software_projects/include/class_lib.h"

int main(int argc, char ** argv) {

  char *vecfile;
  char *clsfile;

  FILE *fs;

  long n;
  long nvar;

  float **vec;
  long *cls;

  float *all;

  float *ave, *std;
  int normflag;

//  long *clind;
//  long ncl;
//
  
  if (argc < 2) {
    printf("Converts vector data in libAGF compatible binary format\n");
    printf("to LIBSVM compatible ASCII format.\n\n");
    printf("usage: agf2svm [-n] basename\n\n");
    printf("    where:\n");
    printf("-n      = option to normalise the data\n");
    printf("basename = base name of binary files:\n");
    printf("           .vec for vector data, .cls for class data\n");
    exit(-1);
  }

  //generate the input file names:
  vecfile=new char[strlen(argv[1])+5];
  strcpy(vecfile, argv[1]);
  strcat(vecfile, ".vec");

  clsfile=new char[strlen(argv[1])+5];
  strcpy(clsfile, argv[1]);
  strcat(clsfile, ".cls");

  //read in the training data:
  fs=fopen(vecfile, "r");
  fread(&nvar, sizeof(nvar), 1, fs);
  fseek(fs, 0, SEEK_END);
  n=(ftell(fs)-4)/nvar/4;

  //printf("%d training vectors found in file %s\n", n, vecfile);

  fseek(fs, 4, SEEK_SET);
  vec=new float *[n];
  all=new float [n*nvar];
  fread(all, sizeof(float), n*nvar, fs);
  for (long i=0; i<n; i++) vec[i]=all+i*nvar;

  fclose(fs);

  //keep these as stand-alone utilities:
  //
  if (normflag) {
    float diff;

    ave=new float[nvar];
    std=new float[nvar];

    for (long i=0; i<nvar; i++) {
      ave[i]=0;
      std[i]=0;
      for (long j=0; j<n; j++) ave[i]+=vec[j][i];
      ave[i]/=n;
      for (long j=0; j<n; j++) {
        diff=vec[j][i]-ave[i];
        vec[j][i]=diff;
        std[i]+=diff*diff;
      }
      std[i]=sqrt(std[i]/(n-1));
      for (long j=0; j<n; j++) vec[j][i]/=std[i];
    }
  }

  fs=fopen(clsfile, "r");

  if (fs != NULL) {
    fseek(fs, 0, SEEK_END);
    assert(n == ftell(fs)/4);

    rewind(fs);
    cls=new long[n];
    fread(cls, sizeof(long), n, fs);
  } else {
    cls=new long [n];
    for (long i=0; i<n; i++) cls[i]=0;
  }

  fclose(fs);

//  clind=sort_classes(vec, nvar, cls, n, ncl);

  for (long i=0; i<n; i++) {
    printf("%d ", cls[i]);
    for (long j=0; j<nvar; j++) printf("%d:%g ", j+1, vec[i][j]);
    printf("\n");
  }

  delete [] all;
  delete [] vec;
  delete [] cls;

}

