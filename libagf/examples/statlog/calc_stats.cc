#include <stdio.h>
#include <math.h>
#include <vector>

#define NCHAR 10000

using namespace std;

int main(int argc, char **argv) {
  FILE *fs;
  char line[NCHAR];
  float val;
  int n, nline;
  int k=0;
  int j;
  vector<float> row;
  double *mean;
  double *std;
  int err=0;

  fs=fopen(argv[1], "r");
  fgets(line, NCHAR, fs);
  while (sscanf(line+k, "%g%n", &val, &j)==1) {
    row.push_back(val);
    k+=j;
    //printf("%g ", val);
  }
  //printf("\n");

  n=row.size();
  mean=new double[n];
  std=new double[n];

  for (int i=0; i<n; i++) {
    mean[i]=row[i];
    std[i]=row[i]*row[i];
  }
  nline=1;

  while (fgets(line, NCHAR, fs)!=NULL) {
    k=0;
    for (int i=0; i<n; i++) {
      if (sscanf(line+k, "%g%n", &val, &j)!=1) {
        err=1;
        break;
      }
      row[i]=val;
      k+=j;
      //printf("%g ", val);
    }
    //printf("\n");
    if (err) break;
    nline++;
    for (int i=0; i<n; i++) {
      mean[i]+=row[i];
      std[i]+=row[i]*row[i];
    }
  }

  for (int i=0; i<n; i++) {
    mean[i]/=nline;
    printf("%lg ", mean[i]);
  }
  printf("\n");

  for (int i=0; i<n; i++) {
    std[i]=sqrt((std[i]-nline*mean[i]*mean[i])/(nline-1));
    printf("%lg ", std[i]);
  }
  printf("\n");

}

