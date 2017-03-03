#include <stdio.h>
#include <stdlib.h>
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
  int **ind;
  int nsum[argc];	//number of cols to sum
  int nout;		//number of output cols
  vector<float> row;
  double *mean;
  double *std;
  int err=0;

  ind=new int *[argc];
  ind[0]=new int[argc];

  nout=0;
  nsum[0]=0;
  for (int i=1; i<argc; i++) {
    if (argv[i][0]=='/') {
      ind[nout+1]=ind[nout]+nsum[nout];
      nout++;
      nsum[nout]=0;
      continue;
    }
    ind[nout][nsum[nout]]=atoi(argv[i]);
    nsum[nout]++;
  }
  nout++;

  fs=stdin;
  fgets(line, NCHAR, fs);
  while (sscanf(line+k, "%g%n", &val, &j)==1) {
    row.push_back(val);
    k+=j;
    //printf("%g ", val);
  }
  //printf("\n");
  n=row.size();

  for (int i=0; i<nout; i++) {
    val=0;
    for (int j=0; j<nsum[i]; j++) {
      val+=row[ind[i][j]];
    }
    printf("%g ", val);
  }
  printf("\n");

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
    for (int i=0; i<nout; i++) {
      val=0;
      for (int j=0; j<nsum[i]; j++) {
        val+=row[ind[i][j]];
      }
      printf("%g ", val);
    }
    printf("\n");
  }

}

