#include <stdio.h>
#include <stdlib.h>

#define MAXLL 1000

int main(int argc, char **argv) {
  FILE *fs=stdin;
  FILE *outfs=stdout;
  char line[MAXLL];
  float **data;		//data read in from file
  float **std;		//standard deviations
  float dum;
  int ncol;		//number of columns
  int offset;		//offset column
  int ndata;		//number of datasets (rows)
  int *maxind;		//max value

  if (argc<4) {
    fprintf(stderr, "write_results offset ncol dataset1 [ dataset2 [ dataset3 ...]]\n");
    return 1;
  }

  offset=atof(argv[1]);
  ncol=atof(argv[2]);

  ndata=argc-3;

  data=new float*[ndata];
  data[0]=new float[ndata*ncol];
  std=new float*[ndata];
  std[0]=new float[ndata*ncol];

  maxind=new int[ndata];

  for (int i=0; i<ndata; i++) {
    data[i]=data[0]+i*ncol;
    std[i]=std[0]+i*ncol;
    for (int j=0; j<offset; j++) fscanf(fs, "%g", &dum);
    maxind[i]=0;
    for (int j=0; j<ncol; j++) {
      fscanf(fs, "%g", data[i]+j);
      if (data[i][j]<data[i][maxind[i]]) maxind[i]=j;
    }
    fgets(line, MAXLL, fs);
    for (int j=0; j<offset; j++) fscanf(fs, "%g", &dum);
    for (int j=0; j<ncol; j++) fscanf(fs, "%g", std[i]+j);
    fgets(line, MAXLL, fs);
  }

  //fscanf(fs, "%10", name);
  for (int i=0; i<ndata; i++) {
    fprintf(outfs, " & %s", argv[i+3]);
  }
  fprintf(outfs, "\\\\\n");
  for (int i=0; i<ncol; i++) {
    for (int j=0; j<ndata; j++) {
      //if (maxind[j]==i) {
      //  fprintf(outfs, "& $\\mathbf{%g \\pm %g}$ ", data[j][i], std[j][i]);
      //} else {
        fprintf(outfs, "& $%.3g \\pm %.1g $ ", data[j][i], std[j][i]);
      //}
    }
    fprintf(outfs, "\\\\\n");
  }

  delete [] data[0];
  delete [] data;
  delete [] std[0];
  delete [] std;

  delete [] maxind;

  fclose(outfs);

}  

