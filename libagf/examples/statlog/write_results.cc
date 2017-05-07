#include <stdio.h>

#define MAXLL 1000

int main(int argc, char **argv) {
  FILE *fs=stdin;
  FILE *outfs=stdout;
  int ncol=6;		//total number of columns
  char line[MAXLL];
  float data[ncol];	//data read in from file
  float std[ncol];	//standard deviations

  //fscanf(fs, "%10", name);
  for (int i=0; i<ncol; i++) fscanf(fs, "%g", data+i);
  fgets(line, MAXLL, fs);
  for (int i=0; i<ncol; i++) fscanf(fs, "%g", std+i);
  fgets(line, MAXLL, fs);
  //print out name of dataset:
  if (argc>=2) fprintf(outfs, "%s ", argv[1]);

  for (int i=0; i<ncol; i++) {
    fprintf(outfs, " & $%12.3g$", data[i]);
    if (std[i]>0) fprintf(outfs, "\\pm%9.2g$", std[i]);
  }
  fprintf(outfs, "\\\\\n");

  fclose(outfs);

  return 0;

}  

