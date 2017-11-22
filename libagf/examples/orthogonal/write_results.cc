#include <stdio.h>

#define MAXLL 1000

int main(int argc, char **argv) {
  FILE *fs=stdin;
  FILE *outfs=stdout;
  int ncol=10;		//total number of columns
  int namelen=10;
  char name[11];
  char line[MAXLL];
  float data[ncol];	//data read in from file
  float std[ncol];	//standard deviations
  int ind;		//index of extremum
  float ext;		//extremum

  //fscanf(fs, "%10", name);
  for (int i=0; i<ncol; i++) fscanf(fs, "%g", data+i);
  //fgets(line, MAXLL, fs);
  for (int i=0; i<ncol; i++) fscanf(fs, "%g", std+i);
  //fgets(line, MAXLL, fs);
  //print out name of dataset:
  if (argc>=2) fprintf(outfs, "%s", argv[1]);

  for (int i=0; i<ncol; i++) {
    float val=data[i];
    float sd=std[i];
    fprintf(outfs, " & $");
    if (ind==i) fprintf(outfs, "\\mathbf{");
    fprintf(outfs, "%12.6g", val);
    if (sd>0) fprintf(outfs, "\\pm%9.2g", sd);
    if (ind==i) fprintf(outfs, "}");
    fprintf(outfs, "$");
  }
  fprintf(outfs, "\\\\\n");

  fclose(outfs);

}  

