#include <stdio.h>
#include <stdlib.h>

#define MAXLL 1000

int main(int argc, char **argv) {
  FILE *fs=stdin;
  FILE *outfs=stdout;
  int nfield=100;		//number of columns
  int ncol=3;		//total number of fields
  char name[11];
  char line[MAXLL];
  float data[nfield];	//data read in from file
  float std[nfield];	//standard deviations
  int ind;		//index of extremum
  float ext;		//extremum

  nfield=atof(argv[1]);
  ncol=atof(argv[2]);

  //fscanf(fs, "%10", name);
  for (int i=0; i<nfield; i++) fscanf(fs, "%g", data+i);
  //fgets(line, MAXLL, fs);
  for (int i=0; i<nfield; i++) fscanf(fs, "%g", std+i);
  //fgets(line, MAXLL, fs);
  //print out name of dataset:
  if (argc>=4) fprintf(outfs, "%s", argv[3]);

  for (int i=0; i<nfield; i++) {
    float val=data[i];
    float sd=std[i];
    fprintf(outfs, " & $");
    if (ind==i) fprintf(outfs, "\\mathbf{");
    fprintf(outfs, "%12.6g", val);
    if (sd>0) fprintf(outfs, "\\pm%9.2g", sd);
    if (ind==i) fprintf(outfs, "}");
    fprintf(outfs, "$");
    if ((i+1) % ncol == 0) {
      fprintf(outfs, "\\\\\n");
      //if (argc>=4) fprintf(outfs, " & ");
    }
  }
  //fprintf(outfs, "\\\\\n");

  fclose(outfs);

}  

