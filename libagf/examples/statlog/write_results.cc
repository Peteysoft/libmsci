#include <stdio.h>

#define MAXLL 1000

int main(int argc, char **argv) {
  FILE *fs=stdin;
  FILE *outfs=stdout;
  int ncol=15;		//total number of columns
  int namelen=10;
  char name[11];
  char line[MAXLL];
  float data[ncol];	//data read in from file
  float std[ncol];	//standard deviations
  int ind;		//index of extremum
  float ext;		//extremum

  //fscanf(fs, "%10", name);
  for (int i=0; i<ncol; i++) fscanf(fs, "%g", data+i);
  fgets(line, MAXLL, fs);
  for (int i=0; i<ncol; i++) fscanf(fs, "%g", std+i);
  fgets(line, MAXLL, fs);
  //print out name of dataset:
  if (argc>=2) fprintf(outfs, "%s & train (s)", argv[1]); else fprintf(outfs, " & train (s)");

  //print out training times:
  ind=0;
  ext=data[0];
  for (int i=1; i<3; i++) {
    if (data[i*2-1]< ext) {
      ext=data[i*2-1];
      ind=i;
    }
  }
  if (ind==0) fprintf(outfs, " & {\\bf - }"); else fprintf(outfs, " & -");
  for (int i=1; i<4; i++) {
    float val=data[i*2-1];
    float sd=std[i*2-1];
    fprintf(outfs, " & $");
    if (ind==i) fprintf(outfs, "\\mathbf{");
    fprintf(outfs, "%12.3g", val);
    if (sd>0) fprintf(outfs, "\\pm%9.2g", sd);
    if (ind==i) fprintf(outfs, "}");
    fprintf(outfs, "$");
  }
  fprintf(outfs, "\\\\\n");

  //print out test times:
  ind=0;
  ext=data[0];			//KNN test is same as trainging time
  for (int i=1; i<4; i++) {
    if (data[i*2] < ext) {
      ext=data[i*2];
      ind=i;
    }
  }
  fprintf(stderr, "%g\n", std[0]);
  for (int i=0; i<4; i++) {
    float val=data[i*2];
    float sd=std[i*2];
    fprintf(outfs, " & $");
    if (ind==i) fprintf(outfs, "\\mathbf{");
    fprintf(outfs, "%12.3g", val);
    if (sd>0) fprintf(outfs, "\\pm%9.2g", sd);
    if (ind==i) fprintf(outfs, "}");
    fprintf(outfs, "$");
  }
  fprintf(outfs, "\\\\\n");

  //print out accuracies:
  fprintf(outfs, " & acc      ");
  ind=0;
  ext=data[7];
  for (int i=1; i<4; i++) {
    if (data[i*2+7] > ext) {
      ext=data[i*2+7];
      ind=i;
    }
  }
  for (int i=0; i<4; i++) {
    float val=data[i*2+7];
    float sd=std[i*2+7];
    fprintf(outfs, " & $");
    if (ind==i) fprintf(outfs, "\\mathbf{");
    fprintf(outfs, "%12.3g", val);
    if (sd>0) fprintf(outfs, "\\pm%9.2g", sd);
    if (ind==i) fprintf(outfs, "}");
    fprintf(outfs, "$");
  }
  fprintf(outfs, "\\\\\n");

  //print out uncertainty coefficients:
  fprintf(outfs, " & U.C.     ");
  ind=0;
  ext=data[8];
  for (int i=1; i<4; i++) {
    if (data[i*2+8] > ext) {
      ext=data[i*2+8];
      ind=i;
    }
  }
  for (int i=0; i<4; i++) {
    float val=data[i*2+8];
    float sd=std[i*2+8];
    fprintf(outfs, " & $");
    if (ind==i) fprintf(outfs, "\\mathbf{");
    fprintf(outfs, "%12.3g", val);
    if (sd>0) fprintf(outfs, "\\pm%9.2g", sd);
    if (ind==i) fprintf(outfs, "}");
    fprintf(outfs, "$");
  }
  fprintf(outfs, "\\\\\n");
  fclose(outfs);

}  

