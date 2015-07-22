#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>

//#include "/home/Petey/software_projects/include/class_lib.h"

int main(int argc, char ** argv) {
  char *confile;
  char *clsfile;
  char line[200];

  FILE *fs;

  long n;

  float *con;
  long *cls;

  float *con1;
  long *cls1;

  float p1, p2;

  float ave, ave1;
  float std, std1;
  float r;

  long nsame;

//  long *clind;
//  long ncl;

  //generate the input file names:
  confile=new char[strlen(argv[1])+5];
  strcpy(confile, argv[1]);
  strcat(confile, ".con");

  clsfile=new char[strlen(argv[1])+5];
  strcpy(clsfile, argv[1]);
  strcat(clsfile, ".cls");

  //read in the training data:
  fs=fopen(confile, "r");
  fseek(fs, 0, SEEK_END);
  n=ftell(fs)/4;

  fseek(fs, 0, SEEK_SET);
  con=new float [n];
  fread(con, sizeof(float), n, fs);
  fclose(fs);

  fs=fopen(clsfile, "r");
  fseek(fs, 0, SEEK_END);
  assert(n == ftell(fs)/4);
  rewind(fs);
  cls=new long [n];
  for (long i=0; i<n; i++) {
    fread(cls+i, sizeof(long), 1, fs);
  }

  //it's a hack, but fix 'nan's and other nuisances:
  for (long i=0; i<n; i++) if (finite(con[i])==0) con[i]=1;

  fclose(fs);

  fs=fopen(argv[2], "r");

  //skip header:
  //fgets(line, 200, fs);
  cls1=new long[n];
  con1=new float[n];
  for (long i=0; i<n; i++) {
    fgets(line, 200, fs);
    sscanf(line, "%d %g %g", cls1+i, &p1);  //, &p2);
    //con1[i]=fabs(p1-p2);
    con1[i]=p1;
  }

  fclose(fs);

  for (long i=0; i<n; i++) {
    printf("%d %d %8.3f %8.3f\n", cls[i], cls1[i], con[i], con1[i]);
  }

  //calculate accuracy:
  ave=0;
  ave1=0;
  nsame=0;
  for (long i=0; i<n; i++) {
    ave+=con[i];
    ave1+=con1[i];
    if (cls[i] == cls1[i]) nsame++;
  }

  ave/=n;
  ave1/=n;

  r=0;
  std=0;
  std1=0;
  for (long i=0; i<n; i++) {
    float diff=con[i]-ave;
    float diff1=con1[i]-ave1;
    r+=diff*diff1;
    std+=diff*diff;
    std1+=diff1*diff1;
  }
  std=sqrt(std/(n-1));
  std1=sqrt(std1/(n-1));
  printf("%f %f %f %f %f\n", ave, ave1, std, std1, r);

  r=r/std/std1/(n-1);

  printf("similarity: %f \%\n", 1.*nsame/n);
  printf("correlation: %f \n", r);

  delete [] con;
  delete [] cls;
  delete [] con1;
  delete [] cls1;

}
