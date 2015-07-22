#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>

//#include "/home/Petey/software_projects/include/class_lib.h"

int main(int argc, char ** argv) {
  FILE *fs;

  long n;

  float *con;

  float *con1;

  float ave, ave1;
  float std, std1;
  float r;

//  long *clind;
//  long ncl;

  //read in the training data:
  fs=fopen(argv[1], "r");
  fseek(fs, 0, SEEK_END);
  n=ftell(fs)/4;

  fseek(fs, 0, SEEK_SET);
  con=new float [n];
  fread(con, sizeof(float), n, fs);
  fclose(fs);

  fs=fopen(argv[2], "r");
  fseek(fs, 0, SEEK_END);
  assert(n == ftell(fs)/4);
  rewind(fs);
  con1=new float [n];
  for (long i=0; i<n; i++) {
    fread(con1+i, sizeof(long), 1, fs);
  }

  //it's a hack, but fix 'nan's and other nuisances:
  for (long i=0; i<n; i++) if (finite(con[i])==0) con[i]=1;
  for (long i=0; i<n; i++) if (finite(con1[i])==0) con1[i]=1;

  fclose(fs);


  for (long i=0; i<n; i++) {
    printf("%8.3f %8.3f\n", con[i], con1[i]);
  }

  //calculate accuracy:
  ave=0;
  ave1=0;
  for (long i=0; i<n; i++) {
    ave+=con[i];
    ave1+=con1[i];
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

  printf("correlation: %f \n", r);

  delete [] con;
  delete [] con1;

}
