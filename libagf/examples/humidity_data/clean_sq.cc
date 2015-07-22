#include <stdio.h>
#include <assert.h>
#include <stdint.h>
#include <math.h>

//fuck all that "structured programming" "software engineering" bullshit!
//we need to get some shit done...

int main() {

  char vecfile[30]="est_bt_cld_land_z0.00.vec";
  char datfile[20]="train_SQ36s.dat";

  FILE *fs;
  int32_t nvar;
  int32_t n1, n;
  int32_t nnew;
  float **x;
  float *y;

  fs=fopen(vecfile, "r");
  fread(&nvar, sizeof(nvar), 1, fs);
  fseek(fs, 0, SEEK_END);
  n1=(ftell(fs)-sizeof(nvar))/sizeof(float)/nvar;
  printf("%d\n", n1);
  fseek(fs, sizeof(nvar), SEEK_SET);

  x=new float *[n1];
  x[0]=new float[n1*nvar];
  for (int i=1; i<n1; i++) x[i]=x[0]+nvar*i;
  
  fread(x[0], sizeof(float), n1*nvar, fs);
  fclose(fs);

  fs=fopen(datfile, "r");
  fseek(fs, 0, SEEK_END);
  n=ftell(fs)/sizeof(float);
  printf("%d\n", n);
  fseek(fs, 0, SEEK_SET);
  assert(n==n1);
  y=new float[n];
  fread(y, sizeof(float), n, fs);
  fclose(fs);

  nnew=0;
  for (int i=0; i<n; i++) {
    if (y[i]>0. && y[i]<1.) {
      printf("%g\n", y[i]);
      y[nnew]=log(y[i]);
      x[nnew]=x[i];
      nnew++;
    }
  }

  fs=fopen("sq_clean_log.vec", "w");
  fwrite(&nvar, sizeof(nvar), 1, fs);
  for (int i=0; i<nnew; i++) fwrite(x[i], sizeof(float), nvar, fs);
  fclose(fs);

  fs=fopen("sq_clean_log.dat", "w");
  fwrite(y, sizeof(float), nnew, fs);
  fclose(fs);

  return 0;

}
  
