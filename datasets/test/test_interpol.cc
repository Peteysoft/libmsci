#include <stdio.h>

#include "time_class.h"
#include "simple_temp.h"
#include "dependent_temp.h"
#include "string_petey.h"

simple<time> *tgrid;
simple<float> *lon;
simple<float> *lat;
simple<float> *theta;
dependent<float> *u;
dependent<float> *v;

void read_wind_file(char *filename) {
  FILE *lun;
  //long magic;
  long ndim;
  long nvar;
  long type;
  string varname;
  float *data;
  long n;
  short y;
  short mon;
  short d;
  short min;
  short h;
  float s;
  time *tdata;
  
  printf("Opening file %s\n", filename);
  lun=fopen(filename, "r");  
  //fread(&magic, sizeof(long), 1, lun);
  fread(&ndim, sizeof(long), 1, lun);
  printf("%d dimensions found in file\n", ndim);
  
  varname.read(lun);
  varname.print();
  fread(&type, sizeof(long), 1, lun);
  fread(&n, sizeof(long), 1, lun);
  printf(": %d elements of type %d\n", n, type);
  data=new float[n];
  fread(data, sizeof(long), n, lun);
  lon=new simple<float>(data, n);
  delete [] data;
  
  varname.read(lun);
  varname.print();
  fread(&type, sizeof(long), 1, lun);
  fread(&n, sizeof(long), 1, lun);
  printf(": %d elements of type %d\n", n, type);
  data=new float[n];
  fread(data, sizeof(long), n, lun);
  lat=new simple<float>(data, n);
  delete [] data;

  varname.read(lun);
  varname.print();
  fread(&type, sizeof(long), 1, lun);
  fread(&n, sizeof(long), 1, lun);
  printf(": %d elements of type %d\n", n, type);
  data=new float[n];
  fread(data, sizeof(long), n, lun);
  theta=new simple<float>(data, n);
  delete [] data;
  
  varname.read(lun);
  varname.print();
  fread(&type, sizeof(long), 1, lun);
  fread(&n, sizeof(long), 1, lun);
  printf(": %d elements of type %d\n", n, type);
  tdata=new time[n];
  for (long i=0; i<n; i++) {
    fread(&y, sizeof(short), 1, lun);
    fread(&mon, sizeof(short), 1, lun);
    fread(&d, sizeof(short), 1, lun);
    fread(&h, sizeof(short), 1, lun);
    fread(&min, sizeof(short), 1, lun);
    fread(&s, sizeof(float), 1, lun);
    tdata[i].init(y, mon, d, h, min, s);
  }
  tgrid=new simple<time>(tdata, n);
  delete [] tdata;
  
  fread(&nvar, sizeof(long), 1, lun);
  printf("%d variables found in file\n");
  
  fread(&type, sizeof(long), 1, lun);
  varname.read(lun);
  varname.print();
  printf(": type %d\n", type);
  u=new dependent<float>(lon, lat, theta, tgrid);
  u->read(lun);
  
  fread(&type, sizeof(long), 1, lun);
  varname.read(lun);
  varname.print();
  printf(": type %d\n", type);
  v=new dependent<float>(lon, lat, theta, tgrid);
  v->read(lun);

  fclose(lun);
}  
  
int main(int argc, char *argv[]) {
  float x;
  float y;
  char tstring[100];
  time t;
  double int1, int2, int3;
  float ui, vi;
    
  read_wind_file(argv[1]);
  
  do {
    scanf("%g", &x);
    scanf("%g", &y);
    scanf("%s", tstring);
  
    t.read_string(tstring);
    
    int1=lon->interp(x);
    int2=lat->interp(y);
    int3=tgrid->interp(t);
    
    u->interpol(ui, int1, int2, 0., int3);
    v->interpol(vi, int1, int2, 0., int3);

    printf("%g %g\n", ui, vi);
  } while (1);
  
  delete lon;
  delete lat;
  delete theta;
  delete tgrid;

  delete u;
  delete v;
    
}