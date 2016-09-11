#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "error_codes.h"
#include "parse_command_opts.h"
#include "ctraj_defaults.h"

using namespace libpetey;
using namespace ctraj;

float correlate(float *v1, float *v2, int32_t n) {
  float ave1, ave2;
  float var1, var2, cov;
  float diff1, diff2;

  ave1=0;
  ave2=0;
  cov=0;
  for (int32_t i=0; i<n; i++) {
    ave1+=v1[i];
    ave2+=v2[i];
  }
  ave1/=n;
  ave2/=n;

  var1=0;
  var2=0;
  for (int32_t i=0; i<n; i++) {
    diff1=v1[i]-ave1;
    var1+=diff1*diff1;
    diff2=v2[i]-ave2;
    var2+=diff2*diff2;
    cov+=diff1*diff2;
  }

  return cov/sqrt(var1*var2);
}

int main(int argc, char **argv) {
  FILE *fs1, *fs2;
  size_t ncon;
  float *field1, *field2;
  int32_t nvar1, nvar2; 
  int32_t index1=-1, index2=-1;
  int32_t n1, n2;
  int32_t ind0=0, N;
  size_t fsize;
  char c;
  size_t recbytes;
  void *optarg[10];
  int flag[10];

  optarg[0]=&index1;
  optarg[1]=&index2;
  //optarg[2]=&ind0=0;
  //optarg[3]=&N;

  argc=parse_command_opts(argc, argv, "12?", "%d%d%", optarg, flag, OPT_WHITESPACE);
  if (argc<0) exit(21);

  if (argc < 3 || flag[2]) {
    FILE *docfs;
    int err;
    if (flag[2]) {
      docfs=stdout;
      err=0;
    } else {
      docfs=stderr;
      err=INSUFFICIENT_COMMAND_ARGS;
    }

    fprintf(docfs, "Usage: correlate_field [-1 index1] [-2 index2] file1 file2\n");
    ctraj_optargs(docfs, "12?");
    return err;
  }

  fs1=fopen(argv[1], "r");
  fs2=fopen(argv[2], "r");

  fread(&nvar1, sizeof(nvar1), 1, fs1);
  fread(&nvar2, sizeof(nvar2), 1, fs2);

  if (nvar1!=nvar2) {
    fprintf(stderr, "correlate_fields: vector lengths don't match; %d vs %d\n", nvar1, nvar2);
    exit(401);
  }

  recbytes=nvar1*sizeof(float);
  fseek(fs1, 0, SEEK_END);
  fsize=ftell(fs1);
  n1=(fsize-sizeof(nvar1))/recbytes;

  if ((fsize-sizeof(nvar1))%recbytes!=0) {
    fprintf(stderr, "%s: not an even number of records;\n", argv[0]);
    fprintf(stderr, "     %ld mod %ld = %ld\n", fsize, recbytes, fsize%recbytes);
  }

  fseek(fs2, 0, SEEK_END);
  fsize=ftell(fs2);
  n2=(fsize-sizeof(nvar1))/recbytes;

  if ((fsize-sizeof(nvar1))%recbytes!=0) {
    fprintf(stderr, "%s: not an even number of records;\n", argv[1]);
    fprintf(stderr, "     %ld mod %ld = %ld\n", fsize, recbytes, fsize%recbytes);
  }

  if (index1 >= n1 || index2 >= n2) {
    fprintf(stderr, "Index out of range; %s: %d; %s: %d\n", argv[0], n1, argv[1], n2);
    exit(2);
  }

  field1=new float[nvar1];
  field2=new float[nvar1];

  if (index1 < 0 && index2 < 0) {
    fseek(fs1, sizeof(nvar1), SEEK_SET);
    fseek(fs2, sizeof(nvar1), SEEK_SET);
    if (n1 > n2) n1=n2;
    for (int32_t i=ind0; i<n1; i++) {
      fread(field1, sizeof(float), nvar1, fs1);
      fread(field2, sizeof(float), nvar1, fs2);
      printf("%g\n", correlate(field1, field2, nvar1));
    }
  } else {
    if (index1<0) index1=0;
    if (index2<0) index2=0;
    fseek(fs1, index1*recbytes+sizeof(nvar1), SEEK_SET);
    fseek(fs2, index2*recbytes+sizeof(nvar1), SEEK_SET);
    fread(field1, sizeof(float), nvar1, fs1);
    fread(field2, sizeof(float), nvar2, fs2);
    printf("%g\n", correlate(field1, field2, nvar1));
  }
  
  delete [] field1;
  delete [] field2;

  fclose(fs1);
  fclose(fs2);

  return 0;
}  

