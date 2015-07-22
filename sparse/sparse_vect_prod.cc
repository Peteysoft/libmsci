#include <stdio.h>
#include <stdlib.h>

#include "parse_command_opts.h"
#include "sparse.h"
#include "error_codes.h"

using namespace std;
using namespace libpetey;
using namespace libsparse;

int main(int argc, char *argv[]) {
  sparse_matrix mat;
  char *vect_file;
  char *mat_file;
  FILE *vect_fs;
  FILE *mat_fs;
  ind_t m, n;
  int32_t nvar;
  int32_t i0, i;
  int32_t N=-1;
  float *vect1, *vect2, *swp;
  size_t ncon; 
  int errcode;
  void *optarg[3];
  int flags[3];

  errcode=0;
  i0=0;

  optarg[0]=&i0;
  optarg[1]=&N;
  argc=parse_command_opts(argc, argv, "in", "%d%d", optarg, flags, 1);
  if (argc < 0) exit(FATAL_COMMAND_OPTION_PARSE_ERROR);

  if (argc < 3) {
    printf("\n  purpose: multiplies a vector with an array of sparse matrices\n");
    printf("\n");
    printf("  usage:  sparse_vect_prod [-n n] [-i i0] infile outfile\n\n");
    printf("      where:\n");
    printf("  infile  = binary input file containing an array of\n");
    printf("            sparse matrices (should be square...)\n");
    printf("  outfile = binary output file containing initial vector\n");
    printf("            and array of output vectors\n");
    printf("  n       = number of matrices to muliply\n\n");
    printf("  i0      = starting point\n\n");
    return 1;
  }

  vect_file=argv[2];
  mat_file=argv[1];

  mat_fs=fopen(mat_file, "r");
  mat.read(mat_fs);
  mat.dimensions(m, n);

  printf("%d %d\n", m, n);

  vect1=new float[n];
  vect2=new float[n];
  
  vect_fs=fopen(vect_file, "r+");
  //if we add a header containing vector size, files are compatible
  //with libagf:
  fread(&nvar, sizeof(nvar), 1, vect_fs);
  if (n!=nvar) {
    fprintf(stderr, "sparse_vect_prod: Matrix inner dimension and vector dimension do not match: %d vs. %d\n", n, nvar);
    exit(411);
  }

  fread(vect2, sizeof(float), n, vect_fs);

  for (i=0; i<i0; i++) {
    if (mat.read(mat_fs)==0) {
      errcode=FILE_READ_ERROR;
      fprintf(stderr, "Error reading file: %s\n", mat_file);
      exit(errcode);
    }
  }

  mat.vect_mult(vect2, vect1);
  fwrite(vect1, sizeof(float), n, vect_fs);

  //for repeated mappings, need square matrix:
  if (n != m) {
    fclose(vect_fs);
    fclose(mat_fs);
    exit(410);
  }

  i=i0;

  while (!feof(mat_fs)) {
    printf("%d\n", i);
    if (flags[1] && i-i0 >= N) break;
    i++;
    if (mat.read(mat_fs)==0) break;
    mat.vect_mult(vect1, vect2);
    fwrite(vect2, sizeof(float), n, vect_fs);
    swp=vect1;
    vect1=vect2;
    vect2=swp;
  }

  delete [] vect1;
  delete [] vect2;

  fclose(mat_fs);
  fclose(vect_fs);
  
}

