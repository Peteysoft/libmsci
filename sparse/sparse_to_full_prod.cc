#include <stdio.h>

#include "sparse.h"
#include "full_util.h"

using namespace libpetey;
using namespace libsparse;

int main(int argc, char *argv[]) {
  FILE *in;
  FILE *out;
  ind_t m, n;
  ind_t mr, nr;		//dimensions of the result
  float **a;
  float **b;
  sparse_matrix read_mat;
  float **exch;
  int32_t maxn=-1, dn=0;
  int32_t i;
  size_t nread;

  if (argc < 3) {
    printf("\n  purpose: takes an array of sparse matrices and \n");
    printf("           multiplies them to one or more \n");
    printf("           full matrices\n\n");
    printf("  usage:  sparse_mat_prod infile outfile [int [n]]\n\n");
    printf("      where:\n");
    printf("  infile  = input binary file containing an array of\n");
    printf("            sparse matrices\n");
    printf("  outfile = binary output file\n");
    printf("  int     = write interval\n");
    printf("  n       = number of matrices to muliply\n\n");
    return 1;
  }

  if (argc > 3) {
    sscanf(argv[3], "%d", &dn);
    if (argc > 4) {
      sscanf(argv[4], "%d", &maxn);
    }
  }

  in=fopen(argv[1], "r");
  read_mat.read(in);
  read_mat.dimensions(m, n);
  if (m>n) {
    mr=m;
    nr=n;
  } else {
    mr=n;
    nr=n;
  }
  a=allocate_matrix<float,ind_t>(mr, nr);
  b=allocate_matrix<float,ind_t>(mr, nr);
  identity_matrix(b, mr, nr);
  read_mat.mat_mult(b, a, nr);

  out=fopen(argv[2], "w");

  i=0;
  while (!feof(in)) {
    if (dn != 0) if (i % dn == 0) {
      write_matrix(out, a, mr, nr);
    }
    i++;
    printf("%d\n", i);
    nread=read_mat.read(in);
    if (nread == 0) break;
    read_mat.mat_mult(a, b, nr);
    exch=a;
    a=b;
    b=exch;
    if (maxn != -1 && i >= maxn) break;
  }

  fclose(in);

  write_matrix(out, a, m, n);
  fclose(out);

  delete_matrix(a);
  delete_matrix(b);

  return 0;
}
