#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>

#include "sparse.h"
#include "error_codes.h"

using namespace std;
using namespace libpetey;
using namespace libsparse;

int main(int argc, char *argv[]) {
  FILE *in;
  FILE *out;
  sparse_matrix *a;
  sparse_matrix *prod;
  sparse_matrix read_mat;
  sparse_matrix *exch;
  ind_t maxn=-1, dn=0;
  ind_t i;
  ind_t m, n;
  size_t nread;
  int errcode;

  char c;

  int Dflag=0;
  int Nflag=0;

  while ((c = getopt(argc, argv, "ND")) != -1) {
    switch (c) {
      case ('N'):
             Nflag=1;
             break;
      case ('D'):
             Dflag=1;
             break;
      case ('?'):
             fprintf(stderr, "Unknown option: -%c -- ignored\n", optopt);
             errcode=COMMAND_OPTION_PARSE_ERROR;
             break;
      default:
             fprintf(stderr, "Error parsing command line\n");
             errcode=FATAL_COMMAND_OPTION_PARSE_ERROR;
             exit(errcode);
    }
  }

  argc-=optind;
  argv+=optind;

  if (argc < 2 && (argc < 1 || (Dflag == 0 && Nflag == 0))) {
    printf("\n  purpose: takes an array of sparse matrices and \n");
    printf("       multiplies them, returning one or more \n");
    printf("       sparse matrices\n\n");
    printf("  usage:  sparse_mat_prod infile [-D] [-N] outfile [int [n]]\n\n");
    printf("     where:\n");
    printf("  infile  = input binary file containing an array of\n");
    printf("            sparse matrices\n");
    printf("  outfile = binary output file\n");
    printf("  int     = write interval\n");
    printf("  n       = number of matrices to muliply\n\n");
    printf("  -D      = output dimensions\n\n");
    printf("  -N      = output number of matrices\n\n");
    return 1;
  }

  if (argc > 2) {
    sscanf(argv[2], "%d", &dn);
    if (argc > 3) {
      sscanf(argv[3], "%d", &maxn);
    }
  }

  in=fopen(argv[0], "r");
  read_mat.read(in);

  if (Dflag) {
    read_mat.dimensions(m, n);
    if (m == n) printf("%d\n", m);
    	else printf("%dX%d\n", m, n);
    return 0;
  }

  read_mat.transpose();
  prod=new sparse_matrix(read_mat);
  a=new sparse_matrix();
  a->extend(read_mat.size());

  out=fopen(argv[1], "w");

  i=0;
  while (!feof(in)) {
    nread=read_mat.read(in);
    if (nread == 0) break;
    i++;
    if (Nflag) continue;
    if (dn != 0) if (i % dn == 0) {
      //fprintf(stderr, "Matrix %d is %7.1f%% full\n", i, prod->storage_ratio()*100);
      prod->transpose();
      prod->write(out);
      prod->transpose();
    }
    prod->dimensions(m, n);
    fprintf(stderr, "%d %ld; product matrix is %7.1f%% full; storage ratio: %8.3g\n", 
		i, (int64_t) prod->size(), 100.*(float) prod->size()/(float) (m*n), prod->storage_ratio());
    //read_mat.print(out);
    prod->mat_mult_t(read_mat, *a);
    exch=prod;
    prod=a;
    a=exch;
    //prod->transpose();
    //prod->print(out);
    //prod->transpose();
    if (maxn != -1 && i >= maxn) break;
  }

  fclose(in);

  if (Nflag) {
    printf("%d\n", i+1);
    return 0;
  }

  prod->dimensions(m, n);
  prod->transpose();
  prod->write(out);
  //prod->print(out);
  fclose(out);

  delete prod;
  delete a;

  return 0;
}
