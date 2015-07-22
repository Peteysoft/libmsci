#include <assert.h>

#include "sparse.h"

#define MAX_LL 200

extern void conj_grad_norm(sparse_matrix *, float *, float *, float, long);

int main(int argc, char *argv[]) {
  char *infile;		//input file
  char *matfile;	//sparse matrix file
  char *outfile;	//output file

  FILE *infs;		//input file stream
  FILE *matfs;		//sparse matrix stream
  FILE *outfs;		//output file stream

  long n;		//location of the matrix in the file
  sparse_matrix a;	//the sparse matrix
  sparse_matrix at;	//transpose of the matrix
  float *x;		//vector to solve for
  ind_type nx;		//number of elements in x vector
  float *b;		//result vector
  float *b_check;	//check result
  ind_type nb;		//number of elements in b vector

  //used for reading files:
  char line[MAX_LL];	//one line from the file
  long nmat;		//number of matrices in file
  ind_type nn, mm;		//dimensions of matrix

  if (argc < 4 ) return 1;

  infile=argv[1];
  matfile=argv[2];
  outfile=argv[3];

  if (argc >= 5) sscanf(argv[4], "%d", &n);
  else n=1;

  //read in the b vector:
  infs=fopen(infile, "r");
  do {
    fgets(line, MAX_LL, infs);
  } while (line[0] == '#');
  sscanf(line, "%d", &nmat);
  do {
    fgets(line, MAX_LL, infs);
  } while (line[0] == '#');
  sscanf(line, "%d %d", &mm, &nn);
  printf("%d %d\n", mm, nn);
  assert(nn==1);
  nb=mm;
  b=new float[nb];
  for (long i=0; i<nb; i++) {
    fscanf(infs, "%g", b+i);
    printf("%g\n", b[i]);
  }

  //read in the sparse matrix:
  matfs=fopen(matfile, "r");
  for (long i=0; i<n; i++) a.read(matfs);
  fclose(matfs);
  a.dimensions(mm, nx);

  assert(mm==nb);

  //read in the first guess of the x-vector if it exists:
  x=new float[nx];
  if (nmat >=2) {
    do {
      fgets(line, MAX_LL, infs);
    } while (line[0] == '#');
    sscanf(line, "%d %d", &mm, &nn);
    assert (nn == 1);
    assert(nx==mm);
    for (long i=0; i<nx; i++) fscanf(infs, "%g", x+i);
  } else {
    //otherwise use the diagonal of the matrix:
    for (long i=0; i<nx; i++) {
      x[i]=a(i, i);
      printf("x[%d]=%g\n", i, x[i]);
    }
//    a.print(stdout);
  }
  fclose(infs);

  a.transpose(at);

  conj_grad_norm(&at, x, b, 0.00001, 100000L);

  //print the result to a file:
  outfs=fopen(outfile, "w");
  fprintf(outfs, "%d\n", 2);
  fprintf(outfs, "%d %d\n", nx, 1);
  for (long i=0; i<nx; i++) fprintf(outfs, "%g\n", x[i]);
  fprintf(outfs, "#\n");

  //check the result and print it also to a file:
  b_check=new float[nb];
  a.vect_mult(x, b_check);
  fprintf(outfs, "%d %d\n", nb, 1);
  for (long i=0; i<nb; i++) {
    fprintf(outfs, "%g\n", b[i]-b_check[i]);
  }
  fclose(outfs);

}

