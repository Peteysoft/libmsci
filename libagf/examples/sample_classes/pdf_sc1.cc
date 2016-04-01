#include "stdio.h"
#include "stdlib.h"

#include "agf_lib.h"
#include "sample_class1_obj.h"

using namespace std;
using namespace libagf;

int main(int argc, char **argv) {
  //file names:
  char *testfile;
  char *outfile;

  //basic data
  sample_class1_obj<real_a> sc;
  real_a **test;		//list of test points
  nel_ta ntest;		//number of test vectors
  dim_ta nvar;		//dimensionality of test data, resp.

  //output
  real_a *pdf;		//probility density

  FILE *fs;		//output file stream
  int nwritten;	//number of data elements written

  if (argc < 3) {
    printf("\nEstimates the probability density (pdf's) of the first sample class\n\n");
    printf("usage:  pdf_agf testfile outfile\n\n");
    printf("     where:\n");
    printf("testfile  = binary file containing list of test points\n");
    printf("outfile   = binary file in which to return the results\n");
    exit(1);
  }

  testfile=argv[1];
  outfile=argv[2];

  test=read_vecfile<real_a>(testfile, ntest, nvar);
  if (test == NULL) {
    fprintf(stderr, "Unable to read input file: %s\n", testfile);
    exit(-2);
  }

  if (nvar != 2) {
    fprintf(stderr, "Sample classes have only two (%d) dimensions, not %d!\n", 2, nvar);
    exit(2);
  }

  pdf=new real_a[ntest]; 
  for (nel_ta i=0; i<ntest; i++) {
    pdf[i]=sc.pdf(test[i][0], test[i][1]);
    printf("%10.4f %10.4f %10.5g\n", test[i][0], test[i][1], pdf[i]);
  }

  fs=fopen(outfile, "w");
  if (fs == NULL) {
    fprintf(stderr, "Unable to open output file for writing: %s\n", outfile);
    exit(-5);
  }
  nwritten=fwrite(pdf, sizeof(real_a), ntest, fs);

  if (nwritten != ntest) {
    fprintf(stderr, "Error writing data: wrote %d out of %d elements\n", nwritten, ntest);
    fclose(fs);
    exit(-6);
  }

  fclose(fs);
  delete [] pdf;

}
  
