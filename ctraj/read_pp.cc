#include <stdio.h>
#include <stdint.h>

#include "error_codes.h"
#include "full_util.h"
#include "time_class.h"
#include "parse_command_opts.h"

#include "ctraj_defaults.h"

#define HEADLEN 64

void swap_endian (int32_t *data, int n) {
  char *fourbyte;
  char swp;
  for (int i=0; i<n; i++) {
    fourbyte=(char *) (data+i);
    swp=fourbyte[0];
    fourbyte[0]=fourbyte[3];
    fourbyte[3]=swp;
    swp=fourbyte[1];
    fourbyte[1]=fourbyte[2];
    fourbyte[2]=swp;
  }
}

int pp_read_all(char *fname, int32_t **headers_all, float ***fields, int nmax) {

  FILE *fs;
  int32_t f77recsep;
  int32_t *header;
  float **data;
  int readcount;
  int nrec;

  fs=fopen(fname, "r");
  if (fs==NULL) {
    fprintf(stderr, "Error: unable to open %s for input\n", fname);
    throw UNABLE_TO_OPEN_FILE_FOR_READING;
  }

  for (nrec=0; nrec<nmax; nrec++) {
    readcount=fread(&f77recsep, sizeof(f77recsep), 1, fs);
    printf("%d %d\n", nrec, readcount);
    if (readcount==0) {
      break;
    }
    swap_endian(&f77recsep, 1);

    //should be 256:
    //printf("%d\n", f77recsep);
    if (f77recsep != HEADLEN*4) {
      fprintf(stderr, "Error: error in record separator--wrong file type\n");
      fprintf(stderr, " actual: %d; expected: %d\n", f77recsep, HEADLEN*4);
      throw FILE_READ_ERROR;
    }

    header=new int32_t[HEADLEN];
    fread(header, sizeof(int32_t), HEADLEN, fs);
    swap_endian(header, HEADLEN);
    headers_all[nrec]=header;

    fread(&f77recsep, sizeof(f77recsep), 1, fs);
    fread(&f77recsep, sizeof(f77recsep), 1, fs);
    swap_endian(&f77recsep, 1);

    //should be 6912*4:
    printf("%d\n", f77recsep);

    if (header[14]!=header[17]*header[18]) {
      fprintf(stderr, "Error in record header: %d*%d != %d\n", header[17], header[18], header[14]);
      throw FILE_READ_ERROR;
    }

    data=allocate_matrix<float, int32_t>(header[17], header[18]);

    fread(data[0], sizeof(float), header[14], fs);
    swap_endian((int32_t *) data, header[14]);
    fields[nrec]=data;

    fread(&f77recsep, sizeof(f77recsep), 1, fs);

  }

  fclose(fs);

  return nrec;

}

#define MAXNREC 1000

int main(int argc, char **argv) {
  float **data[MAXNREC];
  int32_t *header[MAXNREC];
  int nrec;
  
  nrec=pp_read_all(argv[1], header, data, MAXNREC);

  for (int i=0; i<nrec; i++) {
    printf("Record %d:\n", i);
    //print out first 38 header records:
    for (int j=0; j<38; j++) printf("%d %d\n", j, header[i][j]);

    //print out 50-64:
    for (int j=49; j<HEADLEN; j++) printf("%d %g\n", j, ((float *) header[i])[j]);
    printf("\n");
  }
  printf("\n");

  for (int i=0; i<nrec; i++) {
    printf("deleting %dth record\n", i);
    delete [] header[i];
    delete_matrix(data[i]);
  }

}
