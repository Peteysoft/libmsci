#include <math.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>

#include "class_lib.h"

float **global_train;
long *global_clind;
long global_nvar;

int main(int argc, char *argv[]) {

  char *vecfile;		//vector training data
  char *classfile;		//class training data
  char *brdfile;	//binary file sampling border
  char *grdfile;	//binary file containing gradient vectors

  long ntrain;		//number of training data points
  int32_t nvar;		//number of variables

  long *cls;		//training data classes

  FILE *ifs;		//input file stream

  long nclass;		//number of classes (anything above two is ignored though...)

  float *all;		//for reading the training vectors

  if (argc != 2) {
    printf("syntax:  a.out train\n");
    printf("\n");
    printf("  where:\n");
    printf("train   = binary files containing locations of the samples:\n");
    printf("          .vec for vectors;  .cls for classes\n");
    return 1;
  }

  //parse the command line arguments:
  vecfile=new char[strlen(argv[1])+5];
  strcpy(vecfile, argv[1]);
  strcat(vecfile, ".vec");

  classfile=new char[strlen(argv[1])+5];
  strcpy(classfile, argv[1]);
  strcat(classfile, ".cls");

  brdfile=new char[strlen(argv[1])+7];
  strcpy(brdfile, argv[1]);
  strcat(brdfile, ".s.vec");

  grdfile=new char[strlen(argv[1])+7];
  strcpy(grdfile, argv[1]);
  strcat(grdfile, ".s.cls");

  //read in the training data:
  ifs=fopen(vecfile, "r");
  fread(&global_nvar, sizeof(global_nvar), 1, ifs);
  fseek(ifs, 0, SEEK_END);
  ntrain=(ftell(ifs)-4)/global_nvar/4;

  printf("%ld training vectors found in file %s\n", ntrain, vecfile);

  fseek(ifs, 4, SEEK_SET);
  all=new float [ntrain*global_nvar];
  fread(all, sizeof(float), ntrain*global_nvar, ifs);
  global_train=new float *[ntrain];
  for (long i=0; i<ntrain; i++) global_train[i]=&all[i*global_nvar];

  fclose(ifs);

  ifs=fopen(classfile, "r");

  fseek(ifs, 0, SEEK_END);
  assert(ntrain == ftell(ifs)/4);

  rewind(ifs);
  cls=new long [ntrain];
  nclass=0;
  for (long i=0; i<ntrain; i++) {
    fread(&cls[i], sizeof(long), 1, ifs);
    if (cls[i] >= nclass) nclass=cls[i]+1;
  }

  fclose(ifs);

  //sort the samples into classes:
  for (long i=0; i<ntrain; i++) printf("%ld ", cls[i]);
  printf("\n");
  global_clind=sort_classes(global_train, global_nvar, ntrain, cls, nclass);
  for (long i=0; i<ntrain; i++) printf("%ld ", cls[i]);
  printf("\n");
  for (long i=0; i<=nclass; i++) printf("%ld ", global_clind[i]);
  printf("\n");
  
  ifs=fopen(brdfile, "w");
  nvar=global_nvar;
  fwrite(&nvar, 1, sizeof(nvar), ifs);
  for (long i=0; i<ntrain; i++) {
    fwrite(global_train[i], 1, sizeof(float)*global_nvar, ifs);
  }
  fclose(ifs);
  
  ifs=fopen(grdfile, "w");
  fwrite(cls, 1, sizeof(long)*ntrain, ifs);
  fclose(ifs);

}

  
