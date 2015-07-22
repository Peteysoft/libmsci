#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <getopt.h>

#include "error_codes.h"
#include "read_ascii_all.h"

using namespace libpetey;

int main(int argc, char **argv) {
  int fileid=getpid();		//pid alone should be sufficient to prevent collision
  char fname[15];			//in/out file name
  char *command;			//run this command
  FILE *fs;				//in/out file-stream
  int D;				//number of features
  int32_t lab1, lab2, cls;		//labels, class
  float p1, p2;				//conditional probabilities
  float R;
  int err;
  int pflag=0;		//calculate probability estimates
  int Kflag=0;		//don't delete temporary file
  int Nflag=0;		//return number of dimensions
  int c;

  while ((c=getopt(argc, argv, "pKN")) != -1) {
    switch(c) {
      case ('p'):
        pflag=1;
        break;
      case ('K'):
        Kflag=1;
        break;
      case ('N'):
        Nflag=1;
        break;
      default:
        fprintf(stderr, "svm_1shot_binary: Error parsing command line\n");
        exit(FATAL_COMMAND_OPTION_PARSE_ERROR);
    }
  }
  argc-=optind;
  argv+=optind;

  if (argc<1 || argc<2 && Nflag==0) {
    printf("svm_1shot_binary [-p] [-K] model x1 [x2 [x3 ...]]\n");
    if (argc==1) exit(0); else exit(INSUFFICIENT_COMMAND_ARGS);
  }

  //find out dimensionality of the data:
  if (Nflag) {
    char *line=NULL;
    FILE *fs=fopen(argv[0], "r");
    int32_t nvar;
    do {
      line=fget_line(fs, 1);
      if (strcmp(line, "SV")==0) break;
      if (line!=NULL) {
        delete [] line;
      } else {
        fprintf(stderr, "svm_1shot_binary: error in model file format\n");
        exit(FILE_READ_ERROR);
      }
    } while (feof(fs)==0);
    line=fget_line(fs, 1);
    nvar=0;
    for (int i=0; line[i]!='\0'; i++) if (line[i]==':') nvar++;
    delete [] line;
    printf("%d\n", nvar);
    exit(0);
  }


  //write the feature data from the command line to a file:
  D=argc-1;
  sprintf(fname, "in.%6.6d.svm", fileid);
  fs=fopen(fname, "w");
  if (fs==NULL) {
    fprintf(stderr, "svm_1shot_binary: unable to open input file, %s\n", fname);
    exit(UNABLE_TO_OPEN_FILE_FOR_WRITING);
  }
  fprintf(fs, "%d", 1);
  for (int i=0; i<D; i++) fprintf(fs, " %d:%s", i+1, argv[i+1]);
  fclose(fs);

  //create and run the command that calls LIBSVM:
  command=new char [strlen(argv[0])+50];
  sprintf(command, "svm-predict -q -b %d in.%6.6d.svm %s out.%6.6d.svm", pflag, fileid, argv[0], fileid);
  err=system(command);
  if (err!=0) {
    fprintf(stderr, "svm_1shot_binary: system call,\n");
    fprintf(stderr, "  \"%s\"\n", command);
    fprintf(stderr, "  returned error code, %d\n", err);
    system(command);
    exit(err);
  }

  //delete input file:
  if (Kflag==0) {
    sprintf(command, "rm -f in.%6.6d.svm", fileid);
    system(command);
  }
 
  //read data in output file: 
  sprintf(fname, "out.%6.6d.svm", fileid);
  fs=fopen(fname, "r");
  if (fs==NULL) {
    fprintf(stderr, "svm_1shot_binary: unable to open output file, %s\n", fname);
    exit(UNABLE_TO_OPEN_FILE_FOR_READING);
  }
  //if there are probability estimates, need to read those and print the difference:
  if (pflag) {
    fseek(fs, 6, SEEK_CUR);
    err=fscanf(fs, "%d %d %d %g %g", &lab1, &lab2, &cls, &p1, &p2);
    fclose(fs);
    if (err!=5) {
      fprintf(stderr, "svm_1shot_binary: failed to read sufficient data from output file, %s\n", fname);
      exit(FILE_READ_ERROR);
    }
    if (lab1 < lab2) {
      R=p2-p1;
    } else {
      R=p1-p2;
    }
  } else {
    //otherwise just read the class and convert it to sign:
    err=fscanf(fs, "%d", &cls);
    if (cls==1) R=1; else R=-1;
    fclose(fs);
    if (err!=1) {
      fprintf(stderr, "svm_1shot_binary: failed to read sufficient data from output file, %s\n", fname);
      exit(FILE_READ_ERROR);
    }
  }

  //delete output file:
  if (Kflag==0) {
    sprintf(command, "rm -f out.%6.6d.svm", fileid);
    system(command);
  }

  printf("%g", R);
  delete [] command;

  return 0;

} 
