#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include "agf_lib.h"

#define OPT_VER ""

using namespace std;
using namespace libagf;

int main(int argc, char *argv[]) {
  FILE *fs;
  FILE *diagfs;

  char *basename;		//base name of _ALL_ the files

  char *trainbase;		//base name of separated training data
  char *trainvecfile;		//for finding the class borders

  char *trainclsfile;		//	"
  char *testbase;
  char *testfile;		//test data
  char *outfile;		//output classes
  char *confile;		//output confidence ratings
  char *brdfile;		//contains class borders
  int nchar;			//number of characters in filename

  nel_ta nsamp;		//total number of samples
  dim_ta nvar;		//number of variables
  cls_ta nclass;	//number of classes
  nel_ta ntrain;	//number of training data points
  nel_ta ntest;		//number of test data points
  dim_ta nvar1;

  real_a **train;	//training data vectors
  cls_ta *cls;		//training data classes
  real_a **test;		//test data vectors
  cls_ta *result;		//results of classification
  real_a *con;

  real_a *all;		//for memory allocation

  nel_ta ntrue;
  nel_ta nfalse;
  real_a corr;

  char command[200];	//the system command

  agf_command_opts opt_args;
  char normflag[3];

  int exit_code;
  int err_code;

  diagfs=stdout;

  //set defaults and parse command line options:
  opt_args.k=-1;
  opt_args.W2=W_DEFAULT_BORDERS;
  opt_args.n=NBORD_DEFAULT;
  opt_args.tol=TOL_DEFAULT;
  opt_args.ftest=F_DEFAULT;

  exit_code=0;
  exit_code=agf_parse_command_opts(argc, argv, "h:i:I:l:f:k:W:s:t:v:V:n", &opt_args);
  if (exit_code==FATAL_COMMAND_OPTION_PARSE_ERROR) return exit_code;

  if (argc < 1) {
    printf("\n");
    printf("syntax:  class_borders [-n] [-f ftest] [-k k] [-w Wc] [-s n] [-t tol] train\n");
    printf("\n");
    printf("arguments:\n");
    printf("  train    = binary files containing locations of the samples:\n");
    printf("              .vec for vectors;\n");
    printf("              .cls for classes\n");
    printf("\n");
    printf("options:\n");
    printf("  -f ftest = fraction of training data to use for testing (default=%g)\n", opt_args.ftest);
    printf("  -v var1  = lower filter variance bracket\n");
    printf("               --default is to use the total variance/n^(2/D)\n");
    printf("  -V var2  = upper filter variance bracket/initial filter variance\n");
    printf("               --default is to use the total variance of the data\n");
    printf("  -k k     = number of nearest neighbours to use in each calculation\n");
    printf("               --default is to use all the data\n");
    printf("  -W Wc    = objective total weight (default=%g)\n", opt_args.W2);
    printf("  -s n     = number of times to sample the border (default=%d)\n", opt_args.n);
    printf("  -t tol   = desired tolerance (default=%g)\n", opt_args.tol);
    printf("  -n       = option to normalise the data\n");
    printf("\n");
    return INSUFFICIENT_COMMAND_ARGS;
  }

  if (opt_args.normflag == 1) {
    strcpy(normflag, "-n");
  } else {
    strcpy(normflag, "");
  }

  //figure out all the various file names:
  basename=argv[0];
  nchar=strlen(basename);
  
  err_code=agf_read_train(argv[0], train, cls, nsamp, nvar);
  if (err_code!=0) return err_code;
  all=train[0];		//must do this because of the way the data is allocated...

  //count the number of classes:
  nclass=1;
  for (nel_ta i=0; i<nsamp; i++) if (cls[i]>=nclass) nclass++;
  
  //file name for the separated data (training and test):
  trainbase=new char[nchar+5];
  strcpy(trainbase, basename);
  strcat(trainbase, ".trn");
  
  trainvecfile=new char[nchar+9];
  strcpy(trainvecfile, trainbase);
  strcat(trainvecfile, ".vec");
  
  trainclsfile=new char[nchar+9];
  strcpy(trainclsfile, trainbase);
  strcat(trainclsfile, ".cls");
  
  testbase=new char[nchar+6];
  strcpy(testbase, basename);
  strcat(testbase, ".test");
  
  testfile=new char[nchar+10];
  strcpy(testfile, testbase);
  strcat(testfile, ".vec");
  
  //files containing class borders:
  brdfile=new char[nchar+5];
  strcpy(brdfile, basename);
  strcat(brdfile, ".brd");
  
  //output data (the stuff we need!):
  outfile=new char[nchar+10];
  strcpy(outfile, testbase);
  strcat(outfile, ".cls");
  
  confile=new char[nchar+10];
  strcpy(confile, testbase);
  strcat(confile, ".con");

  fprintf(diagfs, "%d training vectors found: %s\n", nsamp, basename);

/*
  for (long i=0; i<nsamp; i++) {
    for (long j=0; j<nvar; j++) printf("%f ", train[i][j]);
    printf("%d\n", cls[i]);
  }
*/

  //randomize the training data:
  fprintf(diagfs, "Randomizing the data...\n");
  randomize_vec(train, nvar, nsamp, cls);

  //calculate the number of test vectors:
  ntest=(nel_ta) (opt_args.ftest*nsamp);
  ntrain=nsamp-ntest;
  fprintf(diagfs, "Using %d vectors for testing\n", ntest);

  if (opt_args.k <= opt_args.W2 || opt_args.k >= ntrain) {
    if (opt_args.k != -1) {
      fprintf(stderr, "Parameter k=%d out of range.  Using all the training data.\n", opt_args.k);
      opt_args.k=-1;
      exit_code=PARAMETER_OUT_OF_RANGE;
    }
  }

  //now we need to output the training and test data to separate files:
  fs=fopen(trainvecfile, "w");
  fwrite(&nvar, 1, sizeof(nvar), fs);
  for (nel_ta i=ntest; i<nsamp; i++) {
    fwrite(train[i], sizeof(real_a), nvar, fs);
  }
  fclose(fs);
  
  fs=fopen(trainclsfile, "w");
  fwrite(cls+ntest, sizeof(cls_ta), ntrain, fs);
  fclose(fs);
  
  fs=fopen(testfile, "w");
  fwrite(&nvar, 1, sizeof(nvar), fs);
  for (nel_ta i=0; i<ntest; i++) {
    fwrite(train[i], sizeof(real_a), nvar, fs);
  }
  fclose(fs);
  
  //begin the classification scheme:
  fprintf(diagfs, "Beginning classification...\n");
  //now we run the program to find the class borders:
  //generate the command:
  sprintf(command, "time class_borders%s %s -h %ld -i %ld -l %g -v %g -V %g -k %d -W %f -s %d -t %g %s %s", 
		OPT_VER, normflag, 
		agf_global_borders_maxiter, agf_global_weights_maxiter,
		agf_global_weights_tol,
		opt_args.var[0], opt_args.var[1], opt_args.k, opt_args.W2, opt_args.n, opt_args.tol, 
		trainbase, brdfile);
  fprintf(diagfs, "%s\n", command);
  fprintf(diagfs, "%d\n", err_code=system(command));
  if (err_code!=0) exit(err_code);
  
  //now we use the class borders so generated to classify the test data:
  sprintf(command, "time classify_b%s %s %s %s %s", OPT_VER, normflag,
		  brdfile, testfile, testbase);
  fprintf(diagfs, "%s\n", command);
  fprintf(diagfs, "%d\n", err_code=system(command));
  if (err_code!=0) exit(err_code);

  //read in the results of our efforts:
  result=new cls_ta[ntest];
  con=new real_a[ntest];
  
  fs=fopen(outfile, "r");
  fread(result, sizeof(cls_ta), ntest, fs);
  fclose(fs);  

  //no confidence rating yet!!  
  fs=fopen(confile, "r");
  fread(con, 1, sizeof(real_a)*ntest, fs);
  fclose(fs);  

  ntrue=0;
  for (nel_ta i=0; i<ntest; i++) if (result[i]==cls[i]) ntrue++;

/*
  nfalse=0;
  corr=0;
  acc_mat=new long *[nclass];
  acc_mat[0]=new long[nclass*nclass];
  for (long i=1; i<nclass; i++) acc_mat[i]=&acc_mat[0][i*nclass];
  for (long i=0; i<nclass*nclass; i++) acc_mat[0][i]=0;

  printf("  number cls ret   con nfalse  ntrue\n");
  for (long i=0; i<ntest; i++) {
    if (cls[i] == result[i]) ntrue++; else nfalse++;
    //printf("%8d %3d %3d %6d %6d\n", i, cls[i], result[i], nfalse, ntrue);
  }
*/

  //evaluate results:
  class_eval(cls, result, ntest);
  check_confidence(cls, result, con, nclass, NCONHIST);

  //clean up:

  //delete the temporary files so they don't clutter the directory:
  sprintf(command, "rm %s", testfile);
  printf("%s\n", command);
  system(command);
  sprintf(command, "rm %s.cls", testbase);
  printf("%s\n", command);
  system(command);
  sprintf(command, "rm %s.con", testbase);
  printf("%s\n", command);
  system(command);
  sprintf(command, "rm %s.vec", trainbase);
  printf("%s\n", command);
  system(command);
  sprintf(command, "rm %s.cls", trainbase);
  printf("%s\n", command);
  system(command);
  sprintf(command, "rm %s.brd", brdfile);
  printf("%s\n", command);
  system(command);
  sprintf(command, "rm %s.bgd", brdfile);
  printf("%s\n", command);
  system(command);

  delete [] result;
  delete [] all;
  delete [] train;
  delete [] cls;
  
  //delete various file names:
  delete [] trainbase;
  delete [] trainvecfile;
  delete [] trainclsfile;
  
  delete [] testbase;
  delete [] testfile;
  delete [] outfile;
  delete [] confile;
  
  delete [] brdfile;
  
//  delete [] clind;
  //delete [] acc_mat[0];
  //delete [] acc_mat;

}


