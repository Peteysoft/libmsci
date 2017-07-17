//
// This software is released under the following terms:
//
// 1. No commercial use.
// 2. Copies and derivative works must carry a copy of this license.
// 3. Attribution must be given to all contributors of both original and derivative works.
//
// Authors:
//
// 2017-07-16 Peter Mills: added license information 
//

//
// What the fuck actually does this thing do? Whatever it is, I don't think it
// works....
//

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>

#include "full_util.h"
#include "linked.h"
#include "agf_lib.h"

using namespace std;
using namespace libagf;
using namespace libpetey;

int main(int argc, char *argv[]) {
  char *vecfile;	//back to the old way...
  char *classfile=NULL;
  char *modelfile;	//optional model data to pass to external command
  char *outfile;	//output file
  char *cpcommand;	//copy input borders files to output
  int cpnormfile=0;	//copy the normalization file?

  agf_command_opts opt_args;

  real_a **train;	//training vectors
  cls_ta *cls;		//training data classes

  nel_ta ntrain;		//number of training data points
  dim_ta nvar=0;	//number of variables
  nel_ta n1;		//number of class labels

  borders_calibrated<real_a, cls_ta> *classifier;

  FILE *fs;		//file stream

  cls_ta nclass;		//number of classes (anything above two is ignored
  			//though...)

  real_a *ave=NULL;
  real_a **mat=NULL;
  real_a *all;
  dim_ta nvar1, nvar2;	//check number of variables against other files

  real_a r0;		//threshold probability

  int exit_code=0;
  int err_code;

  //set defaults:
  opt_args.nt=NCONHIST;
  opt_args.Qtype=CALIBRATION_ORDER;
  opt_args.algtype=0;
  opt_args.rthresh=0;

  //normalization options: -n -S -a
  //supernewton iteration: -h -i -I
  exit_code=agf_parse_command_opts(argc, argv, "a:c:q:Q:r:nuAM", &opt_args);
  if (exit_code==FATAL_COMMAND_OPTION_PARSE_ERROR) return exit_code;

  if (argc < 3) {
    printf("\n");
    printf("purpose:  calibrates the probability estimates of a borders classifier\n");
    printf("\n");
    printf("syntax:   calibrate_borders [-n] [-u] [-a normfile] [-q nhist] [-Q order] \\\n");
    printf("                   [-A [-M]]  border-in train border-out \\\n");
    printf("                   [cls1a cls1b cls1c ... %c cls2a cls2b cls2c ...] \n", PARTITION_SYMBOL);
    printf("\n");
    printf("arguments:\n");
    printf("  border-in   base name of input border classifier:\n");
    printf("                .brd samples the border;\n");
    printf("                .bgd contains gradient vectors\n");
    printf("                .std contains normalization data (unless -a specified)\n");
    printf("                .clb calibration file\n");
    printf("  train       binary files containing locations of the samples:\n");
    printf("                .vec for vectors;\n");
    printf("                .cls for classes\n");
    printf("  border-out  base name of output border classifier\n");
    printf("  clsIJ       for partitioning multiple classes: Jth member of the Ith partition\n");
    printf("\n");
    printf("options:\n");
    printf("  -a normfile file containing normalization data (input/output)\n");
    printf("  -A          input file in ASCII format\n");
    printf("  -c          method: [0]\n");
    printf("                0 = free fit\n");
    printf("                1 = fix threshold (r0) based on -r\n");
    printf("                2 = optimize accuracy\n");
    printf("                3 = optimize uncertainty coefficient\n");
    printf("                4 = optimize correlation coefficient\n");
    printf("  -M          input file in LIBSVM format\n");
    printf("  -q nhist    number of histogram partitions [%d]\n", NCONHIST);
    printf("  -Q order    order of fitted polynomial [%d]\n", CALIBRATION_ORDER);
    printf("  -n          option to normalise the data\n");
    printf("  -r r0       threshold probability\n");
    printf("  -u          store borders data in un-normalized coordinates\n");
    printf("\n");
    printf("*** to calibrated multi-borders classifieres, use the multi_borders command\n");
    printf("    to recursively calibrate each of the binary models, e.g.:\n");
    printf("\n");
    printf("$> multi_borders -Z -- \"calibrate_borders -c 3\" control train fbase output\n");
    printf("\n");
    return INSUFFICIENT_COMMAND_ARGS;
  }

  if (opt_args.asciiflag) {
    vecfile=new char[strlen(argv[1])+1];
    sprintf(vecfile, "%s", argv[1]);
  } else {
    vecfile=new char[strlen(argv[1])+5];
    sprintf(vecfile, "%s.vec", argv[1]);
  }

  //get the training co-ordinate data:
  fs=fopen(vecfile, "r");
  if (fs==NULL) {
    fprintf(stderr, "class_borders: unable to open file, %s, for reading\n", vecfile);
    exit(UNABLE_TO_OPEN_FILE_FOR_READING);
  }
  if (opt_args.asciiflag) {
    if (opt_args.Mflag) {
      ntrain=read_svm(fs, train, cls, nvar, opt_args.missing, opt_args.Uflag);
    } else {
      ntrain=read_lvq(fs, train, cls, nvar, opt_args.Hflag);
    }
  } else {
    nel_ta nvar1;
    train=read_matrix<real_a, nel_ta>(fs, ntrain, nvar1);
    nvar=nvar1;
  }
  if (nvar==-1 || ntrain==-1) {
    fprintf(stderr, "class_borders: Error reading input file, %s\n", vecfile);
    exit(FILE_READ_ERROR);
  }
  fclose(fs);

  //for binary data, read in class data separately:
  if (opt_args.asciiflag==0) {
    classfile=new char[strlen(argv[1])+5];
    sprintf(classfile, "%s.cls", argv[1]);
    cls=read_clsfile<cls_ta>(classfile, n1);
    if (n1 == -1) {
      fprintf(stderr, "class_borders: Error reading data file, %s\n", classfile);
      exit(FILE_READ_ERROR);
    }
    if (cls == NULL) {
      fprintf(stderr, "Unable to open file, %s, for reading.\n", classfile);
      exit(UNABLE_TO_OPEN_FILE_FOR_READING);
    }
    if (n1!=ntrain) {
      fprintf(stderr, "Sample count mismatch: %d in %s, %d in %s.\n", ntrain, vecfile, n1, classfile);
      exit(SAMPLE_COUNT_MISMATCH);
    }
    delete [] classfile;
  }

  printf("%d %d-dimensional training vectors found: %s\n", ntrain, nvar, argv[0]);

  //count the number of classes:
  nclass=1;
  for (nel_ta i=0; i<ntrain; i++) if (cls[i]>=nclass) nclass=cls[i]+1;

  //if there are partitions:
  nel_ta *clind;
  if (argc>3) {
    cls_ta map[nclass];
    cls_ta nncls=1;			//new number of classes
    err_code=parse_partition(argc-3, argv+3, nclass, map);
    if (err_code!=0) {
      fprintf(stderr, "class_borders: error parsing class partition\n");
      exit(err_code);
    }
    apply_partition(cls, ntrain, map);
    for (cls_ta i=0; i<nclass; i++) if (map[i]>=nncls) nncls=map[i]+1;
    nclass=nncls;
  }

  //remove excluded/out-of-range classes:
  all=train[0];
  clind=sort_classes(train, ntrain, cls, nclass);

  if (nclass < 2) {
    fprintf(stderr, "class_borders: Cannot perform classifications with less than two classes!\n");
    return PARAMETER_OUT_OF_RANGE;
  }

  classifier=new borders_calibrated<real_a, cls_ta>(argv[0]);

  //normalization:
  if ((opt_args.uflag || opt_args.normflag) && opt_args.normfile==NULL) {
    cpnormfile=1;
    opt_args.normfile=new char [strlen(argv[0])+5];
    sprintf(opt_args.normfile, "%s.std", argv[0]);
  }

  if (opt_args.normfile!=NULL) {
    mat=read_stats2(opt_args.normfile, ave, nvar1, nvar2);

    err_code=classifier->ltran(mat, ave, nvar1, nvar2, opt_args.uflag);
    if (err_code!=0) exit(err_code);

    if (classifier->n_feat_t() != nvar) {
      fprintf(stderr, "classify_b: Dimensions of classifier (%d) do not match those of test data (%d).\n",
                classifier->n_feat_t(), nvar);
      exit(DIMENSION_MISMATCH);
    }
  } else {
    if (classifier->n_feat_t() != nvar) {
      fprintf(stderr, "classify_b: Dimension of classifier (%d) do not match dimension of test data (%d).\n",
                classifier->n_feat_t(), nvar);
      exit(DIMENSION_MISMATCH);
    }
  }

  classifier->calibrate2(train+clind[0], cls+clind[0], clind[2]-clind[0], opt_args.Qtype, opt_args.nt, opt_args.algtype, opt_args.rthresh);

  //create the new calibration file:
  outfile=new char[strlen(argv[2])+5];
  sprintf(outfile, "%s.clb", argv[2]);
  fs=fopen(outfile, "w");
  classifier->print_calib(fs);
  fclose(fs);

  //copy all the other stuff over:
  cpcommand=new char [strlen(argv[0])+strlen(argv[2])+13];
  if (cpnormfile) {
    sprintf(cpcommand, "cp %s.std %s.std", argv[0], argv[2]);
    err_code=system(cpcommand);
    printf("%s\n", cpcommand);
  }
  sprintf(cpcommand, "cp %s.brd %s.brd", argv[0], argv[2]);
  printf("%s\n", cpcommand);
  err_code=system(cpcommand);
  sprintf(cpcommand, "cp %s.bgd %s.bgd", argv[0], argv[2]);
  printf("%s\n", cpcommand);
  err_code=system(cpcommand);

  exit(0);

  //delete character strings containing file names:
  delete [] vecfile;
  delete [] outfile;
  delete [] cpcommand;

  //delete integer and real_aing point arrays:
  delete[] all;
  delete[] train;
  delete [] cls;
  delete [] clind;

  if (ave!=NULL) delete [] ave;
  if (mat!=NULL) delete_matrix(mat);

  if (opt_args.normfile!=NULL) {
    delete [] opt_args.normfile;
  }

  return 0;

}


