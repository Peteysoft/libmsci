
//Copyright (C) 2007 Peter Mills.  All rights reserved.

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
  char *brdfile;	//binary file sampling border
  char *grdfile;	//binary file containing gradient vectors

  agf_command_opts opt_args;

  real_a **train;	//training vectors
  cls_ta *cls;		//training data classes

  nel_ta ntrain;		//number of training data points
  nel_ta ntrain2;	//revised training samples (after excluding -1 classes)
  dim_ta nvar=0;	//number of variables
  nel_ta n1;		//number of class labels

  real_a **border;	//border vectors
  real_a **gradient;	//gradient vectors

  FILE *fs;		//file stream
  char *command;	//pre-process command

  cls_ta nclass;		//number of classes (anything above two is ignored
  			//though...)
  nel_ta *clind;		//indices of sorted classes
  nel_ta clind2;	//for the two-class classification 
  			//(relative to start of non-excluded data)

  real_a *all;		//for deleting the training vectors

  dim_ta nvar1, nvar2;	//check number of variables against other files
  real_a *std, *ave;	//average, standard deviation...
  real_a vart;

  int exit_code;
  int err_code;

  //set defaults and parse command line options:
  opt_args.k=-1;
  opt_args.W2=W_DEFAULT_BORDERS;

  //AFG options: -v -V -W -k
  //borders options: -s -t -r -T
  //normalization options: -n -S -a
  //supernewton iteration: -h -i -I
  exit_code=agf_parse_command_opts(argc, argv, "a:h:i:I:k:l:N:q:r:s:t:T:v:V:W:nuS:O:KMU", &opt_args);
  if (exit_code==FATAL_COMMAND_OPTION_PARSE_ERROR) return exit_code;

  if (argc < 2) {
    printf("\n");
    printf("syntax:  class_borders [-n] [-s nsv] [-u] [-a normfile] \\\n");
    printf("                   [-k k] [-W Wc] [-s n] [-t tol] [-r r0] [-v var1] [-v var2] \\\n");
    printf("                   [-O command [-K] [-M] model]  train border \\\n");
    printf("                   [cls1a cls1b cls1c ... %c cls2a cls2b cls2c ...] \n", PARTITION_SYMBOL);
    printf("\n");
    printf("arguments:\n");
    printf("  model       optional argument used in conjunction with -O option.  See below.\n");
    printf("  train       binary files containing locations of the samples:\n");
    printf("                .vec for vectors;\n");
    printf("                .cls for classes\n");
    printf("  border      base name of output files:\n");
    printf("                .brd samples the border;\n");
    printf("                .bgd contains gradient vectors\n");
    printf("                .std contains normalization data (unless -a specified)\n");
    printf("  clsIJ       for partitioning multiple classes: Jth member of the Ith partition\n");
    printf("\n");
    printf("options:\n");
    printf("  -a normfile file containing normalization data (input/output)\n");
    printf("  -h hrel     relative displacement in numerical derivative calculations [%g]\n", opt_args.hrel);
    printf("  -i maxit1   maximum number of iterations when searching for class border (%d)\n", (int32_t) agf_global_borders_maxiter);
    printf("  -I maxit2   maximum number of iterations when calculating weights (%d)\n", (int32_t) agf_global_weights_maxiter);
    printf("  -k k        number of nearest neighbours to use in each calculation\n");
    printf("                --default is to use all the data\n");
    printf("  -K          keep temporary files\n");
    printf("  -l tol      tolerance of W (default=%g)\n", (float) agf_global_weights_tol);
    printf("  -M          files in LIBSVM format\n");
    printf("  -O command  command used to return probability estimates from model\n");
    printf("                (if specified, training data is assumed to be in ASCII format)\n");
    printf("  -n          option to normalise the data\n");
    printf("  -N maxit3   maximum number of iterations in supernewton (%d, %d)\n", (int32_t) agf_global_weights_maxiter, (int32_t) agf_global_borders_maxiter);
    printf("  -q niter    number of iterations in numerical derivative calculations [%d]\n", opt_args.nt);
    printf("  -r r0       location of discrimination border (default=0)\n");
    printf("  -s n        number of times to sample the border (default=%d)\n", (int32_t) opt_args.n);
    printf("  -S nsv      perform SVD, keep nsv singular values\n");
    printf("  -t tol      tolerance of border samples (default=%g)\n", (float) opt_args.tol);
    printf("  -T cthresh  class threshold (default=%d)\n", 1);
    printf("  -u          store borders data in un-normalized coordinates\n");
    printf("  -v var1     lower filter variance bracket\n");
    printf("                --default is to use (the total variance of the data)/n^(2/D)\n");
    printf("  -V var2     upper filter variance bracket\n");
    printf("                --default is to use the total variance of the data\n");
    printf("  -W Wc       objective total weight (default=%g)\n", (float) opt_args.W2);
    printf("\n");
    printf("  To train the class borders using an external command, set the -O option to\n");
    printf("to supply a command that, in conjunction with the model parameters, returns\n");
    printf("probability estimates for use in the root-finding algorithm.  The train\n");
    printf("parameter is still required for drawing samples on either side of the border,\n");
    printf("but need not be the same data as used to train the original model.\n");
    printf("\n");
    printf("The external command is called exactly as svm-predict (part of LIBSVM) as\n");
    printf("follows:\n");
    printf("\n");
    printf("  command test model output\n");
    printf("\n");
    printf("where:\n");
    printf("- command is the name of the command as supplied by the -O option, \n");
    printf("- test is a file containing test data generated by the root-finding algorithm\n");
    printf("- model is the name of the model file as supplied in the input parameters\n");
    printf("- output is a file containing output classes and probabilities.\n");
    printf("\n");
    printf("The output file format is the same as svm-predict, while the input format in\n");
    printf("test is the LVQPAK ASCII format unless the -M option is specified in which case\n");
    printf("LIBSVM format is used.  File formats for the train parameter switches to ASCII\n");
    printf("as well.  Output file format for the final, returned border parameter stays the same.\n");
    return INSUFFICIENT_COMMAND_ARGS;
  }

  ave=NULL;
  std=NULL;

  if (opt_args.multicommand!=NULL) {
    opt_args.asciiflag=1;
    modelfile=argv[0];
    argv++;
    argc--;
  }

  if (opt_args.asciiflag) {
    vecfile=new char[strlen(argv[0])+1];
    sprintf(vecfile, "%s", argv[0]);
  } else {
    vecfile=new char[strlen(argv[0])+5];
    sprintf(vecfile, "%s.vec", argv[0]);
  }

  //if we need a normalization file and one hasn't been named, 
  //construct the name:
  if ((opt_args.svd>0 || opt_args.normflag) && opt_args.normfile == NULL) {
    opt_args.normfile=new char[strlen(argv[1])+5];
    sprintf(opt_args.normfile, "%s.std", argv[1]);
  }

  //get the training co-ordinate data, pre-process if necessary
  command=compile_precommand(vecfile, &opt_args);
  if (command==NULL) {
    fs=fopen(vecfile, "r");
  } else {
    fprintf(stderr, "%s\n", command);
    fs=popen(command, "r");
  }
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
  if (command==NULL) {
    fclose(fs);
  } else {
    pclose(fs);
    delete [] command;
  }

  //for binary data, read in class data separately:
  if (opt_args.asciiflag==0) {
    classfile=new char[strlen(argv[0])+5];
    sprintf(classfile, "%s.cls", argv[0]);
    cls=read_clsfile(classfile, n1);
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
  }

  printf("%d %d-dimensional training vectors found: %s\n", ntrain, nvar, argv[0]);
  all=train[0];		//this will help us down the line: since sub-pointers are re-arranged
			//deleting the 0th element may cause an error...

  //count the number of classes:
  nclass=1;
  for (nel_ta i=0; i<ntrain; i++) if (cls[i]>=nclass) nclass=cls[i]+1;
  if (nclass < 2) {
    fprintf(stderr, "class_borders: Cannot perform classifications with less than two classes!\n");
    return PARAMETER_OUT_OF_RANGE;
  }
  if (opt_args.cl_thresh >= nclass) {
    fprintf(stderr, "class_borders: Class threshold greater than number of classes.\n");
    return PARAMETER_OUT_OF_RANGE;
  }

  //if there are partitions:
  if (argc>2) {
    cls_ta map[nclass];
    cls_ta nncls=1;			//new number of classes
    err_code=parse_partition(argc-2, argv+2, nclass, map);
    if (err_code!=0) {
      fprintf(stderr, "class_borders: error parsing class partition\n");
      exit(err_code);
    }
    apply_partition(cls, ntrain, map);
    for (cls_ta i=0; i<nclass; i++) if (map[i]>=nncls) nncls=map[i]+1;
    nclass=nncls;
  }

  brdfile=new char[strlen(argv[1])+5];
  strcpy(brdfile, argv[1]);
  strcat(brdfile, ".brd");

  grdfile=new char[strlen(argv[1])+5];
  strcpy(grdfile, argv[1]);
  strcat(grdfile, ".bgd");

  //dammit: recalculate the variances, even if they've just been calculated, above??
  if (opt_args.var[0] <= 0 || opt_args.var[1] <= 0) {
    //calculate the averages and standard deviations:
    std=new real_a[nvar];
    ave=new real_a[nvar];
    calc_norm(train, nvar, ntrain, ave, std);

    printf("Statistics:\n");
    print_stats(stdout, ave, std, nvar);
    printf("\n");

    //if the initial filter variance is not set, set it to the total
    //variance of the data:
    vart=0;
    for (dim_ta i=0; i<nvar; i++) vart+=std[i]*std[i];
    if (opt_args.var[0] <= 0) {
      opt_args.var[0]=vart/pow(ntrain, 2./nvar);
      printf("Using %10.3g for lower filter variance bracket\n\n", opt_args.var[0]);
    }
    if (opt_args.var[1] <= 0) {
      opt_args.var[1]=vart;
      printf("Using %10.3g for upper filter variance bracket\n\n", opt_args.var[1]);
    }
    delete [] ave;
    delete [] std;
    ave=NULL;
    std=NULL;
  }

  //sort the classes:
  clind=sort_classes(train, ntrain, cls, nclass);
  //clind sorts all class labels larger than -1:
  ntrain2=clind[nclass]-clind[0];
  clind2=clind[opt_args.cl_thresh]-clind[0];

  //check the range of k:
  if (opt_args.k <= opt_args.W2 || opt_args.k >= ntrain2) {
    if (opt_args.k != -1) {
      fprintf(stderr, "class_borders: Parameter k=%d out of range.  Using all the training data.\n", opt_args.k);
      opt_args.k=-1;
      exit_code=PARAMETER_OUT_OF_RANGE;
    }
  }

  //allocate the arrays for holding the results:
  border=allocate_matrix<real_a, nel_ta>(opt_args.n, nvar);
  gradient=allocate_matrix<real_a, nel_ta>(opt_args.n, nvar);

  //find class borders:
  //if (ntrain < SMALL_DATASET) {
  //should really be based on the ratio between possible combinations
  //(in this case approximate) and the desired number:
  //printf("cb: test=%f\n", (2.*opt_args.n)/(ntrain-clind[1])/clind[1]);
  if (opt_args.multicommand==NULL) {
    if ((2.*opt_args.n)/(ntrain2-clind2)/clind2 > 0.25) {
      opt_args.n=find_class_borders_small(train+clind[0], nvar, ntrain2, clind2, opt_args.n,
		opt_args.var, opt_args.k, opt_args.W2, opt_args.tol,
		border, gradient, opt_args.rthresh);
    } else {
      opt_args.n=find_class_borders(train+clind[0], nvar, ntrain2, clind2, opt_args.n,
		opt_args.var, opt_args.k, opt_args.W2, opt_args.tol,
		border, gradient, opt_args.rthresh);
    }
  } else {
    bordparam<real_a> param;		//parameters for sampling function
    int (*sfunc) (void *, real_a *, real_a *);		//sampling function
    //initialize parameters:
    if ((2.*opt_args.n)/(ntrain2-clind2)/clind2 > 0.25) {
      //for small datasets:
      bordparam_init(&param, train+clind[0], nvar, ntrain2, clind2, 1);
      sfunc=&oppositesample_small<real_a>;
    } else {
      //for large datasets:
      bordparam_init(&param, train+clind[0], nvar, ntrain2, clind2);
      sfunc=&oppositesample<real_a>;
    }
    //not pretty, is it?
    opt_args.n=batch_borders<real_a, cls_ta>(opt_args.multicommand, modelfile, 
		sfunc, (void *) &param, 
		opt_args.n, nvar, opt_args.tol, agf_global_borders_maxiter, opt_args.hrel, 
		border, gradient, opt_args.rthresh, opt_args.Mflag, opt_args.Kflag,
		opt_args.nt);
    bordparam_clean(&param);
  }

  //un-normalize the vectors before writing them to a file:
  //(only if -u set)
  if (opt_args.uflag && opt_args.normfile!=NULL) {
    real_a **bord2;
    real_a **grad2;
    real_a **mat;
    gsl_matrix *mat2;
    gsl_vector *s;
    gsl_matrix *vt;
    gsl_vector *work;
   
    //is this efficient enough?
    mat=read_stats2(opt_args.normfile, ave, nvar1, nvar2);
    //too lazy to do this properly at the moment:
    if (mat==NULL || nvar1==-1 || nvar2==-1) {
      fprintf(stderr, "class_borders: error reading normalization file, %s\n", opt_args.normfile);
      exit(FILE_READ_ERROR);
    }
    if (nvar2!=nvar) {
      fprintf(stderr, "class_borders: second dimension of tranformation matrix (%d) does not match that of features data (%d)\n", nvar2, nvar);
      exit(DIMENSION_MISMATCH);
    }

    printf("Solving to convert borders data back to un-transformed coordinates...\n");
    if (nvar<nvar1) fprintf(stderr, "class_borders: warning, inverse coord. transformation may be under-determined\n");

    mat2=gsl_matrix_alloc(nvar1, nvar);
    s=gsl_vector_alloc(nvar);
    vt=gsl_matrix_alloc(nvar, nvar);
    work=gsl_vector_alloc(nvar);

    for (dim_ta i=0; i<nvar1; i++) {
      for (dim_ta j=0; j<nvar; j++) {
        gsl_matrix_set(mat2, i, j, mat[i][j]);
      }
    }
    //use SVD to do the inversion
    //since it will still work with transformation matrices
    //that reduce the dimension--if the dimension reduction has been done
    //effectively, the results might even be somewhat sensible
    gsl_linalg_SV_decomp(mat2, vt, s, work);
    gsl_vector_free(work);

    //custom implementation of SVD inversion:
    //(I think this is for efficiency...)
    bord2=allocate_matrix<real_a, nel_ta>(opt_args.n, nvar1);
    grad2=allocate_matrix<real_a, nel_ta>(opt_args.n, nvar1);
    for (nel_ta i=0; i<opt_args.n; i++) {
      double tmp_b;
      double tmp_g;
      for (dim_ta j=0; j<nvar1; j++) {
        bord2[i][j]=ave[j];
        grad2[i][j]=0;
        for (dim_ta k=0; k<nvar; k++) {
	  double mat2_el;
          double s_k=gsl_vector_get(s, k);
          if (s_k == 0) continue;
          tmp_b=0;
          tmp_g=0;
          for (dim_ta l=0; l<nvar; l++) {
	    double vt_el=gsl_matrix_get(vt, l, k);
            tmp_b+=vt_el*border[i][l];
            tmp_g+=vt_el*gradient[i][l];
          }
	  mat2_el=gsl_matrix_get(mat2, j, k);
          bord2[i][j]+=tmp_b*mat2_el/s_k;
          grad2[i][j]+=tmp_g*mat2_el/s_k;
        }
      }
    }
    nvar=nvar1;
    delete_matrix(border);
    delete_matrix(gradient);

    border=bord2;
    gradient=grad2;

    delete_matrix(mat);
    gsl_matrix_free(mat2);
    gsl_vector_free(s);
    gsl_matrix_free(vt);
  }
  
  //write them to a file:
  fs=fopen(brdfile, "w");
  fwrite(&nvar, sizeof(nvar), 1, fs);
  fwrite(border[0], sizeof(real_a), nvar*opt_args.n, fs);
  fclose(fs);

  fs=fopen(grdfile, "w");
  fwrite(&nvar, sizeof(nvar), 1, fs);
  fwrite(gradient[0], sizeof(real_a), nvar*opt_args.n, fs);
  fclose(fs);

  //clean up:

  //delete character strings containing file names:
  delete [] vecfile;
  if (classfile!=NULL) delete [] classfile;
  delete[] brdfile;
  delete[] grdfile;

  //delete integer and real_aing point arrays:
  delete[] train;
  delete[] all;
  if (ave != NULL) delete[] ave;
  delete[] clind;
  delete [] cls;

  delete_matrix(border);
  delete_matrix(gradient);

  if (opt_args.normfile!=NULL) {
    delete [] opt_args.normfile;
  }

  return exit_code;

}


