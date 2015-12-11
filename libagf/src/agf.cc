//usage: agf [options] mode model test output
//mode = (pdf|classify|interpolate) (agf|knn)
#include <string.h>
#include <assert.h>

#include <gsl/gsl_linalg.h>

#include "linked.h"
#include "full_util.h"
#include "agf_lib.h"
#include "randomize.h"

using namespace libagf;
using namespace libpetey;

int main(int argc, char **argv) {
  FILE *fs;
  FILE *logfs=stderr;	//log messages
  char *vecfile;
  char *ordfile;
  char *resultfile;
  char *confile;
  char *command;	//command for pre-processing
  int err=0;		//return error code
  real_a **train;	//training features
  void *ord;		//ordinates
  dim_ta nvar;		//number of features
  nel_ta ntrain;	//number of training samples
  nel_ta n1;
  cls_ta nclass;	//number of classes
  real_a **test;	//test data
  nel_ta ntest;		//number of test points

  void *result;		//results
  real_a *con;		//confidence/error
  dim_ta nvar1;
  real_a *pdf=NULL;	//cond. prob.
  real_a pcor;		//correction for pdf calcs for norm. coord
  size_t ressize;	//size of one data element in results array

  int action;	//0=classify, 1=interpolation, 2=pdf estimation

  agf_command_opts opt_args;

  char pdformat[10];

  //diagnostic shit:
  agf_diag_param diag_param;
  iter_ta min_nd, max_nd, total_nd;
  real_a min_f, max_f, total_f;
  real_a min_W, max_W, total_W;

  opt_args.W2=W_DEFAULT;		//why is it W2?
  opt_args.k=K_DEFAULT_AGF;

  err=agf_parse_command_opts(argc, argv, "a:I:k:l:N:S:W:v:V:nj", &opt_args);
  if (err==FATAL_COMMAND_OPTION_PARSE_ERROR) exit(err);

  if (argc < 4) {
    printf("\n");
    printf("Syntax:      agf [-n] [-j] [-W Wc] [-v var1] [-V var2] [-k k] [..] \\\n");
    printf("                   action train test output\n");
    printf("\n");
    printf("arguments:\n");
    printf("    action   action to perform:\n");
    printf("               classify = statistical classification\n");
    printf("               interp   = interpolation/regression\n");
    printf("               pdf      = estimate probability densities\n");
    printf("    train    files containing training data:\n");
    printf("               .vec for vectors\n");
    printf("               .cls for classes\n");
    printf("               .dat for floating point ordinates\n");
    printf("    test     file containing vector data for which estimates are desired\n");
    printf("    output   files containing the results of the estimation:\n");
    printf("               .cls for classification results\n");
    printf("               .dat for floating point interpolates/pdf estimates\n");
    printf("               .con for classification confidence ratings\n");
    printf("               .err for interpolation error estimates\n");
    printf("\n");
    printf("options:\n");
    printf("    -a normfile     file containing normalization/transformation data\n");
    printf("    -I/-N maxiter   maximum number of iterations in supernewton (%d)\n", WEIGHTS_MAXITER);
    printf("    -k k     number of nearest neighbours to use in each estimate\n");
    printf("               --default is to use all of the data\n");
    printf("    -l tol   tolerance of W (default=%g)\n", WEIGHTS_TOL);
    printf("    -S nsv   perform SVD, keep nsv singular values\n");
    printf("    -W Wc    objective total weight (default=%6.1f)\n", opt_args.W2);
    printf("    -v var1  first bracket of filter variance\n");
    printf("               --default is to use the total variance/n^(2/D)\n");
    printf("    -V var2  second bracket of filter variance/initial filter variance\n");
    printf("               --default is to use the total variance of the data\n\n");
    printf("flags:\n");
    printf("    -j       print joint instead of cond. prob. to stdout\n");
    printf("    -n       normalize the data\n");
    printf("\n");
    exit(0);
  }

  ran_init();		//need random numbers to resolve ties

  if (opt_args.jointflag) strcpy(pdformat, "%12.6g "); else strcpy(pdformat, "%8.6f ");
 
  //start right of by determining what action to perform: 
  if (strcmp(argv[0], "classify")==0) {
    action=0;
  } else if (strcmp(argv[0], "interp")==0) {
    action=1;
  } else if (strcmp(argv[0], "pdf")==0) {
    action=2;
  } else {
    fprintf(stderr, "agf: action, '%s', not recognized\n", argv[0]);
    exit(PARAMETER_OUT_OF_RANGE);
  }

  vecfile=new char[strlen(argv[1])+5];
  sprintf(vecfile, "%s.vec", argv[1]);

  //if we need a normalization file and one hasn't been named, 
  //construct the name:
  if ((opt_args.svd>0 || opt_args.normflag) && opt_args.normfile == NULL) {
    opt_args.normfile=new char[strlen(argv[3])+5];
    sprintf(opt_args.normfile, "%s.std", argv[3]);
  }

  //get the training co-ordinate data, pre-process if necessary
  train=agf_get_features(argv[1], &opt_args, ntrain, nvar);

  fprintf(logfs, "%d training vectors found: %s.vec\n", ntrain, argv[1]);

  //dammit: recalculate the variances, even if they've just been calculated, above??
  if (opt_args.var[0] <= 0 || opt_args.var[1] <= 0) {
    //calculate the averages and standard deviations:
    real_a std[nvar];
    real_a ave[nvar];
    real_a vart;

    calc_norm(train, nvar, ntrain, ave, std);

    if (opt_args.normflag==0) {
      fprintf(logfs, "Statistics:\n");
      print_stats(logfs, ave, std, nvar);
      fprintf(logfs, "\n");
    }

    //if the initial filter variance is not set, set it to the total
    //variance of the data:
    vart=0;
    for (dim_ta i=0; i<nvar; i++) vart+=std[i]*std[i];
    if (opt_args.var[0] <= 0) {
      opt_args.var[0]=vart/pow(ntrain, 2./nvar);
      fprintf(logfs, "Using %10.3g for lower filter variance bracket\n\n", opt_args.var[0]);
    }
    if (opt_args.var[1] <= 0) {
      opt_args.var[1]=vart;
      fprintf(logfs, "Using %10.3g for upper filter variance bracket\n\n", opt_args.var[1]);
    }
  }

  //check the range of k:
  if (opt_args.k <= opt_args.W2 || opt_args.k >= ntrain) {
    if (opt_args.k != -1) {
      fprintf(stderr, "agf: Parameter k=%d out of range.  Using all the training data.\n", opt_args.k);
      opt_args.k=-1;
      err=PARAMETER_OUT_OF_RANGE;
    }
  }

  //read test data:
  test=read_vecfile(argv[2], ntest, nvar1);
  if (nvar1 == -1) {
    fprintf(stderr, "Error reading input file: %s\n", argv[2]);
    return FILE_READ_ERROR;
  }
  if (ntest == -1) {
    fprintf(stderr, "Error reading input file: %s\n", argv[2]);
    return ALLOCATION_FAILURE;
  }
  if (test == NULL) {
    fprintf(stderr, "Unable to open file for reading: %s\n", argv[2]);
    return UNABLE_TO_OPEN_FILE_FOR_READING;
  }

  //normalize test data:
  real_a **mat;
  if (opt_args.normfile!=NULL) {
    real_a *b;
    dim_ta nvar2, nvar3;
    real_a **testnew;

    mat=read_stats2(opt_args.normfile, b, nvar2, nvar3);
    printf("nvar=%d; nvar3=%d\n", nvar, nvar3);
    assert(nvar3==nvar);			//should be true
    if (nvar2!=nvar1) {
      fprintf(stderr, "agf: incorrect number of dimensions in test data (%d found, %d expected)\n", nvar1, nvar2);
      exit(DIMENSION_MISMATCH);
    }
    for (nel_ta i=0; i<ntest; i++) {
      for (dim_ta j=0; j<nvar2; j++) {
        test[i][j]-=b[j];
      }
    }
    testnew=matrix_mult(test, mat, ntest, nvar1, nvar);
    delete_matrix(test);
    test=testnew;

    delete [] b;
    //delete_matrix(mat);
  }
   
  fprintf(logfs, "%d test vectors found in file %s\n", ntest, argv[1]);

  ordfile=new char[strlen(argv[1])+5];
  resultfile=new char[strlen(argv[3])+5];
  confile=new char[strlen(argv[3])+5];

  //initialize diagnostic values:
  min_nd=agf_global_weights_maxiter+2;
  max_nd=0;
  total_nd=0;

  min_f=1;
  max_f=0;
  total_f=0;

  min_W=1000*opt_args.W2;
  max_W=0;
  total_W=0;

  //action specific file stuff:
  switch (action) {
    case 0: 
      sprintf(resultfile, "%s.cls", argv[3]);
      sprintf(confile, "%s.con", argv[3]);
      sprintf(ordfile, "%s.cls", argv[1]);
      ord=read_clsfile(ordfile, n1);
      ressize=sizeof(cls_ta);
      break;
    case 1:
      sprintf(resultfile, "%s.dat", argv[3]);
      sprintf(confile, "%s.err", argv[3]);
      sprintf(ordfile, "%s.dat", argv[1]);
      ord=read_datfile(ordfile, n1);
      ressize=sizeof(real_a);
      break;
    case 2:
      sprintf(resultfile, "%s.dat", argv[3]);
      ord=NULL;
      ressize=sizeof(real_a);
      //calculate correction value for normalized coords:
      pcor=1;
      if (opt_args.normfile!=NULL) {
        gsl_matrix *V;
        gsl_matrix *U;
        gsl_vector *S;
        gsl_vector *work;
        //compute the singular value decomposition to get the determinant:
        U=gsl_matrix_alloc(nvar1, nvar);
        for (dim_ta i=0; i<nvar1; i++) {
          for (dim_ta j=0; j<nvar; j++) gsl_matrix_set(U, i, j, mat[i][j]);
        }
        V=gsl_matrix_alloc(nvar, nvar);
        S=gsl_vector_alloc(nvar);
        work=gsl_vector_alloc(nvar);
        gsl_linalg_SV_decomp(U, V, S, work);
        for (dim_ta i=0; i<nvar; i++) {
          real_a s_i=gsl_vector_get(S, i);
          if (s_i>0) pcor=pcor*s_i;
        }
        gsl_matrix_free(U);
        gsl_matrix_free(V);
        gsl_vector_free(S);
        gsl_vector_free(work);
      }
      break;
  }

  if (action!=2) {
    if (n1 == -1) {
      fprintf(stderr, "Error reading file: %s\n", ordfile);
      exit(FILE_READ_ERROR);
    }
    if (ord == NULL && action!=2) {
      fprintf(stderr, "Unable to open file, %s, for reading.\n", ordfile);
      exit(UNABLE_TO_OPEN_FILE_FOR_READING);
    }
    if (n1!=ntrain) {
      fprintf(stderr, "Sample count mismatch: %d in %s, %d in %s.\n", ntrain, vecfile, n1, ordfile);
      exit(SAMPLE_COUNT_MISMATCH);
    }
  }

  if (action==0) {
    nclass=1;
    for (nel_ta i=0; i<n1; i++) if (((cls_ta *) ord)[i]>=nclass) nclass=((cls_ta *) ord)[i]+1;
    pdf=new real_a[nclass];
  }

  if (action!=2) con=new real_a[ntest]; else con=NULL;

  result=malloc(ressize*ntest);

  //don't read on...  this is going to be a fucking disaster...
  for (nel_ta i=0; i<ntest; i++) {
    real_a p_x;				//P(x)
    if (opt_args.k<0) {
      switch (action) {
        case 0:
          ((cls_ta *) result)[i]=agf_classify(train, nvar, (cls_ta *) ord, ntrain, nclass, 
			test[i], opt_args.var, opt_args.W2, pdf, &diag_param, opt_args.jointflag);
          if (opt_args.jointflag) {
            p_x=0;
            for (cls_ta j=0; j<nclass; j++) p_x+=pdf[j];
            con[i]=(nclass*pdf[((cls_ta *)result)[i]]/p_x-1)/(nclass-1);
          } else {
            con[i]=(nclass*pdf[((cls_ta *)result)[i]]-1)/(nclass-1);
          }
  
          //print results to standard out:
          for (cls_ta j=0; j<nclass; j++) printf(pdformat, pdf[j]);
          printf(" %4d", ((cls_ta *)result)[i]);
          printf("\n");
          break;
        case 1:
          ((real_a *) result)[i]=adgaf_err(train, nvar, (real_a *) ord, ntrain, test[i],
			opt_args.var, opt_args.W2, con[i], &diag_param);
          break;
        case 2:
          ((real_a *) result)[i]=agf_calc_pdf(train, nvar, ntrain, test[i], 
			opt_args.var, opt_args.W2, &diag_param)/pcor;
          break;
      }
    } else {
      switch (action) {
        case 0: 
          ((cls_ta *) result)[i]=agf_classify(train, nvar, (cls_ta *) ord, ntrain, nclass, 
			test[i], opt_args.var, opt_args.k, opt_args.W2, 
			pdf, &diag_param, opt_args.jointflag);
          if (opt_args.jointflag) {
            p_x=0;
            for (cls_ta j=0; j<nclass; j++) p_x+=pdf[j];
            con[i]=(nclass*pdf[((cls_ta *)result)[i]]/p_x-1)/(nclass-1);
          } else {
            con[i]=(nclass*pdf[((cls_ta *)result)[i]]-1)/(nclass-1);
          }
  
          //print results to standard out: (use lvq-compatible format)
          for (cls_ta j=0; j<nclass; j++) printf(pdformat, pdf[j]);
          printf(" %4d", ((cls_ta *)result)[i]);
          printf("\n");
          break;
        case 1:
          ((real_a *) result)[i]=adgaf_err(train, nvar, (real_a *) ord, ntrain, test[i],
			opt_args.var, opt_args.k, opt_args.W2, con[i], &diag_param);
          break;
        case 2:
          ((real_a *) result)[i]=agf_calc_pdf(train, nvar, ntrain, test[i], 
			opt_args.var, opt_args.k, opt_args.W2, &diag_param)/pcor;
          break;
      }
      if (diag_param.f < min_f) min_f=diag_param.f;
      if (diag_param.f > max_f) max_f=diag_param.f;
      total_f+=diag_param.f;
    }
    //calculate diagnostics:
    if (diag_param.nd < min_nd) min_nd=diag_param.nd;
    if (diag_param.nd > max_nd) max_nd=diag_param.nd;
    total_nd+=diag_param.nd;

    if (diag_param.W < min_W) min_W=diag_param.W;
    if (diag_param.W > max_W) max_W=diag_param.W;
    total_W+=diag_param.W;
  }

  //print out diagnostics:
  FILE *diagfs=logfs;
  fprintf(diagfs, "\n");
  fprintf(diagfs, "diagnostic parameter          %8s   %8s   %8s\n",
                      "min", "max", "average");
  fprintf(diagfs, "iterations in agf_calc_w:      %8d   %8d %10.3g\n",
                   min_nd, max_nd, (real_a) total_nd/(real_a) ntest);
  if (opt_args.k != -1) {
    fprintf(diagfs, "value of f:                    %10.3g %10.3g %10.3g\n",
                   min_f, max_f, total_f/ntest);
  }
  fprintf(diagfs, "value of W:                    %10.2f %10.2f %10.2f\n",
                   min_W, max_W, total_W/ntest);
  fprintf(diagfs, "\n");

  //write the results to a file:
  fs=fopen(resultfile, "w");
  if (fs == NULL) {
    fprintf(fs, "Unable to open file for writing: %s\n", resultfile);
    return UNABLE_TO_OPEN_FILE_FOR_WRITING;
  }
  fwrite(result, ressize, ntest, fs);
  fclose(fs);

  if (action!=2) {
    fs=fopen(confile, "w");
    if (fs == NULL) {
      fprintf(fs, "Unable to open file for writing: %s\n", confile);
      return UNABLE_TO_OPEN_FILE_FOR_WRITING;
    }
    fwrite(con, sizeof(real_a), ntest, fs);
    fclose(fs);
  }

  delete_matrix(train);
  if (action==0) delete [] (cls_ta *) ord;
  if (action==1) delete [] (real_a *) ord;
  delete_matrix(test);
  if (con!=NULL) delete [] con;
  if (pdf!=NULL) delete [] pdf;
  free(result);

  delete [] vecfile;
  delete [] ordfile;
  delete [] resultfile;
  delete [] confile;

  if (opt_args.normfile!=NULL) {
    delete_matrix(mat);
    delete [] opt_args.normfile;
  }

  ran_end();

  return err;

}
