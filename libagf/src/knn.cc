//usage: agf [options] mode model test output
//mode = (pdf|classify|interpolate) (agf|knn)
#include <string.h>
#include <assert.h>

#include <gsl/gsl_linalg.h>

#include "linked.h"
#include "full_util.h"
#include "randomize.h"

#include "agf_lib.h"

using namespace std;
using namespace libpetey;
using namespace libagf;

int main(int argc, char **argv) {
  FILE *fs;
  FILE *logfs=stderr;		//informational messages
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
  real_a pcor;		//correction for pdf estimates
  size_t ressize;	//size of one data element in results array

  int action;	//0=classify, 1=interpolation, 2=pdf estimation

  char pdformat[10];

  agf_command_opts opt_args;

  opt_args.k=K_DEFAULT_KNN;

  err=agf_parse_command_opts(argc, argv, "a:k:m:S:nj", &opt_args);
  if (err==FATAL_COMMAND_OPTION_PARSE_ERROR) exit(err);

  //parse the command line arguments:
  if (argc < 4) {
    printf("\n");
    printf("purpose:      k-nearest neighbours machine learning algorithms\n");
    printf("\n");
    printf("syntax:       knn [-n] [-k k] [..] action train test output\n");
    printf("\n");
    printf("arguments:\n");
    printf("    action    action to perform\n");
    printf("                classify = statistical classification\n");
    printf("                interp   = interpolation/regression\n");
    printf("                pdf      = estimate probability densities\n");
    printf("    train     files containing training data:\n");
    printf("                .vec for vectors\n");
    printf("                .cls for classes\n");
    printf("                .dat for floating point ordinates\n");
    printf("    test      file containing vector data for which estimates are desired\n");
    printf("    output    files containing the results of the estimation:\n");
    printf("                .cls for classification results\n");
    printf("                .dat for floating point interpolates/pdf estimates\n");
    printf("                .con for classification confidence ratings\n");
    printf("                .err for interpolation error estimates\n");
    printf("\n");
    printf("options:\n");
    printf("    -a normfile file containing normalization/transformation data\n");
    printf("    -k k        number of nearest neighbours to use in each estimate\n");
    printf("                  (default=%d)\n", opt_args.k);
    printf("    -m metric   type of metric (classification and interpolation only)\n");
    printf("                  0=Cartesian (default), 1=Manhattan, 2=quartic, 3=max. element\n");
    printf("    -S nsv      perform SVD, keep nsv singular values\n");
    printf("\n");
    printf("flags:\n");
    printf("    -j          print joint instead of cond. prob. to stdout\n");
    printf("    -n          normalize the data\n");
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
  if ((opt_args.normflag || opt_args.svd>0) && opt_args.normfile == NULL) {
    opt_args.normfile=new char[strlen(argv[3])+5];
    sprintf(opt_args.normfile, "%s.std", argv[3]);
  }

  //get the training co-ordinate data, pre-process if necessary
  train=agf_get_features<real_a>(argv[1], &opt_args, ntrain, nvar);

  fprintf(logfs, "%d training vectors found: %s.vec\n", ntrain, argv[1]);

  //check the range of k:
  if (opt_args.k <= 0 || opt_args.k >= ntrain) {
    if (opt_args.k != -1) {
      //fprintf(stderr, "knn: Parameter k=%d out of range.  Using default (%d).\n", opt_args.k, K_DEFAULT_KNN);
      fprintf(stderr, "knn: Parameter k=%d out of range.  Exiting\n", opt_args.k);
      //opt_args.k=K_DEFAULT_KNN;
      exit(PARAMETER_OUT_OF_RANGE);
    }
  }

  //read test data:
  test=read_vecfile<real_a>(argv[2], ntest, nvar1);
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
  if (nvar1 != nvar && opt_args.normfile==NULL) {
    fprintf(stderr, "knn: Dimension of training data (%d) does not match dimension of test data (%d).\n", nvar, nvar1);
    exit(DIMENSION_MISMATCH);
  }

  //normalize test data:
  real_a **mat;
  if (opt_args.normfile!=NULL) {
    real_a *b;
    dim_ta nvar2, nvar3;
    real_a **testnew;

    mat=read_stats2(opt_args.normfile, b, nvar2, nvar3);
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
  }
   
  fprintf(logfs, "%d test vectors found in file %s\n", ntest, argv[1]);

  ordfile=new char[strlen(argv[1])+5];
  resultfile=new char[strlen(argv[3])+5];
  confile=new char[strlen(argv[3])+5];

  //action specific file stuff:
  switch (action) {
    case 0: 
      sprintf(resultfile, "%s.cls", argv[3]);
      sprintf(confile, "%s.con", argv[3]);
      sprintf(ordfile, "%s.cls", argv[1]);
      ord=read_clsfile<cls_ta>(ordfile, n1);
      ressize=sizeof(cls_ta);
      break;
    case 1:
      sprintf(resultfile, "%s.dat", argv[3]);
      sprintf(confile, "%s.err", argv[3]);
      sprintf(ordfile, "%s.dat", argv[1]);
      ord=read_datfile<real_a>(ordfile, n1);
      ressize=sizeof(real_a);
      break;
    case 2:
      sprintf(resultfile, "%s.dat", argv[3]);
      ord=NULL;
      ressize=sizeof(real_a);

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
        work=gsl_vector_alloc(nvar);
        S=gsl_vector_alloc(nvar);
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
    real_a p_x;
    switch (action) {
      case 0: 
        ((cls_ta *) result)[i]=knn(global_metric2_pointer_list[opt_args.metrictype], 
			train, nvar, ntrain, (cls_ta *) ord, nclass, 
			test[i], opt_args.k, pdf, opt_args.jointflag);
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
        ((real_a *) result)[i]=int_knn(global_metric2_pointer_list[opt_args.metrictype],
			train, nvar, ntrain, (real_a *) ord, test[i], opt_args.k, con[i]);
        break;
      case 2:
        ((real_a *) result)[i]=knn_pdf(train, nvar, ntrain, test[i], opt_args.k)/pcor;
        break;
    }
  }

  //write the results to a file:
  fs=fopen(resultfile, "w");
  if (fs == NULL) {
    fprintf(stderr, "Unable to open file for writing: %s\n", resultfile);
    return UNABLE_TO_OPEN_FILE_FOR_WRITING;
  }
  fwrite(result, ressize, ntest, fs);
  fclose(fs);

  if (action!=2) {
    fs=fopen(confile, "w");
    if (fs == NULL) {
      fprintf(stderr, "Unable to open file for writing: %s\n", confile);
      return UNABLE_TO_OPEN_FILE_FOR_WRITING;
    }
    fwrite(con, sizeof(real_a), ntest, fs);
    fclose(fs);
  }

  delete_matrix(train);
  if (action==0) delete [] (cls_ta *) ord;
  if (action==1) delete [] (cls_ta *) ord;
  free(result);
  delete_matrix(test);
  if (con!=NULL) delete [] con;
  if (pdf!=NULL) delete [] pdf;

  delete [] vecfile;
  delete [] ordfile;
  delete [] resultfile;
  delete [] confile;

  if (opt_args.normfile!=NULL) {
    delete [] opt_args.normfile;
    delete_matrix(mat);
  }

  ran_end();

  return err;

}
