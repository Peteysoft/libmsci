
//Most common pre-processing operation--normalization, SVD, removing a feature
//--are 1. linear operations, 2. require the test data to be multiplied by
//a "pre-condtioning" matrix

#include <math.h>
#include <string.h>
#include <stdio.h>

#include <gsl/gsl_linalg.h>
#include <assert.h>

#include "full_util.h"
#include "agf_lib.h"

using namespace std;
using namespace libagf;
using namespace libpetey;

int main(int argc, char *argv[]) {
  char *infile=NULL;			//training data features only
  char *resultfile=NULL;		//transformed training data
  FILE *fs;
  FILE *diagfs;			//print diagnostics to this file stream

  nel_ta ntrain;		//number of training data points
  dim_ta nvar;			//number of variables in original data
  dim_ta nvar2;			//number of variables after feature selection
  dim_ta nvar3;			//number of variables after SVD

  real_a **train;		//training data vectors
  real_a **result;		//the result
  real_a **mat;			//transformation "conditioning" matrix

  cls_ta *cls=NULL;		//class data: not used

  int exit_value;		//return value

  //components of the transformation matrix:
  gsl_matrix *v;		//right singular vectors
  gsl_vector *s;		//singular values
  real_a *std, *ave;		//standard deviations and averages
  cls_ta *ind;			//freature selection indices

  agf_command_opts opt_args;

  //parse out command switches:
  exit_value=agf_parse_command_opts(argc, argv, "01a:AnS:FME:UCH", &opt_args);
  if (exit_value==FATAL_COMMAND_OPTION_PARSE_ERROR) return exit_value;

  //print help screen:
  if (argc < 1 && (opt_args.stdinflag==0 || opt_args.stdoutflag==0)) {
    printf("\n");
    printf("Syntax:	agf_precondition -a normfile [-n] [-S nsv] [-F] \\\n");
    printf("                          [-A [-M [-E]] [-C] [-H]] {-0 | input} {-1 | output} \\\n");
    printf("                          [ind1 [ind2 [ind3...]]]\n");
    printf("\n");
    printf("arguments:\n");
    printf("  normfile   input/output transformation matrix\n");
    printf("  input      input file containing vector data\n");
    printf("  output     output file containing transformed vector data\n");
    printf("  indN       feature selection index\n");
    printf("\n");
    printf("file options:\n");
    printf("  -0         read from stdin\n");
    printf("  -1         write to stdout\n");
    printf("  -A         input and output files are in ASCII (LVQ) format\n");
    printf("  -M         (in conjunction with -A) specifies LIBSVM format\n");
    printf("  -C         (in conjunction with -A) no class data\n");
    printf("  -H         (in conjunction with -A) no header\n");
    printf("  -E         (in conjunction with -A and -M) value for missing data\n");
    printf("\n");
    printf("operations (in order of execution):\n");
    printf("  -F         select features\n");
    printf("  -n         normalize with standard deviations\n");
    printf("  -S svd     singular value decomposition (SVD); keep top nsv singular values\n");
    //printf("    -N take data from stdin, write to stdout\n");
    printf("\n");
    return INSUFFICIENT_COMMAND_ARGS;
  }

  //must specify a normalization file:
  if (opt_args.normfile==NULL) {
    fprintf(stderr, "agf_precondition: must specify normalization file with -a\n");
    exit(INSUFFICIENT_COMMAND_ARGS);
  }

  //where to stick error messages:
  diagfs=stderr;

  //if there are no arguments or -0 flag set, we read from stdin
  if (argc==0 || opt_args.stdinflag) {
    fs=stdin;
  } else {
    infile=argv[0];
    fs=fopen(infile, "r");
    if (fs == NULL) {
      fprintf(stderr, "agf_precondition: Unable to open file, %s, for reading\n", infile);
      exit(UNABLE_TO_OPEN_FILE_FOR_READING);
    }
  }

  //ascii versus binary files:
  if (opt_args.asciiflag) {
    int readflag=opt_args.Hflag+2*opt_args.Cflag;
    if (opt_args.Mflag) {
      ntrain=read_svm(fs, train, cls, nvar, opt_args.missing, opt_args.Uflag);
    } else {
      ntrain=read_lvq(fs, train, cls, nvar, readflag);
    }
    if (ntrain<=0) {
      fprintf(stderr, "agf_precondtion: Read error\n");
      exit(FILE_READ_ERROR);
    }
  } else {
    cls=NULL;
    nel_ta nv1;
    //read in the training data:
    train=read_matrix<real_a, nel_ta>(fs, ntrain, nv1);
    nvar=nv1;
	//really, why the f* have I duplicated this function??
    //train=read_vecfile(infile, ntrain, nvar);
    fprintf(diagfs, "%d %d-dimensional training vectors found in file: %s\n", ntrain, nvar, infile);
  }
  if (opt_args.stdinflag==0) fclose(fs);

  nvar2=nvar;

  //feature selection:
  ind=new dim_ta[nvar];
  if (opt_args.selectflag) {
    real_a **result2;

    if (opt_args.stdinflag && opt_args.stdoutflag) {
      nvar2=argc;
    } else if (opt_args.stdinflag || opt_args.stdoutflag) {
      nvar2=argc-1;
    } else {
      nvar2=argc-2;
    }

    if (nvar2 <= 0) {
      fprintf(stderr, "agf_precondition: no selection terms in argument list\n");
      exit(PARAMETER_OUT_OF_RANGE);
    }
    for (dim_ta i=0; i<nvar2; i++) {
      if (sscanf(argv[i+argc-nvar2], "%d", ind+i)!=1) {
        fprintf(stderr, "agf_precondition: unable to read seleciton index %d\n", i);
        exit(FATAL_COMMAND_OPTION_PARSE_ERROR);
      }
      if (ind[i]>=nvar || ind[i]<0) {
        fprintf(stderr, "agf_precondition: index out of range\n");
        exit(PARAMETER_OUT_OF_RANGE);
      }
    }
    result2=allocate_matrix<real_a, int32_t>(ntrain, nvar2);
    for (nel_ta i=0; i<ntrain; i++) {
      for (dim_ta j=0; j<nvar2; j++) {
        result2[i][j]=train[i][ind[j]];
      }
    }
    delete_matrix(train);
    train=result2;
  } else {
    //if there is no feature selection: pass along the identity matrix...
    nvar2=nvar;
    for (dim_ta i=0; i<nvar2; i++) ind[i]=i;
  }

  //calculate the averages and standard deviations:
  std=new real_a[nvar2];
  ave=new real_a[nvar2];
  if (opt_args.normflag) {
    dim_ta nzerostd=0;
    calc_norm(train, nvar2, ntrain, ave, std);
    //want to remove variables with standard-deviation of zero:
    for (dim_ta i=0; i<nvar2; i++) {
      ind[i-nzerostd]=ind[i];
      if (std[i]==0) {
        nzerostd++;
      }
    }
    for (nel_ta i=0; i<ntrain; i++) {
      for (dim_ta j=0; j<nvar2; j++) {
        if (std[j]!=0) train[i][j]=(train[i][j]-ave[j])/std[j];
      }
    }
    print_stats(diagfs, ave, std, nvar2);
    if (nzerostd!=0) {
      //not terribly efficient, but expedient:
      dim_ta m=0;
      real_a **result2;
      result2=allocate_matrix<real_a, int32_t>(ntrain, nvar2);
      for (dim_ta j=0; j<nvar2; j++) {
        for (nel_ta i=0; i<ntrain; i++) {
          result2[i][m]=train[i][j];
        }
        if (std[j]!=0) m++;
      }
      delete_matrix(train);
      train=result2;
      nvar2-=nzerostd;
    }
    //if (argc>=2) print_stats(diagfs, ave, std, nvar2);
  } else {
    //(kind of a stupid way of doing it... oh well)
    for (dim_ta i=0; i<nvar2; i++) {
      std[i]=1;
      ave[i]=0;
    }
  }

  //singular value decomposition:
  if (opt_args.svd>0) {
    gsl_matrix *u;
    gsl_vector *work;
    nel_ta k;

    //always remove averages (since they are wasted...)
    if (opt_args.normflag==0) {
      real_a dum[nvar2];
      calc_norm(train, nvar2, ntrain, ave, dum);
      for (nel_ta i=0; i<ntrain; i++) {
        for (dim_ta j=0; j<nvar2; j++) {
          train[i][j]=train[i][j]-ave[j];
        }
      }
    }

    if (ntrain>nvar2) {
      u=gsl_matrix_alloc(ntrain, nvar2);
      for (nel_ta i=0; i<ntrain; i++) {
        for (dim_ta j=0; j<nvar2; j++) {
          gsl_matrix_set(u, i, j, train[i][j]);
        }
      }
      k=nvar2;
    } else {
      u=gsl_matrix_alloc(nvar2, ntrain);
      for (nel_ta i=0; i<ntrain; i++) {
        for (dim_ta j=0; j<nvar2; j++) {
          gsl_matrix_set(u, j, i, train[i][j]);
        }
      }
      k=ntrain;
    }

    v=gsl_matrix_alloc(k, k);
    s=gsl_vector_alloc(k);
    work=gsl_vector_alloc(k);

    fprintf(diagfs, "agf_precondition: calling GSL SVD subroutine...\n");
    //if (argc>=2) fprintf(diagfs, "agf_precondition: calling GSL SVD subroutine...\n");
    //gsl_linalg_SV_decomp(u, v, s, work);
    gsl_linalg_SV_decomp_jacobi (u, v, s);

    if (opt_args.svd>0 && opt_args.svd<k) {
      nvar3=opt_args.svd;
    } else {
      nvar3=k;
    }

    result=allocate_matrix<real_a, nel_ta>(ntrain, nvar3);
    if (ntrain>nvar2) {
      for (nel_ta i=0; i<ntrain; i++) {
        for (dim_ta j=0; j<nvar3; j++) {
          result[i][j]=gsl_matrix_get(u, i, j)*gsl_vector_get(s, j);
        }
      }
    } else {
      for (nel_ta i=0; i<ntrain; i++) {
        for (nel_ta j=0; j<nvar3; j++) {
          result[i][j]=gsl_matrix_get(v, j, i)*gsl_vector_get(s, j);
        }
      }
      gsl_matrix_free(v);
      v=gsl_matrix_alloc(ntrain, nvar2);
      gsl_matrix_transpose_memcpy(v, u);
    }
    /*
    if (argc>=2) {
      for (nel_ta i=0; i<k; i++) printf("%g\n", gsl_vector_get(s, i));
      for (nel_ta i=0; i<k; i++) {
        for (nel_ta j=0; j<k; j++) {
          printf("%12.6g ", gsl_matrix_get(v, i, j));
        }
        printf("\n");
      }
    }
    */

    gsl_matrix_free(u);
    gsl_vector_free(work);
  } else {
    //if there is no singular value decomposition, set v to the identity:
    v=gsl_matrix_alloc(nvar2, nvar2);
    gsl_matrix_set_identity(v);
    s=gsl_vector_alloc(nvar2);
    gsl_vector_set_all(s, 1);
    result=copy_matrix(train, ntrain, nvar2);
    nvar3=nvar2;
  }

  if (opt_args.selectflag==0 && opt_args.normflag==0 && opt_args.svd<=0) {
    //features data is transformed strictly from data in an external file:
    if (opt_args.normfile!=NULL) {
      //read in the normalization data and use it transform the data
      //rather than printing it out:
      mat=read_stats2(opt_args.normfile, ave, nvar2, nvar3);
      //mat=read_matrix<real_a, nel_ta>(fs, nvar2, nvar3);
      assert(nvar==nvar2);
      for (nel_ta i=0; i<ntrain; i++) 
		for (dim_ta j=0; j<nvar; j++) train[i][j]-=ave[j];
      result=matrix_mult(train, mat, ntrain, nvar, nvar3);
    } else {
      //congrats, you just wasted some compute cycles...
      result=copy_matrix<real_a, int32_t>(train, ntrain, nvar);
      nvar3=nvar;
    }
  } else {
    //here we multiply everything together: feature selection, SVD and normalization
    mat=zero_matrix<real_a, nel_ta>(nvar, nvar3+1);
    for (dim_ta i=0; i<nvar2; i++) {
      for (dim_ta j=0; j<nvar3; j++) {
        mat[ind[i]][j]=gsl_matrix_get(v, i, j)/std[i];
      }
      mat[ind[i]][nvar3]=ave[i];	//store averages to right of matrix
    }
  }

  //if the second argument is missing, write to stdout:
  if (opt_args.stdoutflag || argc<2) {
    fs=stdout;
  } else {
    resultfile=argv[1];
    fs=fopen(resultfile, "w");
    if (fs==NULL) {
      fprintf(stderr, "Unable to open file for writing: %s\n", resultfile);
      return UNABLE_TO_OPEN_FILE_FOR_WRITING;
    }
  }

  //write the results to a file:
  if (opt_args.asciiflag) {
    print_lvq_svm(fs, result, cls, ntrain, nvar3, opt_args.Mflag, opt_args.Hflag);
  } else {
    fwrite(&nvar3, sizeof(nvar3), 1, fs);
    fwrite(result[0], sizeof(real_a), nvar3*ntrain, fs);
  }
  if (opt_args.stdoutflag==0 && argc>=2) fclose(fs);

  //write transformation matrix to specified file:
  if (opt_args.normfile!=NULL) {
    fs=fopen(opt_args.normfile, "w");
    if (fs==NULL) {
      fprintf(stderr, "Unable to open file for writing: %s\n", opt_args.normfile);
      return UNABLE_TO_OPEN_FILE_FOR_WRITING;
    }
    nvar3++;
    fwrite(&nvar3, sizeof(nvar3), 1, fs);
    fwrite(mat[0], sizeof(real_a), nvar3*nvar, fs);
    fclose(fs);
  }

  //clean up:
  delete [] train[0];
  delete [] train;

  if (cls!=NULL) delete [] cls;

  delete_matrix(mat);
  delete_matrix(result);

  delete [] std;
  delete [] ave;

  delete [] ind;

  gsl_vector_free(s);
  gsl_matrix_free(v);

  if (opt_args.normfile!=NULL) delete [] opt_args.normfile;

  return exit_value;

}


