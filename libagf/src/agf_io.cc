//
// This software is released under the following terms:
//
// 1. No commercial use.
// 2. Copies and derivative works are free to use and modify.
// 3. Attribution must be given to all contributors of both original and derivative works.
//
// Authors:
//
// 2017-07-16 Peter Mills: added license information 
//

//
// I/O routines for libAGF binary formats
//

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>

#include "error_codes.h"

#include "peteys_tmpl_lib.h"
#include "read_ascii_all.h"
#include "full_util.h"
#include "linked.h"

#include "agf_util.h"
#include "agf_io.h"

using namespace std;
using namespace libpetey;

namespace libagf {

#define STATS_HEADER "dim    average  std. dev.\n"

//-returns null pointer on failure
//-dimension is -1 if this is a bad value
//-sample number is -1 if there is an allocation failure
template <typename real>
real ** read_vecfile(const char *filename, nel_ta &m, dim_ta &n) {
  FILE *fs;
  real **data;
  nel_ta n1;			//dimensions have to be same type for matrix routines

  fs=fopen(filename, "r");
  if (fs == NULL) {
    m=0; n=0;
    return NULL;
  }

  //just piggy-back off the matrix routines:
  data=read_matrix<real, nel_ta>(fs, m, n1);
  n=n1;

  fclose(fs);

  return data;
}

template <typename cls_t>
cls_t * read_clsfile(const char *filename, nel_ta &n) {
  FILE *fs;
  cls_t *data;

  fs=fopen(filename, "r");
  if (fs == NULL) {
    n=0;
    return NULL;
  }

  fseek(fs, 0, SEEK_END);
  n=ftell(fs);
  if (n % sizeof(cls_ta) != 0) { //check for consistency
    n=-1;
    fclose(fs);
    return NULL;
  }
  n=n/sizeof(cls_ta);
  fseek(fs, 0, SEEK_SET);

  data=new cls_t [n];

  fread(data, sizeof(cls_t), n, fs);

  fclose(fs);

  return data;
}

template <typename real>
real * read_datfile(const char *filename, nel_ta &n) {
  FILE *fs;
  real *data;

  fs=fopen(filename, "r");
  if (fs == NULL) {
    n=0;
    return NULL;
  }

  fseek(fs, 0, SEEK_END);
  n=ftell(fs);
  if (n % sizeof(real) != 0) { //check for consistency
    n=-1;
    fclose(fs);
    return NULL;
  }
    
  n=n/sizeof(real_a);
  fseek(fs, 0, SEEK_SET);

  data=new real [n];
  if (data==NULL) {
    //doesn't c++ throw bad allocs anyway??
    n=-1;
    fclose(fs);
    return data;
  }

  fread(data, sizeof(real), n, fs);

  fclose(fs);

  return data;
}

#define MAXLL 200

int read_stats(const char *filename, real_a *ave, real_a *std, dim_ta ndim) {
  FILE *fs;
  dim_ta ivar;
  char header[MAXLL];

  fs=NULL;
  fs=fopen(filename, "r");
  if (fs==NULL) return UNABLE_TO_OPEN_FILE_FOR_READING;

  fgets(header, MAXLL, fs);
  for (dim_ta i=0; i<ndim; i++) {
    fscanf(fs, "%d %g %g", &ivar, &ave[i], &std[i]);
  }

  fclose(fs);

  return 0;
}

template <typename real>
real ** read_stats2(const char *filename, real *&ave, dim_ta &m, dim_ta &n) {
  FILE *fs;
  real **mat;
  dim_ta ivar;
  char *header=NULL;
  char **line=NULL;
  long nline;
  int ncon;

  ave=NULL;

  fs=fopen(filename, "r");
  if (fs==NULL) return NULL;

  header=fget_line(fs, 1);
  if (header==NULL) goto fail;
  if (strcmp(header, STATS_HEADER)==0) {
    line=read_ascii_all(fs, &nline, 1);
    if (line==NULL || nline<=0) goto fail;
    n=nline;
    ave=new real[n];
    mat=allocate_matrix<real, nel_ta>(n, n);
    for (dim_ta i=0; i<n; i++) {
      if (strlen(line[i])==0) {
        n=i+1;
        break;
      }
      ncon=sscanf(line[i], "%d %g %g", &ivar, ave+i, mat[i]+i);
      if (ncon!=3) {
        goto fail;
      }
      delete [] line[i];
    }
    m=n;
    delete [] line;
  } else if (header==NULL) {
    goto fail;
  } else {
    fseek(fs, 0, SEEK_SET);
    mat=read_matrix<real, nel_ta>(fs, m, n);
    if (mat==NULL || m<=0 || n<=0) goto fail;
    ave=new real[m];
    //printf("read_stats2: constant term:\n");
    n--;
    for (dim_ta i=0; i<m; i++) {
      ave[i]=mat[i][n];
      //fprintf(stderr, "%g\n", ave[i]);
    }
    //printf("read_stats2: transformation matrix:\n");
    //print_matrix(stdout, mat, m, n);
  }

  delete [] header;

  fclose(fs);

  return mat;

  fail:			//clean up in event of read failure and set error indicators
    fprintf(stderr, "read_stats2: an error occurred reading data file, %s\n", filename);
    if (line!=NULL) {
      for (long i=0; i<nline; i++) delete [] line[i];
      delete [] line;
    }
    if (ave!=NULL) delete [] ave;
    if (mat!=NULL) delete_matrix(mat);
    if (fs!=NULL) fclose(fs);
    if (header!=NULL) delete [] header;
    n=-1;
    m=-1;
  return NULL;
}

int print_stats(FILE *fs, float *ave, float *std, dim_ta ndim) {
  fprintf(fs, STATS_HEADER);
  for (dim_ta i=0; i<ndim; i++) {
    fprintf(fs, "%3d %10.6g %10.6g\n", i, ave[i], std[i]);
  }
  return 0;
}

int print_stats(FILE *fs, double *ave, double *std, dim_ta ndim) {
  fprintf(fs, STATS_HEADER);
  for (dim_ta i=0; i<ndim; i++) {
    fprintf(fs, "%4d %12.8lg %12.8lg\n", i, ave[i], std[i]);
  }
  return 0;
}

template <typename real, typename cls_t>
int agf_read_train(const char *fbase, real **&train, cls_t *&cls, nel_ta &n, dim_ta &nvar) {
  char *vecfile=NULL;
  char *classfile=NULL;
  nel_ta n1;
  int err=0;			//return code

  train=NULL;
  cls=NULL;

  vecfile=new char[strlen(fbase)+5];
  sprintf(vecfile, "%s.vec", fbase);

  classfile=new char[strlen(fbase)+5];
  sprintf(classfile, "%s.cls", fbase);

  //read in the training data:
  train=read_vecfile<real>(vecfile, n, nvar);
  if (nvar <= 0 || n <= 0) {
    fprintf(stderr, "Error reading file: %s\n", vecfile);
    err=FILE_READ_ERROR;
    goto fail;
  }
  if (train==NULL) {
    fprintf(stderr, "Unable to open file, %s, for reading.\n", vecfile);
    err=UNABLE_TO_OPEN_FILE_FOR_READING;
    goto fail;
  }

  cls=read_clsfile<cls_t>(classfile, n1);

  if (n1 <= 0) {
    fprintf(stderr, "Error reading file: %s\n", classfile);
    err=FILE_READ_ERROR;
    goto fail;
  }
  if (cls == NULL) {
    fprintf(stderr, "Unable to open file, %s, for reading.\n", classfile);
    err=UNABLE_TO_OPEN_FILE_FOR_READING;
    goto fail;
  }
  if (n1!=n) {
    fprintf(stderr, "Sample count mismatch: %d in %s, %d in %s.\n", n, vecfile, n1, classfile);
    err=SAMPLE_COUNT_MISMATCH;
    goto fail;
  }

  fail:

    delete [] classfile;
    delete [] vecfile;

  return err;
}	     

template <typename real>
int agf_read_borders(const char *fbase, real **&brd, real **&grd, nel_ta &n, dim_ta &nvar) {
  char *brdfile;
  char *grdfile;
  nel_ta n1;
  dim_ta nvar1;
  int err=0;

  brd=NULL;
  grd=NULL;

  brdfile=new char[strlen(fbase)+5];
  strcpy(brdfile, fbase);
  strcat(brdfile, ".brd");

  grdfile=new char[strlen(fbase)+5];
  strcpy(grdfile, fbase);
  strcat(grdfile, ".bgd");

  //read in the decision surface data:
  brd=read_vecfile<real>(brdfile, n, nvar);
  if (nvar <= 0 || n <= 0) {
    fprintf(stderr, "Error reading file: %s\n", brdfile);
    err=FILE_READ_ERROR;
    goto fail;
  }
  if (brd == NULL) {
    fprintf(stderr, "Unable to open input file: %s\n", brdfile);
    err=UNABLE_TO_OPEN_FILE_FOR_READING;
    goto fail;
  }

  grd=read_vecfile<real>(grdfile, n1, nvar1);

  if (nvar1 <= 0 || n1 <= 0) {
    fprintf(stderr, "Error reading file: %s\n", grdfile);
    err=FILE_READ_ERROR;
    goto fail;
  }
  if (grd == NULL) {
    fprintf(stderr, "Unable to open input file: %s\n", grdfile);
    err=UNABLE_TO_OPEN_FILE_FOR_READING;
    goto fail;
  }
  if (nvar1 != nvar) {
    fprintf(stderr, "Error: dimensions of border and gradient vectors do not agree:\n");
    fprintf(stderr, "       %s: D=%d vs. %s: D=%d\n", brdfile, nvar, grdfile, nvar1);
    err=DIMENSION_MISMATCH;
    goto fail;
  }
  if (n1 != n) {
    fprintf(stderr, "Error: number of samples in border and gradient files do not agree:\n");
    fprintf(stderr, "       %d in %s vs. %d in %s\n", n, brdfile, n1, grdfile);
    err=SAMPLE_COUNT_MISMATCH;
    goto fail;
  }

  fail:
  
    delete [] grdfile;
    delete [] brdfile;

  return err;
}

char *compile_precommand(const char *fname, agf_command_opts *optargs) {
  char *command=NULL;

  if (optargs->normflag || optargs->svd>0 || optargs->normfile!=NULL) {
    //if the user wants some pre-processing we farm it out to "agf_precondition"

    //if (optargs.normfile == NULL) {
    //  optargs.normfile=new char[strlen(argv[3])+5];
    //  sprintf(optargs.normfile, "%s.std", argv[3]);
    //}
    command=new char[strlen(optargs->normfile)+strlen(fname)+100];
    sprintf(command, "%s%s%s -a %s", AGF_COMMAND_PREFIX, 
		AGF_LTRAN_COM, AGF_OPT_VER, optargs->normfile);
    if (optargs->normflag) strcat(command, " -n");
    if (optargs->asciiflag) strcat(command, " -A");
    if (optargs->Mflag) strcat(command, " -M");
    if (optargs->svd>0) {
      sprintf(command+strlen(command), " -S %d", optargs->svd);
    }
    sprintf(command+strlen(command), " %s", fname);
  }
  
  return command;
  
}

//increasingly I'm finding this bit really brain-dead:
template <typename real>
real **agf_get_features(const char *fbase, agf_command_opts *opt_args, nel_ta &n, dim_ta &nvar, flag_a sufflag)
{
  real **train;
  char *vecfile;

  vecfile=new char[strlen(fbase)+5];
  if (sufflag) {
    sprintf(vecfile, "%s", fbase);
  } else {
    sprintf(vecfile, "%s.vec", fbase);
  }

  if (opt_args->normflag || opt_args->svd>0 || opt_args->normfile != NULL) {
    //if the user wants some pre-processing we farm it out to "agf_precondition"
    FILE *fs;
    char *command;
    nel_ta nvar1;

    //if (opt_args.normfile == NULL) {
    //  opt_args.normfile=new char[strlen(argv[3])+5];
    //  sprintf(opt_args.normfile, "%s.std", argv[3]);
    //}
    command=new char[strlen(opt_args->normfile)+strlen(fbase)+50];
    sprintf(command, "%s%s%s -a %s", AGF_COMMAND_PREFIX, 
		AGF_LTRAN_COM, AGF_OPT_VER, opt_args->normfile);
    if (opt_args->normflag) strcat(command, " -n");
    //if (opt_args->asciiflag) strcat(command, " -A");
    //if (opt_args->Mflag) strcat(command, " -M");
    if (opt_args->svd>0) {
      sprintf(command+strlen(command), " -S %d", opt_args->svd);
    }
    sprintf(command+strlen(command), " %s", vecfile);
    fprintf(stderr, "%s\n", command);
    fs=popen(command, "r");
    train=read_matrix<real, nel_ta>(fs, n, nvar1);
    nvar=nvar1;
    pclose(fs);
    delete [] command;
  } else {
    //read in the training data:
    train=read_vecfile<real>(vecfile, n, nvar);
  }
  if (nvar <= 0 || n <= 0) {
    fprintf(stderr, "Error reading file: %s\n", vecfile);
    exit(FILE_READ_ERROR);
  }
  if (train==NULL) {
    fprintf(stderr, "Unable to open file, %s, for reading.\n", vecfile);
    exit(UNABLE_TO_OPEN_FILE_FOR_READING);
  }

  delete [] vecfile;
  return train;

}

template float **read_vecfile<float>(const char *, dim_ta &, nel_ta &);
template double **read_vecfile<double>(const char *, dim_ta &, nel_ta &);

template cls_ta *read_clsfile<cls_ta>(const char *, nel_ta &);

template float *read_datfile<float>(const char *, nel_ta &);
template double *read_datfile<double>(const char *, nel_ta &);

template float **read_stats2<float>(const char *, float*&, dim_ta &, dim_ta &);
template double **read_stats2<double>(const char *, double*&, dim_ta &, dim_ta &);

template int agf_read_train<float, cls_ta>(const char *, float **&, cls_ta *&, nel_ta &, dim_ta &);
template int agf_read_train<double, cls_ta>(const char *, double **&, cls_ta *&, nel_ta &, dim_ta &);

template int agf_read_borders<float>(const char *, float **&, float **&, nel_ta &, dim_ta &);
template int agf_read_borders<double>(const char *, double **&, double **&, nel_ta &, dim_ta &);

template float **agf_get_features<float>(const char *, agf_command_opts *, nel_ta &, dim_ta &, flag_a);
template double **agf_get_features<double>(const char *, agf_command_opts *, nel_ta &, dim_ta &, flag_a);

} //end namespace libagf

