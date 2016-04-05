
//Copyright (C) 2007 Peter Mills.  All rights reserved.

#include <math.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>
#include <sys/timeb.h>

#include <vector>

#include "kextreme.h"
#include "quicksort.h"
//#include "coeffs.h"
//#include "nr.h"
#include "peteys_tmpl_lib.h"
#include "agf_lib.h"

#include "gsl/gsl_rng.h"

#include "agf_util.h"

#include <stdio.h>

using namespace std;
using namespace libpetey;

namespace libagf {

const gsl_rng_type *agf_gsl_rng_type=gsl_rng_mt19937;

//to prevent the gsl library from aborting if it encounters an error:
void agf_gsl_handler(const char *reason, const char * file, int line, int gsl_errno) {
  fprintf(stderr, "Error in %s, line %d: %s\n", file, line, reason);
}

//this whole thing needs to be redone... (does not scale well)
int agf_parse_command_opts(int &argc, char **&argv, const char *optlist, agf_command_opts *opt_args) {
  char c;
  int ncon;
  int errcode;

  int32_t maxitr;
  float tol;

  //the defaults for these parameters are set in the main routines:
  float var_0;
  float var2;
  float Wc;
  float wmax;
  enum_a Qtype;

  //set defaults
  errcode=0;			//no errors

  //"output" file:
  opt_args->ofile=NULL;

  opt_args->multicommand=NULL;	//-O

  //flags:
  opt_args->stdinflag=0;	//-0
  opt_args->stdoutflag=0;	//-1
  opt_args->asciiflag=0;	//-A
  opt_args->Bflag=0;		//-B
  opt_args->Cflag=0;		//-C
  opt_args->errflag=0;		//-e
  opt_args->Gflag=0;		//-G
  opt_args->Hflag=0;		//-H
  opt_args->jointflag=0;	//-j
  opt_args->Kflag=0;		//-K
  opt_args->Lflag=0;		//-L
  opt_args->Mflag=0;		//-L
  opt_args->normflag=0;		//-n
  opt_args->Nflag=0;		//-N
  opt_args->Pflag=0;		//-P
  opt_args->Rflag=0;		//-R
  opt_args->uflag=0;		//-u
  opt_args->Uflag=0;		//-U
  opt_args->xflag=0;		//-x
  opt_args->Yflag=0;		//-Y
  opt_args->zflag=0;		//-z
  opt_args->Zflag=0;		//-Z

  //enumerated selections:
  //opt_args->algtype=-1;		//-c
  opt_args->metrictype=0;	//-m

  //integer parameters:
  //opt_args->k=K_DEFAULT_KNN;	//-k
  //for (int j=0; optlist[j]!=NULL; j++) if (optlist[j]=='w') opt_args->k=K_DEFAULT_AGF;
  //number of trials for AGF optimal
  opt_args->nt=NT_DEFAULT;	//-q

  //floating point parameters:
  //for AGF:
  opt_args->missing=0.;		//-E
  //opt_args->Wc=W_DEFAULT;	//-W
  opt_args->var[0]=-1;		//-v
  opt_args->var[1]=-1;		//-V

  //for AGF borders:
  opt_args->hrel=HREL_DEFAULT;	//-h
  opt_args->tol=TOL_DEFAULT;	//-t
  opt_args->n=NBORD_DEFAULT;	//-s
  opt_args->rthresh=RTHRESH_DEFAULT;	//-r
  opt_args->cl_thresh=1;		//-T

  //for clustering:
  opt_args->pmin=-1.;		//-p

  //for AGF optimal:
  opt_args->W1=-1;		//-w

  //for validation routines:
  opt_args->ftest=F_DEFAULT;	//-f
  opt_args->fflag=0;
  opt_args->div=-1;		//-d

  //string parameters:
  opt_args->normfile=NULL;	//-a

  //pre-processing:
  opt_args->svd=-1;		//-S
  opt_args->selectflag=0;	//-F

  while ((c = getopt(argc, argv, optlist)) != -1) {
    switch (c) {
      case ('0'):
             opt_args->stdinflag=1;
	     break;
      case ('1'):
             opt_args->stdoutflag=1;
	     break;
      case ('a'):
             opt_args->normfile=new char[strlen(optarg)+1];
             strcpy(opt_args->normfile, optarg);
	     break;
      case ('A'):
             opt_args->asciiflag=1;
	     break;
      case ('B'):
             opt_args->Bflag=1;
	     break;
      case ('c'):
             ncon=sscanf(optarg, "%d", &opt_args->algtype);
             if (ncon != 1) {
               fprintf(stderr, "Warning: garbled option argument: -c %s", optarg);
               errcode=COMMAND_OPTION_PARSE_ERROR;
             }
	     break;
      case ('C'):
             opt_args->Cflag=1;
	     break;
      case ('d'):
             ncon=sscanf(optarg, "%d", &opt_args->div);
             if (ncon != 1) {
               fprintf(stderr, "Warning: garbled option argument: -d %s", optarg);
               errcode=COMMAND_OPTION_PARSE_ERROR;
             }
	     break;
      case ('e'):
             opt_args->errflag=1;
	     break;
      case ('E'):
             ncon=sscanf(optarg, "%f", &opt_args->missing);
             if (ncon != 1) {
               fprintf(stderr, "Warning: garbled option argument: -E %s", optarg);
               errcode=COMMAND_OPTION_PARSE_ERROR;
             }
             break;
      case ('f'):
             ncon=sscanf(optarg, "%f", &opt_args->ftest);
             if (ncon != 1) {
               fprintf(stderr, "Warning: garbled option argument: -f %s", optarg);
               errcode=COMMAND_OPTION_PARSE_ERROR;
	     } else if (opt_args->ftest <= 0.) {
               fprintf(stderr, "Warning: parameter ftest=%g out of range, using default\n", opt_args->ftest);
               errcode=PARAMETER_OUT_OF_RANGE;
               opt_args->ftest=F_DEFAULT;
             }
             opt_args->fflag=1;
	     break;
      case ('F'):
             opt_args->selectflag=1;
	     break;
      case ('G'):
             opt_args->Gflag=1;
             break;
      case ('h'):
             ncon=sscanf(optarg, "%g", &opt_args->hrel);
             if (ncon != 1) {
               fprintf(stderr, "Warning: garbled option argument: -h %s", optarg);
               opt_args->hrel=HREL_DEFAULT;
               errcode=COMMAND_OPTION_PARSE_ERROR;
             }
	     break;
      case ('H'):
             opt_args->Hflag=1;
             break;
      case ('i'):
             ncon=sscanf(optarg, "%d", &maxitr);
             if (ncon != 1) {
               fprintf(stderr, "Warning: garbled option argument: -i %s", optarg);
               errcode=COMMAND_OPTION_PARSE_ERROR;
	     } else if (maxitr <= 0.) {
               fprintf(stderr, "Warning: parameter borders_maxiter=%d out of range, using default\n", maxitr);
               errcode=PARAMETER_OUT_OF_RANGE;
             } else {
               agf_global_borders_maxiter=maxitr;
             }
	     break;
      case ('I'):
             ncon=sscanf(optarg, "%d", &maxitr);
             if (ncon != 1) {
               fprintf(stderr, "Warning: garbled option argument: -I %s", optarg);
               errcode=COMMAND_OPTION_PARSE_ERROR;
	     } else if (maxitr <= 0.) {
               fprintf(stderr, "Warning: parameter weights_maxiter=%d out of range, using default\n", maxitr);
               errcode=PARAMETER_OUT_OF_RANGE;
             } else {
               agf_global_weights_maxiter=maxitr;
             }
	     break;
      case ('j'):
             opt_args->jointflag=1;
	     break;
      case ('k'):
             ncon=sscanf(optarg, "%d", &opt_args->k);
             if (ncon != 1) {
               fprintf(stderr, "Warning: garbled command option argument: -k %s", optarg);
               errcode=COMMAND_OPTION_PARSE_ERROR;
	     }
	     break;
      case ('K'):
             opt_args->Kflag=1;
	     break;
      case ('l'):
             ncon=sscanf(optarg, "%f", &tol);
             if (ncon != 1) {
               fprintf(stderr, "Warning: garbled option argument: -l %s", optarg);
               errcode=COMMAND_OPTION_PARSE_ERROR;
	     } else if (tol <= 0.) {
               fprintf(stderr, "Warning: parameter weights_maxiter=%g out of range, using default\n", tol);
               errcode=PARAMETER_OUT_OF_RANGE;
             } else {
               agf_global_weights_tol=tol;
             }
	     break;
      case ('L'):
             opt_args->Lflag=1;
	     break;
      case ('m'):
             ncon=sscanf(optarg, "%d", &opt_args->metrictype);
             if (ncon != 1) {
               fprintf(stderr, "Warning: garbled option argument: -m %s", optarg);
               errcode=COMMAND_OPTION_PARSE_ERROR;
	     } else if (opt_args->metrictype < 0 || opt_args->metrictype >= global_nmetric) {
               fprintf(stderr, "Warning: parameter, metric-type=%d, out of range\n", opt_args->metrictype);
               fprintf(stderr, "         using default Cartesian metric\n");
	       opt_args->metrictype=0;
               errcode=PARAMETER_OUT_OF_RANGE;
             }
	     //is it wise to set the metric right here??
	     global_metric2=global_metric2_pointer_list[opt_args->metrictype];
	     break;
      case ('M'):
             opt_args->Mflag=1;
	     break;
      case ('n'):
             opt_args->normflag=1;
	     break;
      case ('N'):
             ncon=sscanf(optarg, "%d", &maxitr);
             if (ncon != 1) {
               fprintf(stderr, "Warning: garbled option argument: -N %s", optarg);
               errcode=COMMAND_OPTION_PARSE_ERROR;
	     } else if (maxitr <= 0.) {
               fprintf(stderr, "Warning: parameter weights_maxiter=%d out of range, using default\n", maxitr);
               errcode=PARAMETER_OUT_OF_RANGE;
             } else {
               agf_global_borders_maxiter=maxitr;
               agf_global_weights_maxiter=maxitr;
             }
	     break;
      case ('o'):
             opt_args->ofile=new char[strlen(optarg)+1];
             strcpy(opt_args->ofile, optarg);
	     break;
      case ('O'):
             opt_args->multicommand=new char[strlen(optarg)+1];
             strcpy(opt_args->multicommand, optarg);
	     break;
      case ('p'):
             ncon=sscanf(optarg, "%f", &opt_args->pmin);
             if (ncon != 1) {
               fprintf(stderr, "Warning: garbled option argument: -p %s", optarg);
               errcode=COMMAND_OPTION_PARSE_ERROR;
	     } else if (opt_args->pmin <= 0.) {
               fprintf(stderr, "Warning: parameter P_min=%g out of range\n", opt_args->pmin);
	       opt_args->pmin=-1;
               errcode=PARAMETER_OUT_OF_RANGE;
             }
	     break;
      case ('P'):
             opt_args->Pflag=1;
	     break;
      case ('q'):
             ncon=sscanf(optarg, "%d", &opt_args->nt);
             if (ncon != 1) {
               fprintf(stderr, "Warning: garbled option argument: -q %s", optarg);
               errcode=COMMAND_OPTION_PARSE_ERROR;
	     } else if (opt_args->nt < 4) {
               fprintf(stderr, "Warning: parameter nt=%d out of range\n", opt_args->nt);
               errcode=PARAMETER_OUT_OF_RANGE;
               opt_args->nt=NT_DEFAULT;
             }
	     break;
      case ('Q'):
             ncon=sscanf(optarg, "%d", &Qtype);
             if (ncon != 1) {
               fprintf(stderr, "Warning: garbled option argument: -Q %s", optarg);
               errcode=COMMAND_OPTION_PARSE_ERROR;
             } else {
               opt_args->Qtype=Qtype;
             }
             break;
      case ('r'):
             ncon=sscanf(optarg, "%f", &opt_args->rthresh);
             if (ncon != 1) {
               fprintf(stderr, "Warning: garbled option argument: -r %s", optarg);
               errcode=COMMAND_OPTION_PARSE_ERROR;
	     } else if (opt_args->rthresh <= -1. || opt_args->rthresh > 1.) {
               fprintf(stderr, "Warning: parameter r0=%g out of range, using default\n", opt_args->rthresh);
               errcode=PARAMETER_OUT_OF_RANGE;
               opt_args->rthresh=RTHRESH_DEFAULT;
             }
	     break;
      case ('R'):
             opt_args->Rflag=1;
	     break;
      case ('s'):
             ncon=sscanf(optarg, "%d", &opt_args->n);
             if (ncon != 1) {
               fprintf(stderr, "Warning: garbled option argument: -s %s", optarg);
               errcode=COMMAND_OPTION_PARSE_ERROR;
	     } else if (opt_args->n <= 0) {
               fprintf(stderr, "Warning: parameter n=%d out of range, using default=%d\n", opt_args->n, opt_args->n);
               errcode=PARAMETER_OUT_OF_RANGE;
	       opt_args->n=NBORD_DEFAULT;
             }
	     break;
      case ('S'):
             ncon=sscanf(optarg, "%d", &opt_args->svd);
             if (ncon != 1) {
               fprintf(stderr, "Warning: garbled option argument: -S %s", optarg);
               errcode=COMMAND_OPTION_PARSE_ERROR;
             }
	     break;
      case ('t'):
             ncon=sscanf(optarg, "%g", &opt_args->tol);
             if (ncon != 1) {
               fprintf(stderr, "Warning: garbled option argument: -t %s", optarg);
               errcode=COMMAND_OPTION_PARSE_ERROR;
	     } else if (opt_args->tol <= 0) {
               fprintf(stderr, "Warning: parameter tol=%g out of range, using default=%g\n", opt_args->tol, TOL_DEFAULT);
               errcode=PARAMETER_OUT_OF_RANGE;
               opt_args->tol=TOL_DEFAULT;
             }
	     break;
      case ('T'):
             ncon=sscanf(optarg, "%d", &opt_args->cl_thresh);
             if (ncon != 1) {
               fprintf(stderr, "Warning: garbled option argument: -T %s", optarg);
               errcode=COMMAND_OPTION_PARSE_ERROR;
	     } else if (opt_args->cl_thresh < 0) {
               fprintf(stderr, "Warning: parameter cl-thresh=%d out of range, using default=%d\n", opt_args->cl_thresh, 1);
               errcode=PARAMETER_OUT_OF_RANGE;
               opt_args->cl_thresh=1;
             }
	     break;
      case ('u'):
             opt_args->uflag=1;
	     break;
      case ('U'):
             opt_args->Uflag=1;
	     break;
      case ('v'):
             ncon=sscanf(optarg, "%f", &var_0);
             if (ncon != 1) {
               fprintf(stderr, "Warning: garbled option argument: -v %s", optarg);
               errcode=COMMAND_OPTION_PARSE_ERROR;
	     } else if (var_0 <= 0. && var_0!=-1) {
               fprintf(stderr, "Warning: parameter var_0=%g out of range, using default\n", var_0);
               errcode=PARAMETER_OUT_OF_RANGE;
             } else {
               opt_args->var[0]=var_0;
             }
	     break;
      case ('V'):
             ncon=sscanf(optarg, "%f", &var2);
             if (ncon != 1) {
               fprintf(stderr, "Warning: garbled option argument: -v %s", optarg);
               errcode=COMMAND_OPTION_PARSE_ERROR;
	     } else if (var2 <= 0. && var2 != -1) {
               fprintf(stderr, "Warning: parameter var_0=%g out of range, using default\n", var2);
               errcode=PARAMETER_OUT_OF_RANGE;
             } else {
               opt_args->var[1]=var2;
             }
	     break;
      case ('w'):
             ncon=sscanf(optarg, "%f", &wmax);
             if (ncon != 1) {
               fprintf(stderr, "Warning: garbled option argument: -w %s", optarg);
               errcode=COMMAND_OPTION_PARSE_ERROR;
             } else {
               opt_args->W1=wmax;
             }
             break;
      case ('W'):
             ncon=sscanf(optarg, "%g", &Wc);
             if (ncon != 1) {
               fprintf(stderr, "Warning: garbled option argument: -W %s", optarg);
               errcode=COMMAND_OPTION_PARSE_ERROR;
	     } else if (Wc <= 0) {
               fprintf(stderr, "Warning: parameter Wc=%f out of range, using default=%f\n", Wc, opt_args->W2);
               errcode=PARAMETER_OUT_OF_RANGE;
	     } else {
	       opt_args->W2=Wc;
             }
	     break;
      case ('x'):
             opt_args->xflag=1;
	     break;
      case ('X'):
             real_a pratio;
             ncon=sscanf(optarg, "%g", &pratio);
             if (ncon != 1) {
               fprintf(stderr, "Warning: garbled option argument: -X %s", optarg);
               errcode=COMMAND_OPTION_PARSE_ERROR;
	     } else if (pratio <= 0) {
               fprintf(stderr, "Warning: parameter pratio=%f out of range, using default=%f\n", pratio, opt_args->pratio);
               errcode=PARAMETER_OUT_OF_RANGE;
	     } else {
	       opt_args->pratio=pratio;
             }
	     break;
      case ('Y'):
             opt_args->Yflag=1;
	     break;
      case ('z'):
             opt_args->zflag=1;
	     break;
      case ('Z'):
             opt_args->Zflag=1;
	     break;
      case ('?'):
             fprintf(stderr, "Unknown option: -%c -- ignored\n", optopt);
	     errcode=COMMAND_OPTION_PARSE_ERROR;
	     break;
      default:
	     fprintf(stderr, "Error parsing command line\n");
	     errcode=FATAL_COMMAND_OPTION_PARSE_ERROR;
	     break;
    }
  }

  argc-=optind;
  argv+=optind;

  return errcode;
}

//given averages and std. devs., normalize a set of vectors:
template <class real>
void norm_vec(real **mat, dim_ta D, nel_ta n, real *ave, real *std) {
  for (dim_ta j=0; j<D; j++) {
    for (nel_ta i=0; i<n; i++) mat[i][j]=(mat[i][j]-ave[j])/std[j];
  }
}

//un-normalize the vectors:
template <class real>
void unnorm_vec(real **mat, dim_ta D, nel_ta n, real *ave, real *std) {
  for (dim_ta j=0; j<D; j++) {
    for (nel_ta i=0; i<n; i++) mat[i][j]=mat[i][j]*std[j]+ave[j];
  }
}

template void norm_vec<float>(float **mat, dim_ta D, nel_ta n, float *ave, float *std);
template void norm_vec<double>(double **mat, dim_ta D, nel_ta n, double *ave, double *std);

//calculate the averages and standard deviations of a set of vectors:
template <class real>
void calc_norm(real **mat, dim_ta D, nel_ta n, real *ave, real *std) {
  real anom;		//intermediate output

  for (dim_ta i=0; i<D; i++) {
    ave[i]=0;
    for (nel_ta j=0; j<n; j++) ave[i]+=mat[j][i];
    ave[i]/=n;
    std[i]=0;
    for (nel_ta j=0; j<n; j++) {
      anom=mat[j][i]-ave[i];
      std[i]+=anom*anom;
    }
    std[i]=sqrt(std[i]/(n-1));

  }

}

template void calc_norm<float>(float **mat, dim_ta D, nel_ta n, float *ave, float *std);
template void calc_norm<double>(double **mat, dim_ta D, nel_ta n, double *ave, double *std);

//normalize a set of vectors by the std. deviations of each of the variables:
template <class real>
void norm_vec_std(real **mat1, dim_ta m, nel_ta n1, real **mat2, nel_ta n2) {
  real *ave;
  real *std;
  real anom;		//intermediate output

  //allocate space:
  ave=new real[m];
  std=new real[m];

  for (dim_ta i=0; i<m; i++) {
    ave[i]=0;
    for (nel_ta j=0; j<n1; j++) ave[i]+=mat1[j][i];
    ave[i]/=n1;
    std[i]=0;
    for (nel_ta j=0; j<n1; j++) {
      anom=mat1[j][i]-ave[i];
      std[i]+=anom*anom;
      mat1[j][i]=anom;
    }
    std[i]=sqrt(std[i]/(n1-1));

    //normalize the first matrix:
    for (nel_ta j=0; j<n1; j++) mat1[j][i]=mat1[j][i]/std[i];

    //normalize the second matrix:
    for (nel_ta j=0; j<n2; j++) mat2[j][i]=(mat2[j][i]-ave[i])/std[i];

  }

  //clean up:
  delete [] ave;
  delete [] std;

}

template void norm_vec_std<float>(float **mat1, dim_ta m, nel_ta n1, float **mat2, nel_ta n2);
template void norm_vec_std<double>(double **mat1, dim_ta m, nel_ta n1, double **mat2, nel_ta n2);

//randomizes a set of vectors:
template <class real, class cls_t>
void randomize_vec(real **mat, dim_ta m, nel_ta n, cls_t *cls) {
  int32_t *randata;
  long *index;
  gsl_rng *rann;        //GSL random number generator
  real **newmat;
  cls_t *newcls;

  //allocate space:
  randata=new int32_t[n];
  index=new long[n];
  newmat=new real *[n];
  newcls=new cls_t[n];

  //initialize random number generator:
  rann=gsl_rng_alloc(agf_gsl_rng_type);
  gsl_rng_set(rann, seed_from_clock());

  for (nel_ta i=0; i<n; i++) randata[i]=gsl_rng_get(rann);
  heapsort(randata, index, n);

  for (nel_ta i=0; i<n; i++) {
    newmat[i]=mat[index[i]];
    newcls[i]=cls[index[i]];
  }

  //its not so efficient, but it saves some trouble:
  for (nel_ta i=0; i<n; i++) {
    mat[i]=newmat[i];
    cls[i]=newcls[i];
  }

  //clean up:
  gsl_rng_free(rann);

  delete [] randata;
  delete [] index;
  delete [] newmat;
  delete [] newcls;


}

template void randomize_vec<float,cls_ta>(float **mat, dim_ta m, nel_ta n, cls_ta *cls);
template void randomize_vec<double,cls_ta>(double **mat, dim_ta m, nel_ta n, cls_ta *cls);

//sorts a set of vectors by their classes:
//
//mat is the set of vectors
//m is the dimension of each vector
//n is the number of vectors
//cl are the classes
//ncl are the number of class
//clsind is an longword array of length ncl+1 indexing into the rearranged mat
//	giving the starting point of each class
template <class real, class cls_t>
nel_ta * sort_classes(real **mat, nel_ta n, cls_t *cl, cls_t ncl) {
  real *swap;
  nel_ta *clind;
  cls_t swpcls;

  //2013-09-18 PM: class values of -1 are allowed and will be excluded
  //from the analysis; clind[ncl] should be n, not n-1
  clind=new cls_t[ncl+1];
  clind[0]=0;
  for (cls_t i=0; i<ncl; i++) {
    while (cl[clind[i]] == i-1 && clind[i] < n) {
      clind[i]++;
    }
    for (nel_ta j=clind[i]; j<n; j++) {
      if (cl[j] == i-1) {
	swap=mat[j];
        mat[j]=mat[clind[i]];
	mat[clind[i]]=swap;
	swpcls=cl[j];
	cl[j]=cl[clind[i]];
	cl[clind[i]]=swpcls;
	do {
	  clind[i]++;
	} while (cl[clind[i]] == i-1);
	if (clind[i] >= j) j=clind[i]+1;
      }
    }
    clind[i+1]=clind[i];
  }
  clind[ncl]=n;

  return clind;
}

//lets do this the right way
//it should be an O(n) operation:
//(really, how much difference will it make?
template <class real, class cls_t>
nel_ta * sort_classes_n(real **mat, nel_ta n, cls_t *cl, cls_t ncl) {
  cls_t swpcls;
  nel_ta *clind;
  nel_ta clptr[ncl];
  real **matnew;

  clind=new nel_ta[ncl+1];

  //first we simply count the number of each class label:
  clind[0]=0;
  for (nel_ta i=0; i<n; i++) clind[cl[i]+1]++;
  for (cls_t i=1; i<=ncl; i++) {
    clptr[i-1]=clind[i-1];
    clind[i]+=clind[i-1];
  }

  //need to figure out a way to do it in-place:
  matnew=new real[n];
  for (nel_ta i=0; i<n; i++) {
    matnew[clptr[cl[i]]]=mat[i];
    //clnew[clptr[cl[i]]]=cl[i];
    clptr[cl[i]]++;
  }

  //to rearrange the class labels, just go through clind variable,
  //setting labels as you go (step should be optional...):
  for (cls_ta i=0; i<=ncl; i++) {
    for (nel_ta j=clind[i]; j<clind[i+1]; j++) cl[j]=i-1;
  }

  //ruins the neat allocation scheme (but then so did the old algorithm)
  delete [] mat;
  mat=matnew;

  return clind;

}

template nel_ta * sort_classes<float,cls_ta>(float **mat, nel_ta n, cls_ta *cl, cls_ta ncl);
template nel_ta * sort_classes<double,cls_ta>(double **mat, nel_ta n, cls_ta *cl, cls_ta ncl);

template <class real, class cls_t>
real ** remove_duplicates(real **mat, dim_ta m, nel_ta n, cls_t *cls, real *wt, nel_ta &nnew) {
  long *sind;
  long ind[n];
  vector<real> *mat2;
  real wtnew[n];
  real **result;
  nel_ta lastind;

  mat2=new vector<real> *[n];

  //preparation:
  //-move data into a sortable array (including class values)
  for (nel_ta i=0; i<n; i++) {
    mat2[i]=new vector<real>(m+1);
    for (dim_ta j=0; j<m; j++) mat2[i][j]=mat[i][j];
    mat2[i][m]=cls[i];
    wtnew[i]=0;
  }

  //work, part 1: heapsort the data
  sind=heapsort(mat2, n);

  //work, part 2: group duplicates, average weights for duplicates
  nnew=1;
  lastind=0;
  wtnew[0]=wt[0];
  ind[0]=sind[0];
  for (nel_ta i=1; i<n; i++) {
    if (mat2[sind[i-1]]!=mat2[sind[i]]) {
      ind[i]=sind[i];
      nnew++;
      wtnew[i-1]/=(i-lastind);
      lastind=i-1;
    } else {
      wtnew[lastind]+=wt[i];
    }
  }
  wtnew[nnew-1]/=(nnew-lastind);

  //postprocessing: move data, including new weights, into a new array
  result=new real *[nnew];
  result[0]=new real[nnew*(m+1)];
  for (nel_ta i=0; i<nnew; i++) {
    result[i]=result[0]+i*(m+1);
    for (dim_ta j=0; j<m; j++) result[i][j]=mat[ind[i]][j];
    result[i][m]=wtnew[i];
  }

  //cleanup: delete dynamic data, return results
  delete [] mat2;

  return result;

}

} //end namespace libagf

