#include <stdio.h>
#include <getopt.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <assert.h>
#include <errno.h>

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_multifit.h>

#include "full_util.h"
//#include "peteys_tmpl_lib.h"

#include "agf_lib.h"

using namespace std;
using namespace libpetey;

namespace libagf {

  template <class real>
  real logistic_function(real x) {
    return 1-2/(1+exp(x));
  }

  //arctangent normalized to go from [-1,1] with df/dx|_x=0 = 1
  template <class real>
  real atan_norm(real x) {
    return M_PI*atan(2*x/M_PI)/2;
  }

  template <typename real, typename cls_t>
  agf2class<real, cls_t>::agf2class(real **x, cls_t *cls, nel_ta n, dim_ta nv, real v[2], real Wc, nel_ta k1) {
    cls_t c1[n];
    nel_ta *clind;

    this->order=0;
    this->calcoef=NULL;
    this->sigmoid_func=NULL;

    W=Wc;
    var0[0]=v[0];
    var0[1]=v[1];
    k=k1;
    ntrain=n;
    this->D1=nv;
    this->ncls=2;
    unsort=NULL;
    train=new real*[ntrain];
    for (nel_ta i=0; i<ntrain; i++) {
      train[i]=x[i];
      c1[i]=cls[i];
    }
    clind=sort_classes(train, ntrain, c1, this->ncls);
    cind=clind[1];

    //set variance brackets:
    set_var();

    //zero diagnostics:
    zero_diag();

    delete [] clind;
  }

  template <typename real, typename cls_t>
  agf2class<real, cls_t>::agf2class(agf_classifier<real, cls_t> *other, cls_t *part, cls_t npart) {
    cls_t ncls=1;
    this->ncls=2;
    cls_t *cls;
    nel_ta *clind;

    this->order=0;
    this->calcoef=NULL;
    this->sigmoid_func=NULL;

    W=other->W;
    var0[0]=other->var0[0];
    var0[1]=other->var0[1];
    var[0]=var0[0];
    var[1]=var0[1];
    k=other->k;

    ntrain=other->ntrain;
    this->D1=other->n_feat();
    this->ncls=2;
    this->copy_ltran(other);

    unsort=NULL;

    cls=new cls_t[ntrain];
    train=new real*[ntrain];
    for (nel_ta i=0; i<other->ntrain; i++) {
      if (other->cls[i]<0 || other->cls[i]>=npart) {
        cls[i]=-1;
      } else {
        cls[i]=part[other->cls[i]];
      }
      train[i]=other->train[i];
    }

    for (cls_t i=0; i<npart; i++) if (part[i]>ncls) ncls=part[i];

    clind=sort_classes(train, ntrain, cls, ncls);

    ntrain=ntrain-clind[0];
    trn0=train;
    train+=clind[0];

    cind=clind[1];

    //set variance brackets:
    set_var();

    //zero diagnostics:
    zero_diag();

    delete [] cls;
    delete [] clind;
  }

  template <typename real, typename cls_t>
  agf2class<real, cls_t>::~agf2class() {
    delete [] trn0;
    if (unsort!=NULL) delete [] unsort;
  }

  template <typename real, typename cls_t>
  void agf2class<real, cls_t>::set_var() {
    //check parameter ranges:
    if (var0[0] <= 0 || var0[1] <= 0) {
      //calculate the averages and standard deviations:
      real std[this->D1];
      real ave[this->D1];
      real vart;

      calc_norm(this->train, this->D1, this->ntrain, ave, std);

      //if the initial filter variance is not set, set it to the total
      //variance of the data:
      //variance brackets are then fully "sticky..."
      vart=0;
      for (dim_ta i=0; i<this->D1; i++) vart+=std[i]*std[i];
      if (var0[0] <= 0) {
        var[0]=vart/pow(this->ntrain, 2./this->D1);
        fprintf(stderr, "Using %10.3g for lower filter variance bracket\n\n", var[0]);
      }
      if (var0[1] <= 0) {
        var[1]=vart;
        fprintf(stderr, "Using %10.3g for upper filter variance bracket\n\n", var[1]);
      }
    }
  }

  template <typename real, typename cls_t>
  void agf2class<real, cls_t>::zero_diag() {
    //zero diagnostics:
    min_f=-1;
    max_f=0;
    total_f=0;
    min_W=-1;
    max_W=0;
    total_W=0;
    min_nd=-1;
    max_nd=0;
    total_nd=0;
    ntrial=0;
    ntrial_k=0;
  }

  template <class real, class cls_t>
  real agf2class<real, cls_t>::R(real *x, real *praw) {
    agf_diag_param diag;
    real r;
    if (k>0) {
      r=dcalc(train, this->D1, ntrain, cind, x, var, k, W, &diag);
      if (ntrial_k==0) min_f=diag.f;
      if (diag.f < min_f) min_f=diag.f;
      else if (diag.f > max_f) max_f=diag.f;
      total_f+=diag.f;
      ntrial_k++;
    } else {
      r=dcalc(train, this->D1, ntrain, cind, x, var, W, &diag);
    }
    if (ntrial==0) {
      min_W=diag.W;
      min_nd=diag.nd;
    }

    //calculate diagnostics:
    if (diag.nd < min_nd) min_nd=diag.nd;
    else if (diag.nd > max_nd) max_nd=diag.nd;
    total_nd+=diag.nd;

    if (diag.W < min_W) min_W=diag.W;
    else if (diag.W > max_W) max_W=diag.W;
    total_W+=diag.W;
    ntrial++;

    if (praw!=NULL) {
      *praw=r;
    }

    return r;
  }

  template <class real, class cls_t>
  real agf2class<real, cls_t>::R_deriv(real *x, real *drdx) {
    agf_diag_param diag;
    real r;
    if (k>0) {
      r=dgrad(train, this->D1, ntrain, cind, x, var, k, W, drdx, &diag);
      if (ntrial_k==0) min_f=diag.f;
      if (diag.f < min_f) min_f=diag.f;
      else if (diag.f > max_f) max_f=diag.f;
      total_f+=diag.f;
      ntrial_k++;
    } else {
      r=dgrad(train, this->D1, ntrain, cind, x, var, W, drdx, &diag);
    }
    if (ntrial==0) {
      min_W=diag.W;
      min_nd=diag.nd;
    }

    //calculate diagnostics:
    if (diag.nd < min_nd) min_nd=diag.nd;
    else if (diag.nd > max_nd) max_nd=diag.nd;
    total_nd+=diag.nd;

    if (diag.W < min_W) min_W=diag.W;
    else if (diag.W > max_W) max_W=diag.W;
    total_W+=diag.W;
    ntrial++;

    return r;
  }

  template <class real, class cls_t>
  void agf2class<real, cls_t>::reset_var() {
    //if we don't want variance brackets to be completely "sticky":
    if (var0[0]>0) var[0]=var0[0];
    if (var0[1]>0) var[1]=var0[1];
  }

  template <class real, class cls_t>
  int agf2class<real, cls_t>::param_init(nel_ta ns,
		  int (*sfunc)(void *, real *, real *),
		  bordparam<real> *param) {
    int err;
    reset_var();
    if ((2.*ns)/ntrain/cind > 0.25) {
      //for small datasets:
      err=agfbordparam_init(param, train, this->D1, ntrain, cind, var, k, W, 1);
      sfunc=&oppositesample_small<real>;
    } else {
      //for large datasets:
      err=agfbordparam_init(param, train, this->D1, ntrain, cind, var, k, W);
      sfunc=&oppositesample<real>;
    }
    return err;
      
  }

  template <class real, class cls_t>
  borders_classifier<real, cls_t>::borders_classifier() {
    brd=NULL;
    grd=NULL;
    n=0;
    gd=NULL;
    //since we don't trust C++ to call the parent initializer:
    this->ncls=2;
    this->D=0;
    this->D1=0;
    this->sigmoid_func=NULL;
  }

  template <class real, class cls_t>
  borders_classifier<real, cls_t>::borders_classifier(const char *fbase, int sigtype) {
    int err;
    SIGFUN_TYPE (*sigfun) (SIGFUN_TYPE);
    switch (sigtype) {
      case (-1):
        sigfun=NULL;
	break;
      case (0):
        sigfun=&tanh;
	break;
      case (1):
        sigfun=&erf;
	break;
      case (2):
	sigfun=&logistic_function;
	break;
      case (3):
	sigfun=&atan_norm;
	break;
      default:
	fprintf(stderr, "borders_classifier: code for function to transform decision values (%d) not recognized\n", sigtype);
	fprintf(stderr, "             [0=tanh; 1=erf; 2=2/(1-exp(..))]--using tanh()\n");
        sigfun=&tanh;
    }
    err=init(fbase, sigfun);
    if (err!=0) throw err;
  }

  template <class real, class cls_t>
  borders_classifier<real, cls_t>::borders_classifier(const char *fbase, 
		SIGFUN_TYPE (* sigfun)(SIGFUN_TYPE)) {
    int err=init(fbase, sigfun);
    if (err!=0) throw err;
  }

  template <class real, class cls_t>
  int borders_classifier<real, cls_t>::init(const char *fbase, 
	 	SIGFUN_TYPE (* sigfun)(SIGFUN_TYPE)) {
    int err;

    this->mat=NULL;
    this->id=-1;

    this->sigmoid_func=sigfun;

    this->name=new char[strlen(fbase)+1];
    strcpy(this->name, fbase);

    err=agf_read_borders(fbase, brd, grd, n, this->D);
    if (err!=0) return err;

    fprintf(stderr, "borders_classifier: %d border samples found in model, %s\n", n, fbase);
    this->D1=this->D;
    this->ncls=2;

    //(don't) calculate and store the lengths of all the gradient vectors for possible
    //use later on: (we'll get them when we need them...)
    gd=NULL;

    this->calcoef=new real[2];
    this->calcoef[0]=0;
    this->calcoef[1]=1;
    this->order=1;

    return err;
  }

  template <class real, class cls_t>
  void borders_classifier<real, cls_t>::train(agf2class<real, cls_t> *agf,
		  nel_ta ns, real tol) {
    bordparam<real> param;
    int (*sfunc) (void *, real *, real *);		//sampling function
    real **xtran;

    if (grd!=NULL) delete [] grd;
    if (brd!=NULL) delete [] brd;
    if (gd!=NULL) delete [] gd;

    this->sigmoid_func=&tanh;
    //transfer parameters from agf to this:
    this->ncls=agf->n_class();
    this->D1=agf->n_feat();
    this->copy_ltran(agf);

    //initialize parameters for borders training:
    //(it's kind of weird because the "agf" object instantiation isn't used
    //at all in the process... it just passes it's parameters to the somewhat
    //archaic C-style procedural algorithms)
    agf->param_init(ns, sfunc, &param);

    //allocate the arrays for holding the results:
    brd=allocate_matrix<real, nel_ta>(ns, this->D1);
    grd=allocate_matrix<real, nel_ta>(ns, this->D1);

    //find class borders:
    n=sample_class_borders(&agfrfunc<real>, sfunc, &param, 
		ns, this->D1, tol, agf_global_borders_maxiter, 
		brd, grd);
    //maybe convert to linux system errors:
    if (n<=0) throw ENODATA;

    //clean up:
    agfbordparam_clean(&param);

    gd=NULL;
  }

  //train borders classifier from a SVM model:
  template <class real, class cls_t>
  void borders_classifier<real, cls_t>::train(svm2class<real, cls_t> *svm,
		  real **x, cls_t *cls, dim_ta nvar, nel_ta ntrain,
		  nel_ta ns, real tol, int tflag) {
    bordparam<real> param;
    real *xsort[ntrain];		//sort training data
    cls_t csel[ntrain];
    nel_ta *clind;			//indices for sorted classes
    nel_ta ntrain2;
    int (*sfunc) (void *, real *, real *);		//sampling function
    real **xtran;

    if (grd!=NULL) delete [] grd;
    if (brd!=NULL) delete [] brd;
    if (gd!=NULL) delete [] gd;

    this->sigmoid_func=&tanh;
    //copy parameters:
    this->ncls=svm->n_class();
    this->D1=svm->n_feat();
    //copy linear coord. trans.:
    if (tflag && this->copy_ltran(svm)) {
      //is this necessary? shouldn't the training data be transformed already?
      if (this->D!=nvar) throw DIMENSION_MISMATCH;
      xtran=allocate_matrix<real, int32_t>(ntrain, this->D1);
      for (nel_ta i=0; i<ntrain; i++) {
        for (dim_ta j=0; j<this->D1; j++) {
          xtran[i][j]=0;
          for (dim_ta k=0; k<this->D; k++) {
            real diff=x[i][k]-this->b[k];
            xtran[i][j]+=diff*this->mat[k][j];
	  }
        }
      }
    } else {
      if (this->D1!=nvar) throw DIMENSION_MISMATCH;
      this->D=this->D1;
      xtran=x;
    }

    //remove side effects:
    for (nel_ta k=0; k<ntrain; k++) {
      csel[k]=cls[k];
      xsort[k]=xtran[k];
    }
    //sort the classes:
    clind=sort_classes(xsort, ntrain, csel, this->ncls);
    ntrain2=clind[2]-clind[0];

    //initialize parameters:
    if (ns/(ntrain2-clind[1])/clind[1] > 0.25) {
      //for small datasets:
      bordparam_init(&param, xsort+clind[0], nvar, ntrain2, clind[1]-clind[0], 1);
      sfunc=&oppositesample_small<real>;
    } else {
      //for large datasets:
      bordparam_init(&param, xsort+clind[0], nvar, ntrain2, clind[1]-clind[0]);
      sfunc=&oppositesample<real>;
    }
    param.rparam=svm;

    //allocate the arrays for holding the results:
    brd=allocate_matrix<real, nel_ta>(ns, nvar);
    grd=allocate_matrix<real, nel_ta>(ns, nvar);

    //find class borders:
    n=sample_class_borders(&svmrfunc<real, cls_t>, sfunc, &param, 
		ns, this->D1, tol, agf_global_borders_maxiter, 
		brd, grd);
    if (n<=0) throw ENODATA;

    delete [] clind;
    if (tflag && this->mat!=NULL) {
      delete [] xtran[0];
      delete [] xtran;
    }
    bordparam_clean(&param);

    gd=NULL;
  }

  template <class real, class cls_t>
  borders_classifier<real, cls_t>::~borders_classifier() {
    delete_matrix(brd);
    delete_matrix(grd);
    if (gd!=NULL) delete [] gd;
  }

  template <class real, class cls_t>
  void borders_classifier<real, cls_t>::calc_grad_len() {
    //calculate and store the lengths of all the gradient vectors for possible
    //use later on:
    gd=new real[n];
    for (nel_ta i=0; i<n; i++) {
      gd[i]=0;
      for (dim_ta j=0; j<this->D; j++) gd[i]+=grd[i][j]*grd[i][j];
      gd[i]=sqrt(gd[i]);
    }
  }

  //transform border vectors and gradients: (are stored in original coords)
  template <class real, class cls_t>
  int borders_classifier<real, cls_t>::ltran_model(real **mat1, real *b1, dim_ta d1, dim_ta d2) {
    int err2=0;

    real **brd2;
    real **grd2;

    gsl_matrix *u;
    gsl_vector *s;
    gsl_matrix *vt;
    gsl_vector *work;

    if (d1!=this->D1) {
      fprintf(stderr, "borders_classifier: first dimension (%d) of trans. mat. does not agree with that of borders data (%d)\n", d1, this->D1);
      return DIMENSION_MISMATCH;
    }
    fprintf(stderr, "borders_classifier: Normalising the border samples...\n");

    //apply constant factor:
    for (nel_ta i=0; i<n; i++) {
      for (dim_ta j=0; j<d1; j++) {
        brd[i][j]=brd[i][j]-b1[j];
      }
    }

    brd2=matrix_mult(brd, mat1, n, d1, d2);
    grd2=allocate_matrix<real, int32_t>(n, d2);

    //gradients do NOT transform the same as the vectors:
    u=gsl_matrix_alloc(d1, d2);
    for (dim_ta i=0; i<d1; i++) {
      for (dim_ta j=0; j<d2; j++) {
        gsl_matrix_set(u, i, j, mat1[i][j]);
      }
    }
    s=gsl_vector_alloc(d2);
    vt=gsl_matrix_alloc(d2, d2);
    work=gsl_vector_alloc(d2);
    gsl_linalg_SV_decomp(u, vt, s, work);
    gsl_vector_free(work);
    for (nel_ta i=0; i<n; i++) {
      double tmp_g;
      for (dim_ta j=0; j<d2; j++) {
        grd2[i][j]=0;
        for (dim_ta k=0; k<d2; k++) {
          double vt_el;
          double s_k=gsl_vector_get(s, k);
          if (s_k == 0) continue;
          tmp_g=0;
          for (dim_ta l=0; l<d1; l++) {
            double u_el=gsl_matrix_get(u, l, k);
            tmp_g+=u_el*grd[i][l];
          }
          vt_el=gsl_matrix_get(vt, j, k);
          grd2[i][j]+=tmp_g*vt_el/s_k;
        }
      }
    }

    if (this->mat == NULL) {
      //the classifier now has a different number of features:
      this->D1=d2;
      this->D=d2;
    } else {
      //assert(mat1==this->mat && b1==this->b);
      //from the outside, the classifier looks like it has the same number of
      //features as before normalization:
      assert(this->D1==d2);		//(doesn't really make sense given test above...)
    }

    gsl_matrix_free(u);
    gsl_vector_free(s);
    gsl_matrix_free(vt);

    delete_matrix(brd);
    delete_matrix(grd);
    brd=brd2;
    grd=grd2;

    //have to recalculate all the gradient lengths, dammit:
    if (gd!=NULL) {
      delete [] gd;
      calc_grad_len();
    }

    return err2;
  }

  template <class real, class cls_t>
  real borders_classifier<real, cls_t>::decision(real *x) {
    nel_ta k;		//index of nearest border sample
    real d;		//distance from border sample 
    			//(may be useful for continuum generalization)
    return border_classify0(brd, grd, this->D1, n, x, k, d);
  }

  template <class real, class cls_t>
  void borders_classifier<real, cls_t>::print(FILE *fs, char *fbase, int depth) {
    for (int i=0; i<depth; i++) fprintf(fs, "  ");
    if (fbase==NULL) {
      fprintf(fs, "%s", this->name);
    } else {
      fprintf(fs, "%s", fbase);
    }
  }

  template <class real, class cls_t>
  int borders_classifier<real, cls_t>::load(FILE *fs, int vflag) {
    int32_t n1, n2;
    int32_t nvar1, nvar2;
    int err=0;
    brd=scan_matrix<real, int32_t>(fs, n1, nvar1);
    grd=scan_matrix<real, int32_t>(fs, n2, nvar2);
    if (brd==NULL || grd==NULL) err=FILE_READ_ERROR;
    if (n1!=n2 || nvar1!=nvar2) err=SAMPLE_COUNT_MISMATCH;
    this->D1=nvar1;
    n=n1;
    this->D=this->D1;
    if (vflag) this->sigmoid_func=NULL; else this->sigmoid_func=&tanh;
    gd=NULL;
    return err;
  }

  template <class real, class cls_t>
  int borders_classifier<real, cls_t>::load(FILE *fs) {
    return load(fs, 0);
  }

  template <class real, class cls_t>
  int borders_classifier<real, cls_t>::save(FILE *fs) {
    print_matrix(fs, brd, n, this->D1);
    print_matrix(fs, grd, n, this->D1);
    return 0;
  }

  template <class real, class cls_t>
  borders_calibrated<real, cls_t>::borders_calibrated() {
    this->brd=NULL;
    this->grd=NULL;
    this->n=0;
    this->gd=NULL;
    //since we don't trust C++ to call the parent initializer:
    this->ncls=2;
    this->D=0;
    this->D1=0;
    this->sigmoid_func=&tanh;
    this->order=1;
    this->calcoef=new real[2];
    this->calcoef[0]=0;
    this->calcoef[1]=1;
  }

  template <class real, class cls_t>
  borders_calibrated<real, cls_t>::borders_calibrated(const char *fbase) {
    FILE *fs;
    real **mat;
    int32_t m, n;
    char *fname=new char[strlen(fbase)+5];
    int err=this->init(fbase, &tanh);
    if (err!=0) throw err;
    this->calcoef=new real[2];
    this->calcoef[0]=0;
    this->calcoef[1]=1;
    this->order=1;
    sprintf(fname, "%s.clb", fbase);
    fs=fopen(fname, "r");
    if (fs==NULL) {
      fprintf(stderr, "borders_calibrated: calibration file, %s, not found.\n", fname);
      delete [] fname;
      return;
    }
    mat=scan_matrix<real, cls_t>(fs, m, n);
    if (m!=1 || mat==NULL) {
      fprintf(stderr, "borders_calibrated: error reading calibration file, %s .\n", fname);
      delete [] fname;
      fclose(fs);
      return;
    }
    delete [] this->calcoef;
    this->order=n-1;
    this->calcoef=mat[0];
    delete [] mat;
    fclose(fs);
    delete [] fname;
  }

  template <class real, class cls_t>
  borders_calibrated<real, cls_t>::~borders_calibrated() {
  }

  template <class real, class cls_t>
  void borders_calibrated<real, cls_t>::train(real **train, cls_t *cls, nel_ta ntrain, int type, real *param) {
    calibrate(train, cls, ntrain, param[0], param[1]);
  }

  template <class real, class cls_t>
  void borders_calibrated<real, cls_t>::calibrate(real **train, cls_t *cls, nel_ta ntrain, int O, int nhist, int optimize, real r0) {
    real **tab;			//table of accuracies versus probabilities
    int nfit;
    cls_t nct, ncr;		//number of classes (should be 2...)
    real r[ntrain];		//difference in conditional probabilities
    int nbad=0;			//remove bad values
    gsl_matrix *A;		//A^T A x = A^T b		
    gsl_vector *x;
    gsl_vector *b;
    //for the fitting:
    gsl_matrix *cov;
    gsl_multifit_linear_workspace *work;
    double chisq;
    real *rsort;
    nel_ta *nonecum;
    real cterm=0;			//constant term
    int cflag=0;		//constant term supplied?

    //classify each training sample using this classifier:
    printf("border_calibrated: classifiying training data\n");
    for (int i=0; i<ntrain; i++) {
      r[i]=this->R_t(train[i]);
      //r[i]=this->borders_classifier<real, cls_t>::R_t(train[i]);
    }

    //build a histogram of probabilities:
    printf("border_calibrated: building histogram of accuracies\n");
    tab=con_acc_table2(cls, r, ntrain, nhist);

    for (int i=0; i<2*nhist; i++) {
      printf("%g %g\n", tab[0][i], tab[1][i]);
      tab[0][i]=atanh(tab[0][i]);
      tab[1][i]=atanh(tab[1][i]);
    }

    //remove bad values and prepare for fitting:
    printf("border_calibrated: removing bad values\n");
    for (int i=0; i<2*nhist; i++) {
      if (isfinite(tab[1][i])!=1) {
        nbad++;
      } else {
        tab[0][i-nbad]=tab[0][i];
        tab[1][i-nbad]=tab[1][i];
      }
    }
    nfit=2*nhist-nbad;

    printf("border_calibrated: optimizing skill\n");
    //the coup de grace:
    if (optimize>1) {
      cflag=1;
      rsort=new real[ntrain];
      nonecum=new nel_ta[ntrain];
      sortr_cumulate_ones(cls, r, ntrain, nonecum, rsort);
    }
    switch (optimize) {
      case(1):
        cflag=1;
	cterm=-atanh(r0);
	break;
      case(2):
	cterm=-atanh(optimize_binary_skill_rig(nonecum, rsort, ntrain, &acc_bin));
	break;
      case(3):
	cterm=-atanh(optimize_binary_skill_rig(nonecum, rsort, ntrain, &uc_bin));
	break;
      case(4):
	cterm=-atanh(optimize_binary_skill_rig(nonecum, rsort, ntrain, &corr_bin));
	break;
    }
    if (optimize>1) {
      delete [] rsort;
      delete [] nonecum;
    }

    //allocate vectors and matrices that define the fitting problem:
    //best fit for A x = b
    //where a_ij = atanh^j(r_i), tanh(b_i) is actual accuracy and x are the fitting coefficients
    printf("border_calibrated: preparing for fitting\n");
    this->order=O;
    A=gsl_matrix_alloc(nfit, this->order+1-cflag);
    x=gsl_vector_alloc(this->order+1-cflag);
    b=gsl_vector_alloc(nfit);

    //fill vectors and matrices:
    for (int i=0; i<nfit; i++) {
      printf("%g %g\n", tab[0][i], tab[1][i]);
      for (int j=cflag; j<=this->order; j++) {
        //here we don't care so much about efficiency because we only do the
	//fitting once and it's a small problem...
        gsl_matrix_set(A, i, j-cflag, pow(tab[0][i]+cterm, j));
      }
      gsl_vector_set(b, i, tab[1][i]);
    }

    //allocate extra stuff we need to perform the fitting:
    work=gsl_multifit_linear_alloc(nfit, this->order+1-cflag);
    cov=gsl_matrix_alloc(this->order+1-cflag, this->order+1-cflag);

    //perform fitting:
    printf("border_calibrated: performing fitting\n");
    gsl_multifit_linear(A, b, x, cov, &chisq, work);

    //pull coefficients out from GSL vector type:
    printf("border_calibrated: storing coefficients\n");
    delete [] this->calcoef;
    this->calcoef=new real[this->order+1];
    for (int i=cflag; i<=this->order; i++) this->calcoef[i]=gsl_vector_get(x, i-cflag);
    if (cflag) this->calcoef[0]=cterm;

    //clean up:
    printf("border_calibrated: cleaning up\n");
    delete [] tab[0];
    delete [] tab;

    gsl_multifit_linear_free(work);
    gsl_matrix_free(A);
    gsl_matrix_free(cov);
    gsl_vector_free(x);
    gsl_vector_free(b);
  }

  template <class real, class cls_t>
  void borders_calibrated<real, cls_t>::calibrate2(real **train, cls_t *cls, nel_ta ntrain, int O, int nhist, int optimize, real r0) {
    real r[ntrain];		//difference in conditional probabilities
    real *rsort;
    nel_ta *nonecum;

    this->order=O;
    //classify each training sample using this classifier:
    printf("borders_calibrated: classifiying training data\n");
    for (int i=0; i<ntrain; i++) {
      r[i]=this->R_t(train[i]);
      //r[i]=this->borders_classifier<real, cls_t>::R_t(train[i]);
    }

    //the coup de grace:
    if (optimize>1) {
      rsort=new real[ntrain];
      nonecum=new nel_ta[ntrain];
      sortr_cumulate_ones(cls, r, ntrain, nonecum, rsort);
    } else if (optimize==0) {
      r0=2;
    }
    switch (optimize) {
      case(2):
	r0=optimize_binary_skill_rig(nonecum, rsort, ntrain, &acc_bin);
	break;
      case(3):
	r0=optimize_binary_skill_rig(nonecum, rsort, ntrain, &uc_bin);
	break;
      case(4):
	r0=optimize_binary_skill_rig(nonecum, rsort, ntrain, &corr_bin);
	break;
    }
    if (optimize>1) {
      delete [] rsort;		//this is wasted!
      delete [] nonecum;
    }

    printf("borders_calibrated: calling calibrate_r\n");
    delete [] this->calcoef;
    this->calcoef=calibrate_r(cls, r, ntrain, this->order, nhist, r0);
  }

  template <class real, class cls_t>
  void borders_calibrated<real, cls_t>::print_calib(FILE *fs) {
    print_matrix(fs, &this->calcoef, 1, this->order+1);
  }

  /*
  template <class real, class cls_t>
  real borders_calibrated<real, cls_t>::R(real *x, real *praw) {
    real r;
    nel_ta k;		//intermediate values in the calculation
    real tr=0;		//transformed decision value
    real rpi=1;		//r to power i
    real d;		//may be useful for continuum generalization

    r=border_classify0(this->brd, this->grd, this->D1, this->n, x, k, d);
    if (this->id>=0 && praw!=NULL) {
      praw[this->id]=r;
    }

    for (int i=0; i<=order; i++) {
      tr+=this->calcoef[i]*rpi;
      rpi*=r;
    }
      
    return (*this->sigmoid_func) (tr);
  }
  */

  template float logistic_function<float>(float);
  template double logistic_function<double>(double);

  template float atan_norm<float>(float);
  template double atan_norm<double>(double);

  template class agf2class<float, cls_ta>;
  template class agf2class<double, cls_ta>;

  template class borders_classifier<float, cls_ta>;
  template class borders_classifier<double, cls_ta>;

  template class borders_calibrated<float, cls_ta>;
  template class borders_calibrated<double, cls_ta>;

}

