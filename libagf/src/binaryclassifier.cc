#include <stdio.h>
#include <getopt.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <assert.h>

#include <gsl/gsl_linalg.h>

#include "full_util.h"
//#include "peteys_tmpl_lib.h"

#include "agf_lib.h"

using namespace std;
using namespace libpetey;

namespace libagf {

  template <class real, class cls_t>
  binaryclassifier<real, cls_t>::binaryclassifier() {
    this->name=NULL;
    id=-1;
    this->ncls=2;
    //xtran=NULL;
  }

  template <class real, class cls_t>
  binaryclassifier<real, cls_t>::binaryclassifier(char *nm) {
    this->name=new char[strlen(nm)+1];
    strcpy(this->name, nm);
    id=-1;
    this->ncls=2;
    //xtran=NULL;
  }

  template <class real, class cls_t>
  binaryclassifier<real, cls_t>::~binaryclassifier() {
    delete [] this->name;
  }

  template <class real, class cls_t>
  real binaryclassifier<real, cls_t>::R(real *x, real *praw) {
    //want to make this do a direct classifications using the options listed
    //in name!
  }

  template <class real, class cls_t>
  real binaryclassifier<real, cls_t>::R_t(real *x, real *praw) {
    real r;
    real *xtran;
    xtran=this->do_xtran(x);
    r=R(xtran, praw);
    if (this->mat!=NULL) delete [] xtran;
    return r;
  }

  //next five function definitions are defaults which rely on R being defined:
  template <class real, class cls_t>
  cls_t binaryclassifier<real, cls_t>::classify(real *x, real &p, real *praw) {
    cls_t result;
    real r=R(x, praw);
    return convertR(r, p);
  }

  template <class real, class cls_t>
  cls_t binaryclassifier<real, cls_t>::classify(real *x, real *p, real *praw) {
    cls_t result;
    real r=R(x, praw);
    return convertR(r, p);
  }

  template <class real, class cls_t>
  void binaryclassifier<real, cls_t>::batchR(real **x, real *r, nel_ta n, dim_ta nvar) {
    for (nel_ta i=0; i<n; i++) r[i]=R(x[i]);
  }

  template <class real, class cls_t>
  void binaryclassifier<real, cls_t>::batchR_t(real **x, real *r, nel_ta n, dim_ta nvar) {
    real *xtran[n];
    for (nel_ta i=0; i<n; i++) xtran[i]=this->do_xtran(x[i]);
    batchR(xtran, r, n, nvar);
    if (this->mat!=NULL) for (nel_ta i=0; i<n; i++) delete [] xtran[i];
  }

  template <class real, class cls_t>
  void binaryclassifier<real, cls_t>::batch_classify(real **x, cls_t *c, real *p, nel_ta n, dim_ta nvar) {
    real r[n];
    batchR(x, r, n, nvar);
    for (nel_ta i=0; i<n; i++) {
      c[i]=convertR(r[i], p[i]);
    }
  }

  template <class real, class cls_t>
  void binaryclassifier<real, cls_t>::batch_classify(real **x, cls_t *c, real **p, nel_ta n, dim_ta nvar) {
    real r[n];
    batchR(x, r, n, nvar);
    for (nel_ta i=0; i<n; i++) {
      c[i]=convertR(r[i], p[i]);
    }
  }

  template <class real, class cls_t>
  void binaryclassifier<real, cls_t>::print(FILE *fs, char *fbase, int depth) {
    for (int i=0; i<depth; i++) fprintf(fs, "  ");
    if (fbase==NULL) {
      fprintf(fs, "\"%s\"", this->name);
    } else {
      fprintf(fs, "%s", fbase);
    }
  }

  template <class real, class cls_t>
  int binaryclassifier<real, cls_t>::commands(multi_train_param &param, 
		cls_t **clist, char *fbase) {
    char *tmpname;
    if (param.partcom==NULL) {
      //print command name etc.:
      fprintf(param.commandfs, "%s %s %s %s", param.commandname, this->name,
		param.train, fbase);
    } else {
      //if a command has been set to partition class labels, run that first to generate
      //a temporary file--use the temporary file for training...
      tmpname=new char [strlen(fbase)+30];
      sprintf(tmpname, "%s.%u.tmp", fbase, param.session_id);
      fprintf(param.commandfs, "%s %s", param.partcom, param.train);
      if (param.concom==NULL) {
        fprintf(param.commandfs, " %s", tmpname);
      }
    }

    //print class partitions:
    for (cls_t i=0; clist[0]+i!=clist[1]; i++) {
      fprintf(param.commandfs, " %d", clist[0][i]);
    }
    fprintf(param.commandfs, " %c", PARTITION_SYMBOL);
    for (cls_t i=0; clist[1]+i!=clist[2]; i++) {
      fprintf(param.commandfs, " %d", clist[1][i]);
    }

    //extra features add a lot of clutter!
    if (param.partcom!=NULL) {
      //if there is a command for file conversion, pipe to that first, then to the temp. file:
      if (param.concom!=NULL) {
        fprintf(param.commandfs, " | %s > %s\n", param.concom, tmpname);
      } else {
        fprintf(param.commandfs, "\n");
      }
      //if applicable, run binary classifier on temporary and delete the temporary file;
      fprintf(param.commandfs, "%s %s %s %s\n", param.commandname, this->name, tmpname, fbase);
      if (param.Kflag==0) fprintf(param.commandfs, "rm -f %s\n", tmpname);
      delete [] tmpname;
    } else {
      fprintf(param.commandfs, "\n");
    }


    return 2;
  }

  template <class real, class cls_t>
  void binaryclassifier<real, cls_t>::set_id(cls_t *id1) {
    id=*id1;
    (*id1)++;
  }

  template <class real, class cls_t>
  general2class<real, cls_t>::general2class(const char *mod, const char *com) {
    FILE *fs;

    command=new char [strlen(com)+1];
    strcpy(command, com);
    this->name=new char [strlen(mod)+1];
    strcpy(this->name, mod);
    this->id=-1;

    syscall=new char [strlen(command)+strlen(mod)+5];
    sprintf(syscall, "%s -N %s", command, mod);

    printf("%s\n", syscall);
    fs=popen(syscall, "r");
    fscanf(fs, "%d", &this->D);
    pclose(fs);

    printf("model, %s, has %d variables\n", mod, this->D);

    delete [] syscall;
    insert_pt=strlen(command)+strlen(this->name)+1;
    syscall=new char [insert_pt+this->D*14+40];
    sprintf(syscall, "%s %s", command, this->name);

    this->mat=NULL;
    //this->xtran=NULL;
    this->ncls=2;
    this->D1=this->D;

    Kflag=0;
    Mflag=0;
  }

  template <class real, class cls_t>
  general2class<real, cls_t>::general2class(const char *mod, const char *com, int Mf, int kf) {
    FILE *fs;

    Mflag=Mf;
    Kflag=kf;
    this->id=-1;

    command=new char [strlen(com)+1];
    strcpy(command, com);
    this->name=new char [strlen(mod)+1];
    strcpy(this->name, mod);

    insert_pt=strlen(command)+strlen(this->name)+1;
    syscall=new char [insert_pt+40];
    sprintf(syscall, "%s %s", command, this->name);

    this->mat=NULL;
    //this->xtran=NULL;
    this->ncls=2;
  }

  template <class real, class cls_t>
  general2class<real, cls_t>::~general2class() {
    delete [] command;
    //workspace variable (*** need to eliminate ***):
    delete [] syscall;
    //delete [] this->name;
  }

  template <class real, class cls_t>
  real general2class<real, cls_t>::R(real *x, real *praw) {
    real r;
    batchR(&x, &r, 1, this->D1);
    if (praw!=NULL && this->id>=0) praw[this->id]=r;
    return r;
  }

  template <class real, class cls_t>
  void general2class<real, cls_t>::batchR(real **x, real *r, nel_ta n, dim_ta nvar) {
    cls_t cls[n];
    real **p;

    if (n==0) return;

    p=new real *[n];
    p[0]=new real[n*2];
    for (int i=1; i<n; i++) p[i]=p[0]+i*2;

    batch_classify(x, cls, p, n, nvar);
    for (nel_ta i=0; i<n; i++) r[i]=p[i][1]-p[i][0];

    delete [] p[0];
    delete [] p;
  }

  template <class real, class cls_t>
  void general2class<real, cls_t>::batch_classify(real **x, cls_t *cls, real *p, nel_ta n, dim_ta nvar) {
    real **p2;

    if (n==0) return;

    p2=new real *[n];
    p2[0]=new real[n*2];
    for (int i=1; i<n; i++) p2[i]=p2[0]+i*2;

    batch_classify(x, cls, p2, n, nvar);
    for (nel_ta i=0; i<n; i++) {
      if (cls[i]<=1 && cls[i]>=0) p[i]=p2[i][cls[i]]; else p[i]=0;
    }

    delete [] p2[0];
    delete [] p2;
  }

  template <class real, class cls_t>
  void general2class<real, cls_t>::batch_classify(real **x, cls_t *cls, real **p, nel_ta n, dim_ta nvar) {
    FILE *fs;
    unsigned long uid;
    cls_t ncls;
    int err;

    char infile[25];
    char outfile[25];

    char format[20];
    char fcode[4];

    if (n==0) return;

    //need different format strings depending upon what floating types
    //we're using:
    get_format_code<real>(fcode);

    uid=seed_from_clock();

    //write directly to the file to save a little space:
    sprintf(infile, "in.%u.tmp", uid);
    fs=fopen(infile, "w");
    if (Mflag==0) fprintf(fs, "%d\n", nvar);
    for (nel_ta i=0; i<n; i++) {
      cls[i]=0;			//clear class data
      //we write to the input file:
      if (Mflag) {
        sprintf(format, " %%d:%%.12%s", fcode);
        fprintf(fs, "%d", cls[i]);
        for (dim_ta j=0; j<nvar; j++) fprintf(fs, format, j+1, x[j]);
        fprintf(fs, "\n");
      } else {
        sprintf(format, "%%.12%s ", fcode);
        for (dim_ta j=0; j<nvar; j++) fprintf(fs, format, x[j]);
        fprintf(fs, "%d\n", cls[i]);
      }
    }
    fclose(fs);

    sprintf(outfile, "out.%u.tmp", uid);
    sprintf(syscall, "%s %s %s %s", command, infile, this->name, outfile);

    printf("%s\n", syscall);
    err=system(syscall);
    if (err!=0) {
      fprintf(stderr, "general2class->batch_classify: external command,\n  %s\nreturned error code, %d, exiting\n", syscall, err);
      exit(err);
    }

    //delete input file:
    if (Kflag==0) {
      remove(infile);
      //sprintf(syscall, "rm %s", infile);
      //printf("%s\n", syscall);
      //system(syscall);
    }

    //now we read from the output file:
    //printf("Reading file, %s\n", outfile);
    fs=fopen(outfile, "r");
    if (fs==NULL) {
      fprintf(stderr, "general2class->batch_classify: error opening output file, exiting\n");
      exit(FILE_READ_ERROR);
    }
    err=read_svmout(fs, cls, p, ncls, n);
    fclose(fs);

    //delete output file:
    if (Kflag==0) {
      remove(outfile);
      //sprintf(syscall, "rm %s", outfile);
      //printf("%s\n", syscall);
      //system(syscall);
    }

    if (ncls!=2) {
      fprintf(stderr, "general2class: found %d classes in output file, expecting 2\n", ncls);
      exit(PARAMETER_OUT_OF_RANGE);
    }

    if (err!=n) {
      fprintf(stderr, "general2class: error reading output file, %s, exiting\n", outfile);
      exit(FILE_READ_ERROR);
    }
      
  }

  template <class real, class cls_t>
  void general2class<real, cls_t>::print(FILE *fs, char *fbase, int depth) {
    for (int i=0; i<depth; i++) fprintf(fs, "  ");
    if (fbase==NULL) {
      fprintf(fs, "%s", this->name);
    } else {
      fprintf(fs, "%s", fbase);
    }
  }

  template <class real, class cls_t>
  int general2class<real, cls_t>::commands(multi_train_param &param, 
		cls_t **clist, char *fbase) {
    //accelerator mode:
    if (strcmp(command, "")==0) {
      //if command name is the empty string, we don't include it:
      fprintf(param.commandfs, "%s %s %s %s", param.commandname, this->name, 
		      param.train, fbase);
    } else {
      fprintf(param.commandfs, "%s -O \"%s\" %s %s %s", param.commandname, command, 
		this->name, param.train, fbase);
    }

    //still need to print the class partions:
    for (cls_t i=0; clist[0]+i!=clist[1]; i++) {
      fprintf(param.commandfs, " %d", clist[0][i]);
    }
    fprintf(param.commandfs, " %c", PARTITION_SYMBOL);
    for (cls_t i=0; clist[1]+i!=clist[2]; i++) {
      fprintf(param.commandfs, " %d", clist[1][i]);
    }
    fprintf(param.commandfs, "\n");

    return 2;
  }

  //if we are using svm-predict, can't even figure it out by passing different
  //sized test data since missing dimensions are just treated as zero
  //(svm-predict won't complain if dimensions are missing)
  template <class real, class cls_t>
  dim_ta general2class<real, cls_t>::n_feat() {
    return -1;
  }

  template <class real>
  real logistic_function(real x) {
    return 1-2/(1+exp(x));
  }

  //arctangent normalized to go from [-1,1] with df/dx|_x=0 = 1
  template <class real>
  real atan_norm(real x) {
    return M_PI*atan(2*x/M_PI)/2;
  }

  template <class real, class cls_t>
  agf2class<real, cls_t>::agf2class() {
    brd=NULL;
    grd=NULL;
    n=0;
    gd=NULL;
  }

  template <class real, class cls_t>
  agf2class<real, cls_t>::agf2class(const char *fbase, int sigtype) {
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
	fprintf(stderr, "agf2class: code for function to transform decision values (%d) not recognized\n", sigtype);
	fprintf(stderr, "             [0=tanh; 1=erf; 2=2/(1-exp(..))]--using tanh()\n");
        sigfun=&tanh;
    }
    err=init(fbase, sigfun);
    if (err!=0) throw err;
  }

  template <class real, class cls_t>
  agf2class<real, cls_t>::agf2class(const char *fbase, 
		SIGFUN_TYPE (* sigfun)(SIGFUN_TYPE)) {
    int err=init(fbase, sigfun);
    if (err!=0) throw err;
  }

  template <class real, class cls_t>
  int agf2class<real, cls_t>::init(const char *fbase, 
	 	SIGFUN_TYPE (* sigfun)(SIGFUN_TYPE)) {
    int err;

    this->mat=NULL;
    this->id=-1;

    sigmoid_func=sigfun;

    this->name=new char[strlen(fbase)+1];
    strcpy(this->name, fbase);

    err=agf_read_borders(fbase, brd, grd, n, this->D);
    if (err!=0) return err;

    fprintf(stderr, "agf2class: %d border samples found in model, %s\n", n, fbase);
    this->D1=this->D;
    this->ncls=2;

    //(don't) calculate and store the lengths of all the gradient vectors for possible
    //use later on: (we'll get them when we need them...)
    gd=NULL;

    return err;
  }

  template <class real, class cls_t>
  agf2class<real, cls_t>::agf2class(svm2class<real, cls_t> *svm,
		  real **x, cls_t *cls, dim_ta nvar, nel_ta ntrain,
		  nel_ta ns, real tol, int tflag) {
    bordparam<real> param;
    real *xsort[ntrain];		//sort training data
    cls_t csel[ntrain];
    nel_ta *clind;			//indices for sorted classes
    nel_ta ntrain2;
    int (*sfunc) (void *, real_a *, real_a *);		//sampling function
    real **xtran;

    this->ncls=svm->n_class();
    this->D1=svm->n_feat();
    if (tflag) {
      //is this necessary? shouldn't the training data be transformed already?
      if (this->copy_ltran(svm)) {
        assert(this->D==nvar);
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
      }
    } else {
      assert(this->D1==nvar);
      this->D=this->D1;
      xtran=x;
    }

    //remove side effects:
    for (nel_ta k=0; k<ntrain; k++) {
      csel[k]=cls[k];
      xsort[k]=x[k];
    }
    //sort the classes:
    clind=sort_classes(xsort, ntrain, csel, this->ncls);
    ntrain2=clind[2]-clind[0];

    //initialize parameters:
    if (ns/(ntrain2-clind[1])/clind[1] > 0.25) {
      //for small datasets:
      bordparam_init(&param, xsort+clind[0], nvar, ntrain2, clind[1], 1);
      sfunc=&oppositesample_small<real_a>;
    } else {
      //for large datasets:
      bordparam_init(&param, xsort+clind[0], nvar, ntrain2, clind[1]);
      sfunc=&oppositesample<real_a>;
    }
    param.rparam=svm;

    //allocate the arrays for holding the results:
    brd=allocate_matrix<real, nel_ta>(ns, nvar);
    grd=allocate_matrix<real, nel_ta>(ns, nvar);

    //find class borders:
    n=sample_class_borders(&svmrfunc<real, cls_t>, sfunc, &param, 
		ns, this->D1, tol, agf_global_borders_maxiter, 
		brd, grd);

    delete [] clind;
    if (tflag && this->mat!=NULL) {
      delete [] xtran[0];
      delete [] xtran;
    }

    gd=NULL;
  }

  template <class real, class cls_t>
  agf2class<real, cls_t>::~agf2class() {
    delete_matrix(brd);
    delete_matrix(grd);
    //delete [] this->name;
    if (gd!=NULL) delete [] gd;
  }

  template <class real, class cls_t>
  void agf2class<real, cls_t>::calc_grad_len() {
    //calculate and store the lengths of all the gradient vectors for possible
    //use later on:
    gd=new real[n];
    for (nel_ta i=0; i<n; i++) {
      gd[i]=0;
      for (dim_ta j=0; j<this->D; j++) gd[i]+=grd[i][j]*grd[i][j];
      gd[i]=sqrt(gd[i]);
    }
  }

  //if flag, then border vectors are stored un-normalized
  template <class real, class cls_t>
  int agf2class<real, cls_t>::ltran_model(real **mat1, real *b1, dim_ta d1, dim_ta d2) {
    int err2=0;

    real **brd2;
    real **grd2;

    gsl_matrix *u;
    gsl_vector *s;
    gsl_matrix *vt;
    gsl_vector *work;

    if (d1!=this->D1) {
      fprintf(stderr, "agf2class: first dimension (%d) of trans. mat. does not agree with that of borders data (%d)\n", d1, this->D1);
      return DIMENSION_MISMATCH;
    }
    fprintf(stderr, "agf2class: Normalising the border samples...\n");

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
  real agf2class<real, cls_t>::R(real *x, real *praw) {
    real r;
    nel_ta k;		//intermediate values in the calculation
    real d;		//may be useful for continuum generalization

    r=border_classify0(brd, grd, this->D1, n, x, k, d);
    if (this->id>=0 && praw!=NULL) {
      praw[this->id]=r;
      //printf("r=%g\n" , praw[this->id]);
    }
    if (sigmoid_func != NULL) r=(*sigmoid_func) (r);

    return r;
  }

  template <class real, class cls_t>
  void agf2class<real, cls_t>::print(FILE *fs, char *fbase, int depth) {
    for (int i=0; i<depth; i++) fprintf(fs, "  ");
    if (fbase==NULL) {
      fprintf(fs, "%s", this->name);
    } else {
      fprintf(fs, "%s", fbase);
    }
  }

  template <class real, class cls_t>
  int agf2class<real, cls_t>::load(FILE *fs) {
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
    return err;
  }

  template <class real, class cls_t>
  int agf2class<real, cls_t>::save(FILE *fs) {
    print_matrix(fs, brd, n, this->D1);
    print_matrix(fs, grd, n, this->D1);
    return 0;
  }

  template float logistic_function<float>(float);
  template double logistic_function<double>(double);

  template float atan_norm<float>(float);
  template double atan_norm<double>(double);

  template class binaryclassifier<real_a, cls_ta>;
  template class general2class<real_a, cls_ta>;
  template class agf2class<real_a, cls_ta>;

}

