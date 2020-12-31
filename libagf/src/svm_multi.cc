#include <string.h>
#include <math.h>
#include <assert.h>

#include <gsl/gsl_linalg.h>

#include "constrained.h"
#include "read_ascii_all.h"
#include "error_codes.h"
#include "gsl_util.h"
#include "full_util.h"

#include "agf_lib.h"

using namespace libpetey;

namespace libagf {

  template <class real>
  int solve_cond_prob_1v1(real **r, int ncls, real *p) {
    int err=0;
    gsl_matrix *Q;
    gsl_vector *x;
    gsl_vector *b;

    //Wu and Lin 2004, JMLR 5:975 give linear system to solve for multi-class
    //probabilities in 1 vs. 1 case
    Q=gsl_matrix_alloc(ncls+1, ncls+1);
    x=gsl_vector_alloc(ncls+1);
    b=gsl_vector_alloc(ncls+1);

    //fill the matrix and the solution vector:
    //off-diagonals:
    for (int i=0; i<ncls; i++) {
      for (int j=i+1; j<ncls; j++) {
        real val=-r[i][j]*(1-r[i][j]);
        gsl_matrix_set(Q, i, j, val);
	gsl_matrix_set(Q, j, i, val);
      }
    }

    //diagonals plus solution vector:
    for (int i=0; i<ncls; i++) {
      gsl_vector_set(b, i, 0);
      real val=0;
      for (int j=0; j<i; j++) val+=(1-r[j][i])*(1-r[j][i]);
      for (int j=i+1; j<ncls; j++) val+=r[i][j]*r[i][j];
      gsl_matrix_set(Q, i, i, val);
      //normality constraint:
      gsl_matrix_set(Q, i, ncls, 1);
      gsl_matrix_set(Q, ncls, i, 1);
    }
    //rest of normality constraint:
    gsl_vector_set(b, ncls, 1);
    gsl_matrix_set(Q, ncls, ncls, 0);

    //print_gsl_matrix(stdout, Q);
    //gsl_vector_fprintf(stdout, b, "%lg ");
    //printf("\n");

    //use SVD solver:
    err=gsl_lsq_solver(Q, b, x);

    //gsl_vector_fprintf(stdout, x, "%lg ");
    //printf("\n");

    //re-assign GSL result back to standard floating point array:
    for (int i=0; i<ncls; i++) p[i]=gsl_vector_get(x, i);

    //clean up:
    gsl_matrix_free(Q);
    gsl_vector_free(x);
    gsl_vector_free(b);

    return err;
  }

  template <class real, class cls_t>
  onevone<real, cls_t>::onevone() {
    this->ncls=0;
    this->D=0;
    label=NULL;
    voteflag=0;
  }

  template <class real, class cls_t>
  onevone<real, cls_t>::~onevone() {
    delete [] label;
  }

  template <class real, class cls_t>
  cls_t onevone<real, cls_t>::classify(real *x, real *p, real *praw) {
    real **praw0;
    int k=0;

    praw0=classify_raw(x);
    if (voteflag) {
      for (int i=0; i<this->ncls; i++) p[i]=0;
      for (int i=0; i<this->ncls; i++) {
        for (int j=i+1; j<this->ncls; j++) {
          //sign is reversed relative to LIBSVM implementation:
          if (praw0[i][j]<0) p[i]++; else p[j]++;
	}
      }
    } else {
      solve_cond_prob_1v1(praw0, this->ncls, p);
    }

    k=0;
    if (praw!=NULL) {
      for (int i=0; i<this->ncls; i++) {
        for (int j=i+1; j<this->ncls; j++) {
          praw[k]=praw0[i][j];
	  k++;
	}
      }
    }
/*
    for (int i=0; i<this->ncls; i++) {
      for (int j=0; j<i+1; j++) praw0[i][j]=0;
    }

    print_matrix(stdout, praw0, this->ncls, this->ncls);
*/
    delete [] praw0[0];
    delete [] praw0;

    return label[choose_class(p, this->ncls)];
  }

  template <class real, class cls_t>
  cls_t onevone<real, cls_t>::class_list(cls_t *cls) {
    for (cls_t i=0; i<this->ncls; i++) cls[i]=label[i];
    return this->ncls;
  }

  template <class real, class cls_t>
  svm_multi<real, cls_t>::svm_multi() {
    this->ncls=0;
    this->D=0;
    nsv_total=0;
    nsv=NULL;
    sv=NULL;
    coef=NULL;
    rho=NULL;
    probA=NULL;
    probB=NULL;
    this->label=NULL;
    this->voteflag=0;
  }

  //initialize from LIBSVM model file:
  template <class real, class cls_t>
  svm_multi<real, cls_t>::svm_multi(char *file, int vflag) {
    FILE *fs=fopen(file, "r");
    if (fs==NULL) {
      fprintf(stderr, "svm2class: failed to open model file, %s\n", file);
      throw UNABLE_TO_OPEN_FILE_FOR_READING;
    }
    this->voteflag=vflag;
    probA=NULL;
    probB=NULL;
    load(fs);
    fclose(fs);
  }

  //initialize from LIBSVM model file:
  template <class real, class cls_t>
  int svm_multi<real, cls_t>::load(FILE *fs) {
    char *line=NULL;
    char **substr=NULL;
    int nsub;
    int nparam;		//ncls*(ncls-1)/2: number of binary classifiers
    int nsv1;		//should agree with nsv...
    int pfound=0;	//number of parameters for probability estimation read in (should be 2)
    char format[12];	//format code
    char fcode[4];	//floating point format code
    int lineno=0;

    get_format_code<real>(fcode);
    sprintf(format, "%%%s%%n", fcode);

    //this->mat=NULL;
    //this->b=NULL;
    //this->id=-1;

    rho=0;
    param=new real[3];

    do {
      if (line!=NULL) delete [] line;
      if (substr!=NULL) delete [] substr;
      line=fget_line(fs, 1);
      lineno++;
      substr=split_string_destructive(line, nsub);
      //printf("svm2class: nsub=%d\n", nsub);
      //for (int i=0; i<nsub; i++) printf("%s ", substr[i]);
      //printf("\n");
      if (nsub == 0) continue;
      if (strcmp(substr[0], "SV")!=0 && nsub<2) {
        fprintf(stderr, "svm_multi: error in initialization file; unrecognized keywordi (%s)/not enough parameters\n", substr[0]);
        throw FILE_READ_ERROR;
      }
      if (strcmp(substr[0], "svm_type")==0) {
        if (strcmp(substr[1], "c_svc")!=0 && strcmp(substr[1], "nu-svc")!=0) {
          fprintf(stderr, "svm_multi: not a classifier SVM: %s\n", substr[1]);
	  throw PARAMETER_OUT_OF_RANGE;
        }
      } else if (strcmp(substr[0], "nr_class")==0) {
        this->ncls=atoi(substr[1]);
	if (this->ncls < 2) {
          fprintf(stderr, "svm_multi: one class classifiers not accepted\n");
	  throw PARAMETER_OUT_OF_RANGE;
        }
	nparam=this->ncls*(this->ncls-1)/2;
      } else if (strcmp(substr[0], "kernel_type")==0) {
        if (strcmp(substr[1], "linear")==0) {
	  kernel=&linear_basis<real>;
	  kernel_deriv=&linear_basis_deriv<real>;
	  param=NULL;
        } else if (strcmp(substr[1], "polynomial")==0) {
          kernel=&polynomial_basis<real>;
          kernel_deriv=&polynomial_basis_deriv<real>;
	} else if (strcmp(substr[1], "rbf")==0) {
          kernel=&radial_basis<real>;
          kernel_deriv=&radial_basis_deriv<real>;
	} else if (strcmp(substr[1], "sigmoid")==0) {
          kernel=&sigmoid_basis<real>;
          kernel_deriv=&sigmoid_basis_deriv<real>;
	} else {
          fprintf(stderr, "svm2class: basis function, %s, not recognized\n", substr[1]);
	  throw PARAMETER_OUT_OF_RANGE;
        }
      } else if (strcmp(substr[0], "gamma")==0) {
        param[0]=atof(substr[1]);
	//printf("gamma=%g\n", param[0]);
      } else if (strcmp(substr[0], "coef0")==0) {
        param[1]=atof(substr[1]);
      } else if (strcmp(substr[0], "degree")==0) {
        param[2]=atof(substr[1]);
      } else if (strcmp(substr[0], "total_sv")==0) {
        nsv_total=atoi(substr[1]);
      } else if (strcmp(substr[0], "nr_sv")==0) {
        if (nsub<this->ncls+1) {
          fprintf(stderr, "svm_multi: error in initialization file: not enough parameters (nr_sv)\n");
	  throw FILE_READ_ERROR;
	}
	nsv=new nel_ta[this->ncls];
	for (int i=0; i<this->ncls; i++) nsv[i]=atoi(substr[i+1]);
	//for (int i=0; i<this->ncls; i++) printf("%d ", nsv[i]);
	//printf("\n");
      } else if (strcmp(substr[0], "rho")==0) {
        if (nsub<nparam+1) {
          fprintf(stderr, "svm_multi: error in initialization file: not enough parameters (rho)\n");
	  throw FILE_READ_ERROR;
	}
        rho=new real[nparam];
	for (int i=0; i<nparam; i++) rho[i]=atof(substr[i+1]);
      } else if (strcmp(substr[0], "probA")==0) {
        if (nsub<nparam+1) {
          fprintf(stderr, "svm_multi: error in initialization file: not enough parameters (probA)\n");
	  throw FILE_READ_ERROR;
	}
        probA=new real[nparam];
	for (int i=0; i<nparam; i++) {
          probA[i]=atof(substr[i+1]);
          //printf("probA[%d]=%g\n", i, probA[i]);
	}
	pfound++;
      } else if (strcmp(substr[0], "probB")==0) {
        if (nsub<nparam+1) {
          fprintf(stderr, "svm_multi: error in initialization file: not enough parameters (probB)\n");
	  throw FILE_READ_ERROR;
	}
        probB=new real[nparam];
	for (int i=0; i<nparam; i++) {
          probB[i]=atof(substr[i+1]);
          //printf("probB[%d]=%g\n", i, probB[i]);
	}
	pfound++;
      } else if (strcmp(substr[0], "label")==0) {
        if (nsub<this->ncls+1) {
          fprintf(stderr, "svm_multi: error in initialization file: not enough parameters (label)\n");
	  throw FILE_READ_ERROR;
	}
        this->label=new cls_t[this->ncls];
	for (int i=0; i<this->ncls; i++) this->label[i]=atoi(substr[i+1]);
      }
    } while (strcmp(substr[0], "SV")!=0);
    delete [] line;
    delete [] substr;

    dim_ta *ind[nsv_total];		//dimension indices
    real *raw[nsv_total];		//raw features data
    int nf[nsv_total];
    coef=new real*[this->ncls-1];
    coef[0]=new real[nsv_total*(this->ncls-1)];
    for (int i=1; i<this->ncls-1; i++) coef[i]=coef[0]+i*nsv_total;
    this->D=0;
    fprintf(stderr, "svm_multi: reading in %d support vectors\n", nsv_total);
    for (int i=0; i<nsv_total; i++) {
      int nread;		//number of item scanned
      int pos=0;		//position in line
      int rel;			//relative position in line

      line=fget_line(fs);
      for (int j=0; j<this->ncls-1; j++) {
        nread=sscanf(line+pos, format, coef[j]+i, &rel);
	//printf("%g ", coef[j][i]);
	if (nread!=1) {
          fprintf(stderr, "svm_multi: error reading coefficients, line %d\n", lineno+i);
	  throw FILE_READ_ERROR;
	}
	pos+=rel;
      }
      nf[i]=scan_svm_features(line+pos, ind[i], raw[i]);
      if (nf[i]<=0) {
        fprintf(stderr, "svm_multi: error reading support vectors, line %d\n", lineno+i);
	throw FILE_READ_ERROR;
      }
      for (int j=0; j<nf[i]; j++) {
        if (ind[i][j]>this->D1) this->D1=ind[i][j];
	//printf("%d:%g ", ind[i][j], raw[i][j]);
      }
      //printf("\n");
      delete [] line;
    }
    //transfer to more usual array and fill in missing values:
    real missing=0;
    sv=new real*[nsv_total];
    sv[0]=new real[nsv_total*this->D1];
    for (int i=0; i<nsv_total; i++) {
      sv[i]=sv[0]+i*this->D1;
      for (int j=0; j<this->D1; j++) sv[i][j]=missing;
      for (int j=0; j<nf[i]; j++) sv[i][ind[i][j]-1]=raw[i][j];
      delete [] ind[i];
      delete [] raw[i];
    }
    fprintf(stderr, "svm_multi: read in %d support vectors\n", nsv_total);
    //printf("param[0]=%g\n", param[0]);

    //some book-keeping:
    if (probA==NULL || probB==NULL) this->voteflag=1;
    this->D=this->D1;

    start=new nel_ta[this->ncls];
    start[0]=0;
    for (int i=1; i<this->ncls; i++) start[i]=start[i-1]+nsv[i-1];

    //cls_t cls[nsv_total];
    //print_lvq_svm(stdout, sv, cls, nsv_total, this->D1, 1);
    return 0;
  }

  //copy constructor:
  template <class real, class cls_t>
  svm_multi<real, cls_t>::svm_multi(svm_multi<real, cls_t> *other) {
    int nmod;			//number of binary classifiers
    this->ncls=other->ncls;
    nmod=this->ncls*(this->ncls-1)/2;
    this->D1=other->D1;
    this->D=other->D;
    //number of support vectors:
    nsv_total=other->nsv_total;
    nsv=new nel_ta[this->ncls];
    for (int i=0; i<this->ncls; i++) nsv[i]=other->nsv[i];
    //support vectors themselves:
    sv=new real *[nsv_total];
    sv[0]=new real[nsv_total*this->D1];
    for (int i=0; i<nsv_total; i++) {
      sv[i]=sv[0]+i*this->D1;
      for (int j=0; j<this->D1; j++) sv[i][j]=other->sv[i][j];
    }
    //coefficients:
    coef=new real *[this->ncls-1];
    coef[0]=new real[nsv_total*(this->ncls-1)];
    for (int i=0; i<this->ncls-1; i++) {
      coef[i]=coef[0]+i*nsv_total;
      for (int j=0; j<nsv_total; j++) coef[i][j]=other->coef[i][j];
    }
    //constant term:
    rho=new real[nmod];
    for (int i=0; i<nmod; i++) rho[i]=other->rho[i];
    //kernel functions:
    param=new real[3];
    for (int i=0; i<nmod; i++) param[i]=other->param[i];
    kernel=other->kernel;
    kernel_deriv=other->kernel_deriv;
    //coefficients for determining probabilities:
    if (other->probA!=NULL && other->probB!=NULL) {
      probA=new real[nmod];
      probB=new real[nmod];
      for (int i=0; i<nmod; i++) {
        probA[i]=other->probA[i];
        probB[i]=other->probB[i];
      }
    }
    //class labels (why aren't they in order??):
    this->label=new cls_t[this->ncls];
    for (int i=0; i<this->ncls; i++) this->label[i]=other->label[i];
    //book-keeping:
    start=new nel_ta[this->ncls];
    start[0]=0;
    for (int i=1; i<this->ncls; i++) start[i]=start[i-1]+nsv[i-1];
    this->voteflag=other->voteflag;
  }

  template <class real, class cls_t>
  svm_multi<real, cls_t>::~svm_multi() {
    delete [] sv[0];
    delete [] sv;
    delete [] coef[0];
    delete [] coef;
    delete [] rho;
    delete [] nsv;
    if (probA!=NULL) delete [] probA;
    if (probB!=NULL) delete [] probB;
    delete [] start;
    delete [] param;
  }

  //raw decision values, one for each pair of classes:
  template <class real, class cls_t>
  real ** svm_multi<real, cls_t>::classify_raw(real *x) {
    real **result;
    real kv[nsv_total];
    int si, sj;
    int p=0;

    //calculate kernels for each support vector:
    for (int i=0; i<nsv_total; i++) {
      kv[i]=(*kernel)(x, sv[i], this->D1, param);
    }

    result=new real *[this->ncls];
    //wastes space, but it's simpler this way:
    result[0]=new real[this->ncls*this->ncls];

    //multiply kernels by the coefficients:
    for (int i=0; i<this->ncls; i++) {
      result[i]=result[0]+i*this->ncls;
      si=start[i];
      for (int j=i+1; j<this->ncls; j++) {
	sj=start[j];
        result[i][j]=0;
	//coefficients are organized rather awkardly:
	for (int k=0; k<nsv[i]; k++) result[i][j]+=coef[j-1][si+k]*kv[si+k];
	for (int k=0; k<nsv[j]; k++) result[i][j]+=coef[i][sj+k]*kv[sj+k];
	result[i][j] -= rho[p];
	//sign is reversed relative to LIBSVM implementation:
        if (this->voteflag) {
          result[i][j]=-result[i][j];
	} else {
          //approximate probabilities:
          result[i][j]=1-1./(1+exp(probA[p]*result[i][j]+probB[p]));
	}
	p++;
      }
    }

    return result;
  }

  //raw decision value for a given pair of classes:
  template <class real, class cls_t>
  real svm_multi<real, cls_t>::R(real *x, cls_t i, cls_t j, real *praw) {
    real result;
    real kv;
    int si, sj;
    cls_t swp;
    int p=0;
    int sgn=1;

    if (i==j) throw PARAMETER_OUT_OF_RANGE;

    if (i>j) {
      swp=i;
      i=j;
      j=swp;
      sgn=-1;
    }

    si=start[i];
    sj=start[j];

    result=0;
    for (int k=0; k<nsv[i]; k++) {
      kv=(*kernel)(x, sv[si+k], this->D1, param);
      result+=coef[j-1][si+k]*kv;
      //printf("kernel(%d)=%g\n", si+k, kv);
    }
    for (int k=0; k<nsv[j]; k++) {
      kv=(*kernel)(x, sv[sj+k], this->D1, param);
      result+=coef[i][sj+k]*kv;
      //printf("kernel(%d)=%g\n", sj+k, kv);
    }
    //p=this->ncls*(this->ncls-1)/2-(this->ncls-i)*(this->ncls-i-1)/2+j;
    p=i*(2*this->ncls-i-1)/2+j-i-1;
    //printf("%d %d %d\n", i, j, p);
    result -= rho[p];
    if (praw!=NULL) praw[0]=sgn*result;
    if (this->voteflag==0) {
      result=2./(1+exp(probA[p]*result+probB[p]))-1;
    }

    return sgn*result;
  }

  //raw decision value for a given pair of classes:
  template <class real, class cls_t>
  real svm_multi<real, cls_t>::R_deriv(real *x, cls_t i, cls_t j, real *drdx) {
    real result;
    real kv;
    int si, sj;
    cls_t swp;
    int sgn=1;
    int p=0;
    real t1, t2;
    real deriv[this->D1];
    real drdx1[this->D1];

    if (i==j) throw PARAMETER_OUT_OF_RANGE;

    if (i>j) {
      swp=i;
      i=j;
      j=swp;
      sgn=-1;
    }

    si=start[i];
    sj=start[j];

    result=0;
    for (dim_ta k=0; k<this->D1; k++) drdx1[k]=0;
    for (int k=0; k<nsv[i]; k++) {
      kv=(*kernel_deriv)(x, sv[si+k], this->D1, param, deriv);
      result+=coef[j-1][si+k]*kv;
      for (dim_ta m=0; m<this->D1; m++) drdx1[m]+=coef[j-1][si+k]*deriv[m];
    }
    for (int k=0; k<nsv[j]; k++) {
      kv=(*kernel_deriv)(x, sv[sj+k], this->D1, param, deriv);
      result+=coef[i][sj+k]*kv;
      for (dim_ta m=0; m<this->D1; m++) drdx1[m]+=coef[i][sj+k]*deriv[m];
    }
    //p=this->ncls*(this->ncls-1)/2-(this->ncls-i)*(this->ncls-i-1)/2+j;
    p=i*(2*this->ncls-i-1)/2+j-i-1;
    result -= rho[p];
    if (this->voteflag==0) {
      t1=exp(result*probA[p]+probB[p]);
      t2=probA[p]*t1/(1+t1)/(1+t1);		//derivative of sigmoid function
      for (dim_ta k=0; k<this->D1; k++) drdx1[k]*=-2*t2;
      result=2/(1+t1)-1;
    }

    //printf("drdx=");
    for (dim_ta k=0; k<this->D1; k++) {
      drdx[k]=sgn*drdx1[k];
      //printf("%g ", drdx[k]);
    }
    //printf("\nr=%g\n", result);

    return sgn*result;
  }

  template <class real, class cls_t>
  int svm_multi<real, cls_t>::ltran_model(real **mat1, real *b1, dim_ta d1, dim_ta d2) {
    real **sv2;
    if (d1!=this->D1) {
      fprintf(stderr, "svm_multi: first dimension (%d) of trans. mat. does not agree with that of borders data (%d)\n", d1, this->D1);
      return DIMENSION_MISMATCH;
    }
    //apply constant factor:
    for (nel_ta i=0; i<nsv_total; i++) {
      for (dim_ta j=0; j<this->D1; j++) {
        sv[i][j]=sv[i][j]-b1[j];
      }
    }

    sv2=matrix_mult(sv, mat1, nsv_total, d1, d2);

    if (this->mat == NULL) {
      //the classifier now has a different number of features:
      this->D1=d2;
      this->D=d2;
    } else {
      //assert(mat1==this->mat && b1==this->b);
      //from the outside, the classifier looks like it has the same number of
      //features as before normalization:
      assert(this->D1==d2);
    }

    delete [] sv[0];
    delete [] sv;
    sv=sv2;

    return 0;
  }

  template <class real, class cls_t>
  int svm_multi<real, cls_t>::save(FILE *fs) {
    char fcode[4];
    char format[10];
    int nparam=this->ncls*(this->ncls-1)/2;	//number of binary classifiers

    get_format_code<real>(fcode);
    sprintf(format, " %%%s", fcode);

    fprintf(fs, "svm_type c_svc\n");
    //kernel type and parameters:
    fprintf(fs, "kernel_type ");
    if (kernel==&linear_basis<real>) {
      fprintf(fs, "linear\n");
    } else if (kernel==&polynomial_basis<real>) {
      fprintf(fs, "polynomial\n");
      fprintf(fs, "gamma");
      fprintf(fs, format, param[0]);
      fprintf(fs, "\n");
      fprintf(fs, "coef0");
      fprintf(fs, format, param[1]);
      fprintf(fs, "\n");
      fprintf(fs, "degree");
      fprintf(fs, format, param[2]);
      fprintf(fs, "\n");
    } else if (kernel==&radial_basis<real>) {
      fprintf(fs, "rbf\n");
      fprintf(fs, "gamma");
      fprintf(fs, format, param[0]);
      fprintf(fs, "\n");
    } else if (kernel=&sigmoid_basis<real>) {
      fprintf(fs, "sigmoid\n");
      fprintf(fs, "gamma");
      fprintf(fs, format, param[0]);
      fprintf(fs, "\n");
      fprintf(fs, "coef0");
      fprintf(fs, format, param[1]);
      fprintf(fs, "\n");
    }
    //number of classes:
    fprintf(fs, "nr_class %d\n", this->ncls);
    //total number of support vectors:
    fprintf(fs, "total_sv %d\n", nsv_total);
    //constant parameter:
    fprintf(fs, "rho");
    for (int i=0; i<nparam; i++) fprintf(fs, format, rho[i]);
    fprintf(fs, "\n");
    //class labels:
    fprintf(fs, "label");
    for (int i=0; i<this->ncls; i++) fprintf(fs, " %d", this->label[i]);
    fprintf(fs, "\n");
    //for calculating probabilities:
    if (probA!=NULL) {
      fprintf(fs, "probA");
      for (int i=0; i<nparam; i++) fprintf(fs, format, probA[i]);
      fprintf(fs, "\n");
      fprintf(fs, "probB");
      for (int i=0; i<nparam; i++) fprintf(fs, format, probB[i]);
      fprintf(fs, "\n");
    }
    //support vectors for each class:
    fprintf(fs, "nr_sv");
    for (int i=0; i<this->ncls; i++) fprintf(fs, " %d", nsv[i]);
    fprintf(fs, "\n");
    //support vectors and corresponding coefficients:
    fprintf(fs, "SV\n");
    for (int i=0; i<nsv_total; i++) {
      for (int j=0; j<this->ncls-1; j++) fprintf(fs, format, coef[j][i]);
      for (int j=0; j<this->D1; j++) {
        fprintf(fs, " %d:", j+1);
        fprintf(fs, format, sv[i][j]);
      }
      fprintf(fs, "\n");
    }
    return 0;
    
  }

  template <class real, class cls_t>
  void svm_multi<real, cls_t>::subsample(real frac, int cclassflag) {
    nel_ta *nsv_new;
    nel_ta t1, t2;

    long *rind;

    nsv_new=new nel_ta[this->ncls];
    for (int i=0; i<this->ncls; i++) {
      nsv_new[i]=frac*nsv[i];
    }

    t1=0;
    t2=0;
    for (cls_t i=0; i<this->ncls; i++) {
      rind=randomize(nsv[i]);
      for (nel_ta j=0; j<nsv_new[i]; j++) {
        sv[j+t1]=sv[rind[j]+t2];
	for (cls_t k=0; k<this->ncls-1; k++) {
          coef[k][j+t1]=coef[k][rind[j]+t2];
	}
      }
      t1+=nsv_new[i];
      t2+=nsv[i];
      delete [] rind;
    }

    nsv_total=t1;

    delete [] nsv;
    nsv=nsv_new;

  }

  template <class real, class cls_t>
  borders1v1<real, cls_t>::borders1v1(char *file, int vflag) {
    FILE *fs=fopen(file, "r");
    if (fs==NULL) throw UNABLE_TO_OPEN_FILE_FOR_READING;
    this->voteflag=vflag;
    load(fs);
    fclose(fs);
  }

  template <class real, class cls_t>
  int borders1v1<real, cls_t>::load(FILE *fs) {
    //how do we want to store the borders?
    int err=0;
    //read in number of classes:
    char *line=fget_line(fs, 1);
    if (strcmp(line, "1v1")!=0) {
      fprintf(stderr, "borders1v1: initialization file wrong type (%s)\n", line);
      throw PARAMETER_OUT_OF_RANGE;
    }
    delete [] line;
    line=fget_line(fs);
    sscanf(line, "%d", &this->ncls);
    delete [] line;
    int nmod=(this->ncls-1)*this->ncls/2;
    this->label=new cls_t[this->ncls];
    //read in class labels:
    line=fget_line(fs);
    char **sub;
    int nsub;
    sub=split_string_destructive(line, nsub);
    if (nsub<this->ncls) {
      fprintf(stderr, "borders1v1: Not enough labels found in initialization file\n");
      fprintf(stderr, "  %d vs. %d\n", this->ncls, nsub);
      throw DIMENSION_MISMATCH;
    }
    for (cls_t i=0; i<this->ncls; i++) this->label[i]=atoi(sub[i]);
    delete [] line;
    delete [] sub;
    line=fget_line(fs);
    sub=split_string_destructive(line, nsub);
    if (nsub<nmod) {
      fprintf(stderr, "borders1v1: Not enough \"polarities\" found in initialization file\n");
      fprintf(stderr, "  %d vs. %d\n", this->ncls, nsub);
      throw DIMENSION_MISMATCH;
    }
    pol=new int[nmod];
    for (int i=0; i<nmod; i++) pol[i]=atoi(sub[i]);
    delete [] line;
    delete [] sub;
    //read in individual binary classifiers:
    classifier=new borders_classifier<real, cls_t>*[nmod];
    for (int i=0; i<nmod; i++) {
      classifier[i]=new borders_calibrated<real, cls_t>();
      err=classifier[i]->load(fs, this->voteflag);
      if (err!=0) throw err;
    }

    //book-keeping--check dimensions and make sure they are all the same:
    dim_ta D1, D;
    this->D1=classifier[0]->n_feat();
    this->D=classifier[0]->n_feat_t();
    for (int i=1; i<nmod; i++) {
      D1=classifier[i]->n_feat();
      D=classifier[i]->n_feat_t();
      if (D1!=this->D1 || D!=this->D) {
        fprintf(stderr, "borders1v1: Dimension of classifier %d does not match that of first one in file\n", i);
        fprintf(stderr, "  %d vs. %d\n", this->D1, D1);
        throw DIMENSION_MISMATCH;
      }
    }
    return 0;
  }

  template <class real, class cls_t>
  int borders1v1<real, cls_t>::save(FILE *fs) {
    int err=0;
    int nmod=this->ncls*(this->ncls-1)/2;
    fprintf(fs, "1v1\n");
    fprintf(fs, "%d\n", this->ncls);
    for (cls_t i=0; i<this->ncls; i++) fprintf(fs, "%d ", this->label[i]);
    fprintf(fs, "\n");
    for (int i=0; i<nmod; i++) {
      fprintf(fs, "%d ", pol[i]);
    }
    fprintf(fs, "\n");
    for (int i=0; i<nmod; i++) {
      err=classifier[i]->save(fs);
      if (err!=0) return err;
    }
    return err;
  }

  template <class real, class cls_t>
  borders1v1<real, cls_t>::borders1v1(svm_multi<real, cls_t> *svm,
		  real **x, cls_t *cls, dim_ta nvar, nel_ta ntrain,
		  nel_ta ns, real tol, int tflag) {
    int nmod;
    int m=0;
    svm2class<real, cls_t> *svmbin;
    cls_t csel[ntrain];			//selected classes
    real **xtran;			//transformed training data
    int err;				//error code

    this->ncls=svm->n_class();
    this->D1=svm->n_feat();
    this->label=new cls_t[this->ncls];
    svm->class_list(this->label);
    if (tflag && this->copy_ltran(svm)) {
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
    } else {
      assert(this->D1==nvar);
      this->D=this->D1;
      xtran=x;
    }

    nmod=this->ncls*(this->ncls-1)/2;
    classifier=new borders_classifier<real, cls_t>*[nmod];

    for (int i=0; i<this->ncls; i++) {
      for (int j=i+1; j<this->ncls; j++) {
	//select out pair of classes:
        for (nel_ta k=0; k<ntrain; k++) {
          if (cls[k]==this->label[i]) {
            csel[k]=0;
          } else if (cls[k]==this->label[j]) {
            csel[k]=1;
          } else {
            csel[k]=-1;
          }
        }
	fprintf(stderr, "Training %d vs. %d\n", this->label[i], this->label[j]);
        svmbin=new svm2class<real, cls_t>(svm, i, j);
	classifier[m]=new borders_calibrated<real, cls_t>();
	classifier[m]->train(svmbin, xtran, csel, this->D1, ntrain, ns, tol);
	m++;
	delete svmbin;
	continue;

        //lets test the result:
	printf("%d vs %d comparison:\n", this->label[i], this->label[j]);
	for (int ti=0; ti<100; ti++) {
          real r1, r2;
	  int ind=ranu()*ntrain;
	  r1=svmbin->R(xtran[ind]);
	  r2=classifier[m]->R(xtran[ind]);
	  printf("%g %g\n", r1, r2);
	}
      }
    }

    pol=new int[nmod];
    for (int i=0; i<nmod; i++) pol[i]=-1;

    if (tflag && this->mat!=NULL) {
      delete [] xtran[0];
      delete [] xtran;
    }
  }

  template <class real, class cls_t>
  borders1v1<real, cls_t>::~borders1v1() {
    for (int i=0; i<this->ncls*(this->ncls-1)/2; i++) {
      delete classifier[i];
    }
    delete [] classifier;
    delete [] pol;
  }

  template <class real, class cls_t>
  real **borders1v1<real, cls_t>::classify_raw(real *x) {
    int k;
    real **result=new real *[this->ncls];
    result[0]=new real[this->ncls*this->ncls];
    k=0;
    for (int i=0; i<this->ncls; i++) {
      result[i]=result[0]+i*this->ncls;
      for (int j=i+1; j<this->ncls; j++) {
        result[i][j]=pol[k]*classifier[k]->R(x);
	if (this->voteflag!=1) result[i][j]=(result[i][j]+1)/2;
	k++;
      }
    }
    return result;
  }

  template <class real, class cls_t>
  int borders1v1<real, cls_t>::ltran_model(real **mat1, real *b1, dim_ta d1, dim_ta d2) {
    int nmod;
    nmod=this->ncls*(this->ncls-1)/2;
    for (int i=0; i<nmod; i++) {
      classifier[i]->ltran_model(mat1, b1, d1, d2);
    }
    if (this->mat == NULL) {
      //the classifier now has a different number of features:
      this->D1=d2;
      this->D=d2;
    } else {
      //assert(mat1==this->mat && b1==this->b);
      //from the outside, the classifier looks like it has the same number of
      //features as before normalization:
      assert(this->D1==d2);
    }

    return 0;
  }

  template class svm_multi<float, cls_ta>;
  template class svm_multi<double, cls_ta>;
  template class borders1v1<float, cls_ta>;
  template class borders1v1<double, cls_ta>;

} //end namespace libagf

