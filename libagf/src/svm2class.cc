#include <math.h>
#include <string.h>

#include "read_ascii_all.h"
#include "error_codes.h"
#include "full_util.h"
#include "agf_lib.h"

using namespace libpetey;

namespace libagf {
  //basis functions ("kernels"):
  template <class real>
  real linear_basis(real *x, real *y, dim_ta n, void *param) {
    real dot=0;
    for (dim_ta i=0; i<n; i++) dot+=x[i]*y[i];
    return dot;
  }

  template <class real>
  real polynomial_basis(real *x, real *y, dim_ta n, void *param) {
    real *p2=(real *) param;
    real gamma=p2[0];
    real coef0=p2[1];
    real degree=p2[2];
    real dot=linear_basis(x, y, n, NULL);
    return pow(gamma*dot+coef0, degree);
  }

  template <class real>
  real radial_basis(real *x, real *y, dim_ta n, void *param) {
    real gamma=((real *) param)[0];
    real d2=metric2(x, y, n);
    return exp(-gamma*d2);
  }

  template <class real>
  real sigmoid_basis(real *x, real *y, dim_ta n, void *param) {
    real *p2=(real *) param;
    real gamma=p2[0];
    real coef0=p2[1];
    real dot=linear_basis(x, y, n, NULL);
    return tanh(gamma*dot+coef0);
  }

  //basis functions ("kernels"):
  template <class real>
  real linear_basis_deriv(real *x, real *y, dim_ta n, void *param, real *deriv) {
    real dot=0;
    for (dim_ta i=0; i<n; i++) {
      dot+=x[i]*y[i];
      deriv[i]=y[i];
    }
    return dot;
  }

  template <class real>
  real polynomial_basis_deriv(real *x, real *y, dim_ta n, void *param, real *deriv) {
    real *p2=(real *) param;
    real gamma=p2[0];
    real coef0=p2[1];
    real degree=p2[2];
    real dot=linear_basis_deriv(x, y, n, NULL, deriv);
    for (dim_ta i=0; i<n; i++) {
      deriv[i]*=pow(gamma*dot+coef0, degree-1)*degree*gamma;
    }
    return pow(gamma*dot+coef0, degree);
  }

  template <class real>
  void radial_basis_deriv_num(real *x, real *y, dim_ta n, real dx, void *param, real *deriv) {
    real x1[n];
    real x2[n];
    real r1, r2;
    for (dim_ta i=0; i<n; i++) {
      x1[i]=x[i];
      x2[i]=x[i];
    }
    for (dim_ta i=0; i<n; i++) {
      x1[i]-=dx;
      x2[i]+=dx;
      r1=radial_basis(x1, y, n, param);
      r2=radial_basis(x2, y, n, param);
      deriv[i]=(r2-r1)/dx/2;
      x1[i]=x[i];
      x2[i]=x[i];
    }
  }

  template <class real>
  real radial_basis_deriv(real *x, real *y, dim_ta n, void *param, real *deriv) {
    real deriv1[n];
    real gamma=((real *) param)[0];
    real t1;
    real d2=0;
    for (dim_ta i=0; i<n; i++) {
      deriv[i]=x[i]-y[i];
      d2+=deriv[i]*deriv[i];
    }
    t1=exp(-gamma*d2);
    for (dim_ta i=0; i<n; i++) deriv[i]*=-2*t1*gamma;

    return t1;

    //sanity check (checks out...):
    radial_basis_deriv_num(x, y, n, real(0.001), param, deriv1);
    for (dim_ta i=0; i<n; i++) printf(" %g", deriv[i]);
    printf("\n");
    for (dim_ta i=0; i<n; i++) printf(" %g", deriv1[i]);
    printf("\n");
  }

  template <class real>
  real sigmoid_basis_deriv(real *x, real *y, dim_ta n, void *param, real *deriv) {
    real *p2=(real *) param;
    real gamma=p2[0];
    real coef0=p2[1];
    real dot=linear_basis_deriv(x, y, n, NULL, deriv);
    real t1=tanh(gamma*dot+coef0);
    for (dim_ta i=0; i<n; i++) deriv[i]*=gamma*(1-t1*t1);
    return t1;
  }

  template <class real, class cls_t>
  svm2class<real, cls_t>::svm2class() {
    sv=NULL;
    nsv=0;
    coef=NULL;
    param=NULL;
    this->mat=NULL;
    this->b=NULL;
    this->id=-1;
    this->D=0;
    this->D1=0;
  }

  template <class real, class cls_t>
  svm2class<real, cls_t>::svm2class(char *modfile) {
    int err=init(modfile);
    if (err!=0) throw err;
  }

  template <class real, class cls_t>
  svm2class<real, cls_t>::~svm2class() {
    delete [] sv[0];
    delete [] sv;
    delete [] param;
    delete [] coef;
  }

  template <class real, class cls_t>
  int svm2class<real, cls_t>::init(char *modfile) {
    FILE *fs=fopen(modfile, "r");
    char *line=NULL;
    char **substr=NULL;
    int nsub;
    int nsv1;		//should agree with nsv...
    int pfound=0;	//number of parameters for probability estimation read in (should be 2)

    this->mat=NULL;
    this->b=NULL;
    this->id=-1;

    if (fs==NULL) {
      fprintf(stderr, "svm2class: failed to open model file, %s\n", modfile);
      return UNABLE_TO_OPEN_FILE_FOR_READING;
    }

    rho=0;
    param=new real[3];

    do {
      if (line!=NULL) delete [] line;
      if (substr!=NULL) delete [] substr;
      line=fget_line(fs, 1);
      substr=split_string_destructive(line, nsub);
      //printf("svm2class: nsub=%d\n", nsub);
      //for (int i=0; i<nsub; i++) printf("%s ", substr[i]);
      //printf("\n");
      if (strcmp(substr[0], "svm_type")==0) {
        if (strcmp(substr[1], "c_svc")!=0 && strcmp(substr[1], "nu-svc")!=0) {
          fprintf(stderr, "svm2class: not a classifier SVM in file, %s (%s)\n", modfile, substr[1]);
	  fclose(fs);
	  return PARAMETER_OUT_OF_RANGE;
        }
      } else if (strcmp(substr[0], "nr_class")==0) {
        this->ncls=atoi(substr[1]);
	if (this->ncls != 2) {
          fprintf(stderr, "svm2class: only binary classifiers accepted (file, %s)\n", modfile);
	  fclose(fs);
	  return PARAMETER_OUT_OF_RANGE;
        }
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
          fprintf(stderr, "svm2class: basis function, %s, not recognized (file, %s)\n", substr[1], modfile);
	  fclose(fs);
	  return PARAMETER_OUT_OF_RANGE;
        }
      } else if (strcmp(substr[0], "gamma")==0) {
        param[0]=atof(substr[1]);
      } else if (strcmp(substr[0], "coef0")==0) {
        param[1]=atof(substr[1]);
      } else if (strcmp(substr[0], "degree")==0) {
        param[2]=atof(substr[1]);
      } else if (strcmp(substr[0], "total_sv")==0) {
        nsv1=atoi(substr[1]);
      } else if (strcmp(substr[0], "rho")==0) {
        rho=atof(substr[1]);
      } else if (strcmp(substr[0], "probA")==0) {
        probA=atof(substr[1]);
	pfound++;
      } else if (strcmp(substr[0], "probB")==0) {
        probB=atof(substr[1]);
	pfound++;
      } else if (strcmp(substr[0], "label")==0) {
        label1=atoi(substr[1]);
        label2=atoi(substr[2]);
      }
    } while (strcmp(substr[0], "SV")!=0);
    delete [] line;
    delete [] substr;

    if (pfound!=2) {
      fprintf(stderr, "svm2class: probability estimates must be supported (in file, %s)\n", modfile);
      fclose(fs);
      return FILE_READ_ERROR;
    }

    nsv=read_svm(fs, sv, coef, this->D);

    if (nsv!=nsv1) {
      fprintf(stderr, "svm2class: warning, number of support vectors read in does not match stated (%d vs. %d) (file: %s)\n", nsv, nsv1, modfile);
    }

    this->D1=this->D;

    fclose(fs);
  }

  template <class real, class cls_t>
  real svm2class<real, cls_t>::R(real *x, real *praw) {
    real sum=0;
    real r;
    real *x1=this->do_xtran(x);

    for (nel_ta i=0; i<nsv; i++) {
      sum+=(*kernel)(x1, sv[i], this->D1, param)*coef[i];
    }
    sum-=rho;
    r=1-2/(1+exp(sum*probA+probB));

    if (this->id>=0 && praw!=NULL) praw[this->id]=sum;
    if (this->mat!=NULL) delete [] x1;
    printf("R: x=");
    for (dim_ta i=0; i<this->D; i++) printf(" %g", x[i]);
    printf("\n");

    return r;
  }

  //some book-keeping:
  template <class real, class cls_t>
  cls_t svm2class<real, cls_t>::class_list(cls_t *list) {
    list[0]=label1;
    list[1]=label2;
    return this->ncls;
  }

  template <class real, class cls_t>
  void svm2class<real, cls_t>::R_deriv_num(real *x, real dx, real *drdx) {
    real x1[this->D];
    real x2[this->D];
    real r1, r2;
    for (dim_ta i=0; i<this->D; i++) {
      x1[i]=x[i];
      x2[i]=x[i];
    }
    for (dim_ta i=0; i<this->D; i++) {
      x1[i]-=dx;
      x2[i]+=dx;
      r1=R(x1);
      r2=R(x2);
      //printf("r1=%g; r2=%g\n", r1, r2);
      //printf("x1[%d]=%g; x2[%d]=%g\n", i, x1[i], i, x2[i]);
      drdx[i]=(r2-r1)/dx/2;
      //printf("%g ", drdx[i]);
      x1[i]=x[i];
      x2[i]=x[i];
    }
    printf("\n");
    printf("R_d_num: x=");
    for (dim_ta i=0; i<this->D; i++) printf(" %g", x[i]);
    printf("\n");
  }

  template <class real, class cls_t>
  real svm2class<real, cls_t>::R_deriv(real *x, real *drdx) {
    real sum=0;
    real t1, t2;		//temporaries
    real deriv[this->D1];
    real drdx1[this->D1];
    real drdx2[this->D];
    real r;
    real *x1=this->do_xtran(x);

    for (dim_ta j=0; j<this->D1; j++) drdx1[j]=0;
    for (nel_ta i=0; i<nsv; i++) {
      sum+=coef[i]*(*kernel_deriv)(x1, sv[i], this->D1, param, deriv);
      for (dim_ta j=0; j<this->D1; j++) drdx1[j]+=coef[i]*deriv[j];
    }
    sum-=rho;

    t1=exp(sum*probA+probB);
    t2=probA*t1/(1+t1)/(1+t1);		//derivative of sigmoid function
    for (dim_ta j=0; j<this->D1; j++) drdx1[j]*=2*t2;
    r=1-2/(1+t1);

    if (this->mat!=NULL) {
      delete [] x1;
      //grad_x(f(Ax))=grad_y(f(y)*A
      vector_mult(this->mat, drdx1, drdx, this->D, this->D1);
    } else {
      for (dim_ta j=0; j<this->D; j++) drdx[j]=drdx1[j];
    }

    //basic sanity check:
    R_deriv_num(x, 0.00001, drdx2);
    printf("%g : ", r);
    //printf("%g : ", sum);
    for (dim_ta i=0; i<this->D; i++) printf(" %g", drdx[i]);
    printf("\n");
    printf("%g : ", R(x));
    for (dim_ta i=0; i<this->D; i++) printf(" %g", drdx2[i]);
    printf("\n");

    //return sum;
    return r;

  }

  template <class real, class cls_t>
  real svmrfunc(real *x, void *param, real *deriv) {
    bordparam<real> *p1=(bordparam<real> *) param;
    svm2class<real, cls_t> *p2=(svm2class<real, cls_t> *) p1->rparam;
    return p2->R_deriv(x, deriv);
  }

  template class svm2class<real_a, cls_ta>;

  template real_a svmrfunc<real_a, cls_ta>(real_a *, void *, real_a *);

}
