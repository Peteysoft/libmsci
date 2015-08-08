#include <math.h>

#include "read_ascii_all.h"
#include "error_codes.h"
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
    real dot=linear_basis(x, y, n, param);
    return pow(gamma*dot+coef0, degree);
  }

  template <class real>
  real polynomial_basis(real *x, real *y, dim_ta n, void *param) {
    real gamma=((real *) param)[0];
    real d2=agf_metric2(x, y, n);
    return exp(-gamma*d2);
  }

  template <class real>
  real sigmoid_basis(real *x, real *y, dim_ta n, void *param) {
    real *p2=(real *) param;
    real gamma=p2[0];
    real coef0=p2[1];
    real dot=linear_basis(x, y, n, param);
    return tanh(gamma*dot+coef0);
  }

  template <class real>
  svm2class<real>::svm2class(char *modfile) {
    int err=init(modfile);
    if (err!=0) throw err;
  }

  template <class real>
  int svm2class<real>::init(char *modfile) {
    FILE fs=fopen(modfile, "r");
    char *line;
    char *substr;
    int nsub;
    int nsv1;		//should agree with nsv...
    int pfound=0;	//number of parameters for probability estimation read in (should be 2)

    if (fs==NULL) {
      fprintf(stderr, "svm2class: failed to open model file, %s\n", modfile);
      return UNABLE_TO_OPEN_FILE_FOR_READING;
    }

    rho=0;
    param=new real[3];

    do {
      line=fget_line(fs, 1);
      substr=split_string_destructive(line, nsub);
      if (strcmp(substr[0], "svm_type")==0) {
        if (strcmp(substr[1], "c_svc")!=0 && strcmp(substr[1], "nu-svc")!=0) {
          fprintf(stderr, "svm2class: not a classifier SVM in file, %s\n", modfile);
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
	  kernel=&linear_basis;
	  param=NULL;
        } else if (strcmp(substr[1], "polynomial")==0) {
          kernel=&polynomial_basis;
	} else if (strcmp(substr[1], "rbf")==0) {
          kernel=&radial_basis;
	} else if (strcmp(substr[1], "sigmoid")==0) {
          kernel=&radial_basis;
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
      }
      delete [] line;
      delete [] substr;
    } while (strcmp(substr[0], "SV")!=0);

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
    real sum;
    real r;
    real *x1=do_ltran(x);

    for (nel_ta i=0; i<nsv; i++) {
      sum+=(*kernel)(x, sv[i], this->D1)*coef[i];
    }
    sum-=rho;
    r=1-2/(1+exp(sum*probA+probB));

    if (this->id>=0 && praw!=NULL) praw[this->id]=sum;
    if (this->mat!=NULL) delete [] x1;

    return r;

  }
