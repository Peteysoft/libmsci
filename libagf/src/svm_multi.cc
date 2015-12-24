#include "read_ascii_all.h"
#include "svm_multi.h"

using namespace libpetey;

namespace libagf {
  template <class real, class cls_t>
  svm_multi<real, cls_t>::svm_multi() {
    this->ncls=0;
    this->D=0;
    nsv_total=0;
    sv=NULL;
    coef=NULL;
    rho=NULL;
    probA=NULL;
    probB=NULL;
    label=NULL;
  }

  template <class real, class cls_t>
  svm_multi<real, cls_t>::svm_multi(FILE *file) {
    FILE *fs=fopen(file, "r");
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
    
    probA=NULL;
    probB=NULL;

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
      lineno++;
      substr=split_string_destructive(line, nsub);
      //printf("svm2class: nsub=%d\n", nsub);
      //for (int i=0; i<nsub; i++) printf("%s ", substr[i]);
      //printf("\n");
      if (strcmp(substr[0], "svm_type")==0) {
        if (strcmp(substr[1], "c_svc")!=0 && strcmp(substr[1], "nu-svc")!=0) {
          fprintf(stderr, "svm_multi: not a classifier SVM in file, %s (%s)\n", modfile, substr[1]);
	  fclose(fs);
	  throw PARAMETER_OUT_OF_RANGE;
        }
      } else if (strcmp(substr[0], "nr_class")==0) {
        this->ncls=atoi(substr[1]);
	if (this->ncls < 2) {
          fprintf(stderr, "svm_multi: one class classifiers not accepted (file, %s)\n", modfile);
	  fclose(fs);
	  return PARAMETER_OUT_OF_RANGE;
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
          fprintf(stderr, "svm2class: basis function, %s, not recognized (file, %s)\n", substr[1], modfile);
	  fclose(fs);
	  throw PARAMETER_OUT_OF_RANGE;
        }
      } else if (strcmp(substr[0], "gamma")==0) {
        param[0]=atof(substr[1]);
      } else if (strcmp(substr[0], "coef0")==0) {
        param[1]=atof(substr[1]);
      } else if (strcmp(substr[0], "degree")==0) {
        param[2]=atof(substr[1]);
      } else if (strcmp(substr[0], "total_sv")==0) {
        nsv_total=atoi(substr[1]);
      } else if (strcmp(substr[0], "rho")==0) {
        if (nsub<nparam+1) {
          fprintf(stderr, "svm2class: error in initialization file: not enough parameters (rho) (file, %s)\n", modfile);
	  fclose(fs);
	  throw FILE_READ_ERROR;
	}
        rho=new real[nparam];
	for (int i=0; i<nparam; i++) rho[i]=atof(substr[i+1]);
      } else if (strcmp(substr[0], "probA")==0) {
        if (nsub<nparam+1) {
          fprintf(stderr, "svm2class: error in initialization file: not enough parameters (probA) (file, %s)\n", modfile);
	  fclose(fs);
	  throw FILE_READ_ERROR;
	}
        probA=new real[nparam];
	for (int i=0; i<nparam; i++) probA[i]=atof(substr[i+1]);
	pfound++;
      } else if (strcmp(substr[0], "probB")==0) {
        if (nsub<nparam+1) {
          fprintf(stderr, "svm2class: error in initialization file: not enough parameters (probB) (file, %s)\n", modfile);
	  fclose(fs);
	  throw FILE_READ_ERROR;
	}
        probB=new real[nparam];
	for (int i=0; i<nparam; i++) probB[i]=atof(substr[i+1]);
	pfound++;
      } else if (strcmp(substr[0], "label")==0) {
        if (nsub<this->ncls+1) {
          fprintf(stderr, "svm2class: error in initialization file: not enough parameters (label) (file, %s)\n", modfile);
	  fclose(fs);
	  throw FILE_READ_ERROR;
	}
        label=new real[this->ncls];
	for (int i=0; i<nparam; i++) label[i]=atoi(substr[i+1]);
      }
    } while (strcmp(substr[0], "SV")!=0);
    delete [] line;
    delete [] substr;

    dim_ta **ind=new dim_ta *[nsv_total];	//dimension indices
    real **raw=new real *[nsv_total];		//raw features data
    int nf[nsv_total];
    coef=new real*[nsv_total];
    coef[0]=new real[nsv_total*nparam];
    this->D=0;
    for (int i=0; i<nsv_total; i++) {
      int nread;		//number of item scanned
      int pos;			//position in line
      int rel;			//relative position in line
      int nf;			//number of features read in
      line=read_line(fs);
      coef[i]=coef[0]+i*nparam;
      for (int j=0; j<nparam; j++) {
        nread=sscanf(line+pos, format, coeff[i]+j, &rel);
	if (nread!=1) {
          fprintf(stderr, "svm_multi: error reading coefficients from %s line %d\n", file, lineno+i);
	  throw FILE_READ_ERROR;
	}
	pos+=rel;
      }
      nf[i]=scan_svm_features(line+pos, ind[i], raw[i]);
      if (nf[i]<=0) {
        fprintf(stderr, "svm_multi: error reading support vectors from %s line %d\n", file, lineno+i);
	throw FILE_READ_ERROR;
      }
      for (int j=0; j<nf[i]; j++) if (ind[i]>this->D) this->D=ind[i];
    }
    //transfer to more usual array and fill in missing values:
    real missing=0;
    sv=new real*[nsv_total];
    sv[0]=new real[nsv_total*this->D];
    for (int i=0; i<nsv_total; i++) {
      sv[i]=sv[0]+i*this->D;
      for (int j=0; j<this->D; j++) sv[i][j]=missing;
      for (int j=0; j<nf[i]; j++) sv[i][dim[i][j]-1]=raw[i][j];
      delete [] dim[i];
      delete [] raw[i];
    }
    delete [] dim;
    delete [] raw;
    delete [] line;

    fclose(fs);
  }

  template <class real, class cls_t>
  svm_multi<real, cls_t>::svm_multi() {
    delete [] sv;
    delete [] coef;
    delete [] rho;
    if (probA!=NULL) delete [] probA;
    if (probB!=NULL) delete [] probB;
    delete [] label;
  }
