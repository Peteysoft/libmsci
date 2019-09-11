#include <string.h>
#include <math.h>

#include "read_ascii_all.h"
#include "error_codes.h"
//#include "agf_defs.h"
#include "agf_fconv.h"
#include "linearclassifier.h"

using namespace libpetey;

namespace libagf {

  template <typename real, typename cls_t>
  linearclassifier<real, cls_t>::linearclassifier() {
    this->ncls=2;
    this->calcoef=new real[2];
    this->calcoef[0]=0;
    this->calcoef[1]=0.5;
    this->order=1;
    header=NULL;
    this->sigmoid_func=NULL;
  }

  template <typename real, typename cls_t>
  linearclassifier<real, cls_t>::linearclassifier(char *fname) {
    FILE *fs;
    int err;

    //book-keeping:
    this->ncls=2;
    this->calcoef=new real[2];
    this->calcoef[0]=0;
    this->calcoef[1]=0.5;
    this->order=1;
    header=NULL;
    this->sigmoid_func=&tanh;

    fs=fopen(fname, "r");
    if (fs==NULL) {
      fprintf(stderr, "linearclassifier: Failed to open file, %s, for reading\n", fname);
      throw UNABLE_TO_OPEN_FILE_FOR_READING;
    }

    err=load(fs);
    if (err!=0) {
      fprintf(stderr, "linearclassifier: File read error: %s\n", fname);
      throw FILE_READ_ERROR;
    }
  }

  template <typename real, typename cls_t>
  linearclassifier<real, cls_t>::~linearclassifier() {
    if (header!=NULL) {
      for (int i=0; header[i]!=NULL; i++) delete [] header[i];
      delete [] header;
    }
    delete [] this->calcoef;
    delete [] coef;
  }

  template <typename real, typename cls_t>
  real linearclassifier<real, cls_t>::decision(real *x) {
    real result=offset;
    for (int i=0; i<this->D; i++) result+=x[i]*coef[i];
    return result;
  }

  template <typename real, typename cls_t>
  int linearclassifier<real, cls_t>::load(FILE *fs) {
    char *line=NULL;
    char **substr=NULL;
    int nsub;
    int label1, label2;
    double bias;
    double offset1;
    double *w;
    int lineno=0;
    int nhead=0;

    //this->mat=NULL;
    //this->b=NULL;
    //this->id=-1;

    header=new char *[100];

    do {
      if (line!=NULL) delete [] line;
      if (substr!=NULL) delete [] substr;
      line=fget_line(fs, 1);
      lineno++;
      header[nhead]=new char[strlen(line)+1];
      strcpy(header[nhead], line);
      nhead++;
      substr=split_string_destructive(line, nsub);
      //printf("svm2class: nsub=%d\n", nsub);
      //for (int i=0; i<nsub; i++) printf("%s ", substr[i]);
      //printf("\n");
      if (nsub == 0) continue;
      if (strcmp(substr[0], "w")!=0 && nsub<2) {
        fprintf(stderr, "linearclassifier: error in initialization file; unrecognized keywordi (%s)/not enough parameters\n", substr[0]);
        return FILE_READ_ERROR;
      }
      if (strcmp(substr[0], "nr_class")==0) {
        this->ncls=atoi(substr[1]);
	if (this->ncls!=2) {
          fprintf(stderr, "linearclassifier: only binary classifiers accepted\n");
	  return PARAMETER_OUT_OF_RANGE;
	}
      } else if (strcmp(substr[0], "label")==0) {
        label1=atoi(substr[1]);
        label2=atoi(substr[2]);
	nhead--;
	delete [] header [nhead];
      } else if (strcmp(substr[0], "nr_feature")==0) {
        this->D=atoi(substr[1]);
        this->D1=this->D;
	nhead--;
	delete [] header[nhead];
      } else if (strcmp(substr[0], "bias")==0) {
        bias=atof(substr[1]);
	nhead--;
	delete [] header[nhead];
      } else if (strcmp(substr[0], "w")==0) {
	nhead--;
	delete [] header[nhead];
        break;
      }
    } while (feof(fs)==0);
    if (line!=NULL) delete [] line;
    delete [] substr;
    header[nhead]=NULL;
    coef=new real[this->D];
    w=new double[this->D];
    for (int i=0; i<this->D; i++) {
      fscanf(fs, "%lg", w+i);		//should add some error catch stuff...
      coef[i]=w[i];
    }
    fscanf(fs, "%lg", &offset1);
    if (bias!=-1) offset=offset1*bias; else offset=0;
    if (label1==0 && label2==1) {
      for (int i=0; i<this->D; i++) coef[i]=-coef[i];
      offset=-offset;
    }
    delete [] w;
    return 0;
  }

  template <typename real, typename cls_t>
  int linearclassifier<real, cls_t>::save(FILE *fs) {
    char format[12];    //format code
    char fcode[4];      //floating point format code
    get_format_code<real>(fcode);
    sprintf(format, "%%%s", fcode);
    for (int i=0; header[i]!=NULL; i++) fprintf(fs, "%s\n", header[i]);
    fprintf(fs, "nr_class %d\n", this->ncls);
    fprintf(fs, "label 0 1\n");
    fprintf(fs, "nr_feature %d\n", this->D);
    fprintf(fs, "bias 1\n");
    fprintf(fs, "w\n");
    for (int i=0; i<this->D; i++) {
      fprintf(fs, format, coef[i]);
      fprintf(fs, "\n");
    }
    fprintf(fs, format, offset);
    return 0;
  }

  template class linearclassifier<float, cls_ta>;
  template class linearclassifier<double, cls_ta>;

} //end namespace libagf
