#include <assert.h>

#include "agf_lib.h"

namespace libagf {
  //multi-class classification adapted for continuum retrievals:
  template <typename real, typename cls_t>
  multiclass_c<real, cls_t>::multiclass_c(const char *file, int tp) {
    cls_t *clist;
    FILE *fs;
    int lineno=0;
    cls_t cnt=0;
    multi_parse_param param;

    fs=fopen(file, "r");
    if (fs==NULL) {
      fprintf(stderr, "multiclass_c: Unable to open control file, %s\n", file);
      exit(UNABLE_TO_OPEN_FILE_FOR_READING);
    }
    param.infs=fs;
    param.lineno=0;
    this->init(param);
    this->set_id(&cnt);

    //requirements are more stringent than with plain classifiers:
    if (this->ncls-1!=cnt) {
      fprintf(stderr, "multiclass_c: %d classes and %d binary classifiers;\n", this->ncls, cnt);
      fprintf(stderr, "  there must be at least one fewer binary classifiers as classes\n");
      exit(SAMPLE_COUNT_MISMATCH);
    }

    //we need class labels to go from 0..nc-1:
    clist=new cls_ta[this->ncls];
    this->class_list(clist);
    for (cls_ta i=0; i<this->ncls; i++) {
      if (clist[i]!=i) {
        fprintf(stderr, "multiclass_c: class labels must go from 0..nc-1\n");
        exit(PARAMETER_OUT_OF_RANGE);
      }
    }

    //now we read the abscissa from the control file:
    y0=new real_a[cnt];
    for (cls_ta i=0; i<cnt; i++) {
      int err=fscanf(fs, "%g", y0+i);
      printf("%d %g\n", i, y0[i]);
      if (err!=1) {
        fprintf(stderr, "multiclass_c: Error reading abscissa from control file, line %d\n", lineno);
        exit(FILE_READ_ERROR);
      }
      if (i>0) {
        if (y0[i]<=y0[i-1]) {
          fprintf(stderr, "multiclass_c: Abscissa must be in ascending order, line %d\n", lineno);
          exit(PARAMETER_OUT_OF_RANGE);
        }
      }
      lineno++;
    }
    fclose(fs);
  }


  template <typename real, typename cls_t>
  multiclass_c<real, cls_t>::~multiclass_c() {
    delete [] y0;
  }


  //continuum retrieval:
  template <typename real, typename cls_t>
  real multiclass_c<real, cls_t>::ret(real *x, real &err) {
    real a12, a22, b1, b2;		//elements of normal equation
    real det;				//determinant
    cls_t cls;
    real result;
    real raw[(this->ncls-1)*2];

    for (cls_t j=0; j<this->ncls; j++) raw[j]=NAN;

    if (this->nonh_flag && this->ncls > 2) {
      real p[this->ncls];
      real *gd=raw+this->ncls-1;

      //some waste here, since we use neither the conditional probabilities,
      //nor the class result:
      cls=this->classify_t(x, p, raw);

      for (int i=0; i<this->ncls; i++) printf("%g ", raw[i]);
      printf("\n");

      //least squares by hand:
      a12=0;
      a22=0;
      b1=0;
      b2=0;
      for (int i=0; i<this->ncls-1; i++) {
        //real gd2=gd[i]*gd[i];
        real gd2=1;
        //let the compiler optimizer take car of this?
        a12+=gd2/raw[i];
        a22+=gd2/raw[i]/raw[i];
        b1+=gd2*y0[i]/raw[i];
        b2+=gd2*y0[i]/raw[i]/raw[i];
      }
      a12=-M_SQRT2*a12;
      b1=-M_SQRT2*b1;
      
      //solve the least squares problem:
      det=2*(this->ncls-1)*a22-a12*a12;
      err=(a22*b1-a12*b2)/det;
      result=(-a12*b1+2*(this->ncls-1)*b2)/det;
    } else {
      cls_t l1, l2;
      real p;
      real gd2=1;
      int ng=0;

      cls=this->classify_t(x, p, raw);

      //least squares by hand:
      a12=0;
      a22=0;
      b1=0;
      b2=0;
      for (int i=0; i<this->ncls-1; i++) {
        if (isfinite(raw[i])!=1) continue;
        //real gd2=gd[i]*gd[i];
        printf(" %g", raw[i]);
        real gd2=1;
        //let the compiler optimizer take car of this?
        a12+=gd2/raw[i];
        a22+=gd2/raw[i]/raw[i];
        b1+=gd2*y0[i]/raw[i];
        b2+=gd2*y0[i]/raw[i]/raw[i];
        ng++;
      }
      printf("\n");
      a12=-M_SQRT2*a12;
      b1=-M_SQRT2*b1;
      
      //solve the least squares problem:
      det=2*ng*a22-a12*a12;
      err=(a22*b1-a12*b2)/det;
      result=(-a12*b1+2*ng*b2)/det;

/*
      if (cls==this->ncls-1) {
        l1=cls-2;
        l2=cls-1;
      } else if (cls==0) {
        l1=0;
        l2=1;
      } else {
        l1=cls-1;
        l2=cls;
      }

      for (int i=0; i<this->ncls; i++) {
        printf(" %g", raw[i]);
        if (i==l1) printf("*");
        if (i==l2) printf("*");
      }
      printf("\n");

      assert(isfinite(raw[l1]));
      assert(isfinite(raw[l2]));
      //printf("%g*%g-%g*%g\n", raw[l2], y0[l1], raw[l1], y0[l2]);
      result=(raw[l2]*y0[l1]-raw[l1]*y0[l2])/(raw[l2]-raw[l1]);
      //result=cls;
      err=(result-y0[l1])/M_SQRT2/raw[l1];
*/
    }
    if (fabs(err) > 2) printf("FLAG\n");

    return result;

  }

  template class multiclass_c<real_a, cls_ta>;


}

