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
  binaryclassifier<real, cls_t> * binclass_init(char *name, int typecode) {
    binaryclassifier<real, cls_t> *result;
    switch (typecode) {
      case(0):
        result=new borders_classifier<real, cls_t>(name);
	break;
      case(1):
        result=new svm2class<real, cls_t>(name);
	break;
      case(2):
        result=new borders_calibrated<real, cls_t>(name);
	break;
      case(3):
        result=new svm2class2<real, cls_t>(name);
	break;
      case(4):
        result=new binaryclassifier<real, cls_t>(name);
	break;
      default:
        result=new borders_classifier<real, cls_t>(name);
	break;
    }
    return result;
  }

  template <class real, class cls_t>
  binaryclassifier<real, cls_t>::binaryclassifier() {
    id=-1;
    this->ncls=2;
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
      fprintf(fs, "%s", this->name);
    } else {
      fprintf(fs, "%s", fbase);
    }
  }

  template <class real, class cls_t>
  void binaryclassifier<real, cls_t>::R_deriv_num(real *x, real dx, real *drdx) {
    real x1[this->D1];
    real x2[this->D1];
    real r1, r2;
    for (dim_ta i=0; i<this->D1; i++) {
      x1[i]=x[i];
      x2[i]=x[i];
    }
    for (dim_ta i=0; i<this->D1; i++) {
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
    for (dim_ta i=0; i<this->D1; i++) printf(" %g", x[i]);
    printf("\n");
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

  //does almost the same as the previous function...
  template <class real, class cls_t>
  cls_t binaryclassifier<real, cls_t>::collect_binary_classifiers(binaryclassifier<real, cls_t> **list) {
    list[0]=this;
    return 1;
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
    if (Mflag) {
      sprintf(format, " %%d:%%.12%s", fcode);
    } else {
      fprintf(fs, "%d\n", nvar);
      sprintf(format, "%%.12%s ", fcode);
    }
    for (nel_ta i=0; i<n; i++) {
      cls[i]=0;			//clear class data
      //we write to the input file:
      if (Mflag) {
        fprintf(fs, "%d", cls[i]);
        for (dim_ta j=0; j<nvar; j++) fprintf(fs, format, j+1, x[i][j]);
        fprintf(fs, "\n");
      } else {
        for (dim_ta j=0; j<nvar; j++) fprintf(fs, format, x[i][j]);
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
      throw err;
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
      throw FILE_READ_ERROR;
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
      throw PARAMETER_OUT_OF_RANGE;
    }

    if (err!=n) {
      fprintf(stderr, "general2class: error reading output file, %s, exiting\n", outfile);
      throw FILE_READ_ERROR;
    }
      
  }

  template <class real, class cls_t>
  int general2class<real, cls_t>::commands(multi_train_param &param, 
		cls_t **clist, char *fbase) {
    char Kstr[3];
    if (Kflag) {
      Kstr[0]='-';
      Kstr[1]='K';
      Kstr[2]='\0';
    } else {
      Kstr[0]='\0';
    }
    //accelerator mode:
    if (strcmp(command, "")==0) {
      //if command name is the empty string, we don't include it:
      fprintf(param.commandfs, "%s %s %s %s %s", param.commandname, Kstr, this->name, 
		      param.train, fbase);
    } else {
      fprintf(param.commandfs, "%s %s -O \"%s\" %s %s %s", param.commandname, Kstr, command, 
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

  template class binaryclassifier<float, cls_ta>;
  template class general2class<float, cls_ta>;

  template class binaryclassifier<double, cls_ta>;
  template class general2class<double, cls_ta>;

  template binaryclassifier<float, cls_ta> *binclass_init<float, cls_ta>(char *, int);
  template binaryclassifier<double, cls_ta> *binclass_init<double, cls_ta>(char *, int);

}

