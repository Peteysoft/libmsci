#include <stdio.h>
#include <getopt.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <assert.h>

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
    //if (mat!=NULL) delete [] xtran;
  }

  //flag is ignored
  template <class real, class cls_t>
  int binaryclassifier<real, cls_t>::ltran(real **mat1, real *b1, dim_ta d1, dim_ta d2, int flag) {
    int err2=0;

    //err2=read_stats(normfile, ave, std, D);
    mat=mat1;
    b=b1;
 
    //we don't know offhand what the dimensionality of the problem is: 
    /*if (d2!=this->D) {
      fprintf(stderr, "agf2class: second dimension of trans. mat. does not that of borders data: %d vs. %d\n", d2, this->D);
      return DIMENSION_MISMATCH;
    }
    //this is very clear:
    D1=this->D;
    */

    D1=d2;
    this->D=d1;

    //from the outside, the classifier looks like it has the same number of
    //features as before normalization:
    //xtran=new real[D1];

    return err2;
  }

  template <class real, class cls_t>
  real * binaryclassifier<real, cls_t>::do_xtran(real *x) {
    real *xtran;
    if (mat!=NULL) {
      real tmp;
      xtran=new real[D1];
      //linearly transform the test point:
      for (dim_ta j=0; j<D1; j++) xtran[j]=0;
      for (dim_ta i=0; i<this->D; i++) {
        tmp=x[i]-b[i];
        for (dim_ta j=0; j<D1; j++) xtran[j]=xtran[j]+tmp*mat[i][j];
      }
      //xtran=left_vec_mult(x, mat, this->D, D1);
    } else {
      xtran=x;
    }
    return xtran;
  }

  template <class real, class cls_t>
  real binaryclassifier<real, cls_t>::R(real *x, real *praw) {
    //want to make this do a direct classifications using the options listed
    //in name!
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
  }

  template <class real, class cls_t>
  real general2class<real, cls_t>::R(real *x, real *praw) {
    real r;
    real *xtran;

    xtran=this->do_xtran(x);

    batchR(&xtran, &r, 1, this->D1);
    if (praw!=NULL && this->id>=0) praw[this->id]=r;

    if (this->mat!=NULL) delete [] xtran;

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

    real *xtran;		//use local variable for transformed x
				//to keep things thread-safe

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
      xtran=this->do_xtran(x[i]);
      cls[i]=0;			//clear class data
      //we write to the input file:
      if (Mflag) {
        sprintf(format, " %%d:%%.12%s", fcode);
        fprintf(fs, "%d", cls[i]);
        for (dim_ta j=0; j<nvar; j++) fprintf(fs, format, j+1, xtran[j]);
        fprintf(fs, "\n");
      } else {
        sprintf(format, "%%.12%s ", fcode);
        for (dim_ta j=0; j<nvar; j++) fprintf(fs, format, xtran[j]);
        fprintf(fs, "%d\n", cls[i]);
      }
      if (this->mat!=NULL) delete [] xtran;
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
    fprintf(param.commandfs, "%s -O \"%s\" %s %s %s\n", param.commandname, command, 
		this->name, param.train, fbase);

    //no need to print the class partions...
    return 2;
  }

  //normflag: pick up normalization matrix
  //uflag: border samples are stored un-normalized
  template <class real, class cls_t>
  agf2class<real, cls_t>::agf2class(const char *fbase) {
    //ave=NULL;
    this->mat=NULL;
    //this->xtran=NULL;
    this->id=-1;

    this->name=new char[strlen(fbase)+1];
    strcpy(this->name, fbase);

    int err=agf_read_borders(fbase, brd, grd, n, this->D);
    if (err!=0) exit(err);

    fprintf(stderr, "agf2class: %d border samples found in model, %s\n", n, fbase);
    this->D1=this->D;
 
    this->ncls=2;

    //calculate and store the lengths of all the gradient vectors for possible
    //use later on:
    gd=new real[n];
    for (nel_ta i=0; i<n; i++) {
      gd[i]=0;
      for (dim_ta j=0; j<this->D; j++) gd[i]+=grd[i][j]*grd[i][j];
      gd[i]=sqrt(gd[i]);
    }
  }

  template <class real, class cls_t>
  agf2class<real, cls_t>::~agf2class() {
    delete_matrix(brd);
    delete_matrix(grd);
  }

  //if flag, then border vectors are stored un-normalized
  template <class real, class cls_t>
  int agf2class<real, cls_t>::ltran(real **mat1, real *b1, dim_ta d1, dim_ta d2, int flag) {
    int err2=0;

    real **brd2;
    real **grd2;

    //err2=read_stats(normfile, ave, std, D);
    this->mat=mat1;
    this->b=b1;
  
    if (flag) {
      if (d1!=this->D) {
        fprintf(stderr, "agf2class: first dimension (%d) of trans. mat. does not agree with that of borders data (%d)\n", d1, this->D);
        return DIMENSION_MISMATCH;
      }
      fprintf(stderr, "agf2class: Normalising the border samples...\n");
      //from the outside, the classifier looks like it has the same number of
      //features as before normalization:
      this->D1=d2;

      //apply constant factor:
      for (nel_ta i=0; i<n; i++) {
        for (dim_ta j=0; j<this->D1; j++) {
          brd[i][j]=brd[i][j]-this->b[j];
        }
      }

      brd2=matrix_mult(brd, this->mat, n, this->D, this->D1);
      grd2=matrix_mult(grd, this->mat, n, this->D, this->D1);
      delete_matrix(brd);
      delete_matrix(grd);
      brd=brd2;
      grd=grd2;
      //have to recalculate all the gradient vectors, dammit:
      for (nel_ta i=0; i<n; i++) {
        gd[i]=0;
        for (dim_ta j=0; j<this->D1; j++) gd[i]+=grd[i][j]*grd[i][j];
        gd[i]=sqrt(gd[i]);
      }
    } else {
      if (d2!=this->D) {
        fprintf(stderr, "agf2class: second dimension of trans. mat. does not that of borders data: %d vs. %d\n", d2, this->D);
        return DIMENSION_MISMATCH;
      }
      //this is very clear:
      this->D1=this->D;
      this->D=d1;
      //from the outside, the classifier looks like it has the same number of
      //features as before normalization:
    }
    //this->xtran=new real[this->D1];
    return err2;
  }

  template <class real, class cls_t>
  real agf2class<real, cls_t>::R(real *x, real *praw) {
    real r;
    real *xtran;

    xtran=this->do_xtran(x);
    if (this->id>=0 && praw!=NULL) {
      nel_ta k;
      real d;
      praw[this->id]=border_classify0(brd, grd, this->D1, n, xtran, k, d);
      //praw[this->id+n]=gd[k];//ain't gonna work...
      r=tanh(praw[this->id]);
      printf("r=%g\n" , praw[this->id]);
    } else {
      r=border_classify(brd, grd, this->D1, n, xtran);
    }
    if (this->mat!=NULL) delete [] xtran;
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

  template class binaryclassifier<real_a, cls_ta>;
  template class general2class<real_a, cls_ta>;
  template class agf2class<real_a, cls_ta>;

}

