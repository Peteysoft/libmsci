#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <assert.h>

#include "linked.h"
#include "full_util.h"
#include "agf_lib.h"

using namespace libpetey;

namespace libagf {

  template <class real, class cls_t>
  direct_classifier<real, cls_t>::direct_classifier() {
    train=NULL;
    cls=NULL;
    this->mat=NULL;
    this->b=NULL;
    this->name=NULL;
  }

  template <class real, class cls_t>
  direct_classifier<real, cls_t>::direct_classifier(cls_t nc, char *opts, char t, int mf) {
    this->ncls=nc;
    train=NULL;
    cls=NULL;
    this->mat=NULL;
    this->b=NULL;
    this->name=NULL;
    this->options=new char[strlen(opts)+1];
    strcpy(this->options, opts);
    type=t;
    Mflag=mf;
  }

  template <class real, class cls_t>
  direct_classifier<real, cls_t>::~direct_classifier() {
    if (train!=NULL) delete [] train;
    if (cls!=NULL) delete [] cls;
    if (this->name != NULL) delete [] this->name;
    if (options != NULL) delete [] this->options;
  }

  template <class real, class cls_t>
  int direct_classifier<real, cls_t>::ltran_model(real **mat, real *b, dim_ta d1,
			dim_ta d2) {
    int err=0;
    real **train2;

    if (d1!=this->D) {
      fprintf(stderr, "direct_classifier: first dimension (%d) of trans. mat. does not agree with that of borders data (%d)\n", d1, this->D);
      return DIMENSION_MISMATCH;
    }
    fprintf(stderr, "direct_classifier: Normalising the border samples...\n");
    if (this->mat==NULL) {
      this->D=d2;
      this->D1=d2;
    } else {
      assert(mat==this->mat && b==this->b);
      //from the outside, the classifier looks like it has the same number of
      //features as before normalization:
      assert(this->D1==d2);
    }

    //apply constant factor:
    for (nel_ta i=0; i<ntrain; i++) {
      for (dim_ta j=0; j<this->D1; j++) {
        train[i][j]=train[i][j]-b[j];
      }
    }

    train2=matrix_mult(train, this->mat, ntrain, this->D, this->D1);
    delete_matrix(train);
    train=train2;

    return err;
  }

  template <class real, class cls_t>
  int direct_classifier<real, cls_t>::init(real **x, cls_t *c, nel_ta nsamp, dim_ta ndim) {
    train=x;
    cls=c;
    ntrain=nsamp;
    this->D=ndim;
    this->ncls=1;
    for (nel_ta i=0; i<ntrain; i++) if (cls[i]>this->ncls) this->ncls=cls[i]+1;
    return 0;
  }

  template <class real, class cls_t>
  void direct_classifier<real, cls_t>::print(FILE *fs, char *fbase, int depth) {
    for (int i=0; i<depth; i++) fprintf(fs, "  ");
    if (fbase==NULL) {
      fprintf(fs, "%s %c \"%s\"", this->name, type, options);
    } else {
      fprintf(fs, "%s %c \"%s\"", fbase, type, options);
    }
  }

  template <class real, class cls_t>
  int direct_classifier<real, cls_t>::commands(multi_train_param &param, 
		  cls_t **clist, char *fbase) {
    //the commands consist solely of partitioning the classes and writing it
    //to a file:
    if (param.partcom==NULL) {
      //print command name etc.:
      if (type=='G') {
        char mstr[3];
        if (Mflag) {
          strcpy(mstr, "-M");
        } else {
          mstr[0]='\0';
        }
        //print command name etc.:
        fprintf(param.commandfs, "%s%s%s %s -A -1", AGF_COMMAND_PREFIX, AGF_PARTCOM,
		AGF_OPT_VER, mstr);
      } else {
        fprintf(param.commandfs, "%s%s%s", AGF_COMMAND_PREFIX, AGF_PARTCOM,
		AGF_OPT_VER);
      }
    } else {
      fprintf(param.commandfs, "%s", param.partcom);
    }

    if (type=='G') {
      fprintf(param.commandfs, " %s", param.train);
    } else {
      fprintf(param.commandfs, " %s %s", param.train, fbase);
    }

    //print class partitions:
    for (cls_t i=0; clist[0]+i!=clist[1]; i++) {
      fprintf(param.commandfs, " %d", clist[0][i]);
    }
    for (cls_t k=1; k<this->ncls; k++) {
      fprintf(param.commandfs, " %c", PARTITION_SYMBOL);
      for (cls_t i=0; clist[k]+i!=clist[k+1]; i++) {
        fprintf(param.commandfs, " %d", clist[k][i]);
      }
    }

    if (type=='G') {
      //if there is a command for file conversion, pipe to that first, then to the file:
      if (param.concom!=NULL) {
        fprintf(param.commandfs, " | %s > %s", param.concom, fbase);
      } else {
        fprintf(param.commandfs, "> %s", fbase);
      }
    }

    fprintf(param.commandfs, "\n");

    return this->ncls;
  }

  template <class real, class cls_t>
  agf_classifier<real, cls_t>::agf_classifier(const char *fbase, const char *options) {
    linked_list<char *> args;
    char **argv;
    char *ocpy;
    long argc;
    int argc2;
    int wsflag=0;
    agf_command_opts opt_args;
    int err;

    logfs=stderr;
    this->type='A';

    this->name=new char[strlen(fbase)+1];
    strcpy(this->name, fbase);
    this->options=new char[strlen(options)+1];
    strcpy(this->options, options);

    err=agf_read_train(fbase, this->train, this->cls, this->ntrain, this->D);
    //fprintf(stderr, "agf_classifier: %d sample found in model %s\n", this->ntrain, fbase);
    if (err!=0) exit(err);

    //count number of classes
    for (nel_ta i=0; i<this->ntrain; i++) {
      if (this->cls[i]>this->ncls) this->ncls=this->cls[i]+1;
    }

    //parse the options:
    ocpy=new char [strlen(options)+1];
    strcpy(ocpy, options);

    if (isspace(ocpy[0])) {
      wsflag=1;
    } else {
      wsflag=0;
      args.add(ocpy);
    }

    for (int i=1; ocpy[i]!='\0'; i++) {
      if (wsflag && isspace(ocpy[i])==0) {
        args.add(ocpy+i);
	wsflag=0;
      } else if (wsflag==0 && isspace(ocpy[i])) {
        ocpy[i]='\0';
	wsflag=1;
      }
    }

    argv=args.make_array(argc);
    argc2=argc;

    //printf("agf_classifier: %d options found;\n", argc2);
    //for (int i=0; i<argc2; i++) printf(" %s", argv[i]);
    //printf("\n");

    opt_args.k=K_DEFAULT_AGF;
    opt_args.W2=W_DEFAULT;
    err=agf_parse_command_opts(argc2, argv, "k:v:V:W:I:N:", &opt_args);
    if (err!=0) exit(err);

    k=opt_args.k;
    W=opt_args.W2;
    var[0]=opt_args.var[0];
    var[1]=opt_args.var[1];

    delete [] ocpy;
    delete [] argv;

    //check parameter ranges:
    if (var[0] <= 0 || var[1] <= 0) {
      //calculate the averages and standard deviations:
      real_a std[this->D];
      real_a ave[this->D];
      real_a vart;

      calc_norm(this->train, this->D, this->ntrain, ave, std);

      //if the initial filter variance is not set, set it to the total
      //variance of the data:
      vart=0;
      for (dim_ta i=0; i<this->D; i++) vart+=std[i]*std[i];
      if (var[0] <= 0) {
        var[0]=vart/pow(this->ntrain, 2./this->D);
        fprintf(logfs, "Using %10.3g for lower filter variance bracket\n\n", var[0]);
      }
      if (var[1] <= 0) {
        var[1]=vart;
        fprintf(logfs, "Using %10.3g for upper filter variance bracket\n\n", var[1]);
      }

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

    //check the range of k:
    if (k <= W || k >= this->ntrain) {
      if (k != -1) {
        fprintf(stderr, "agf: Parameter k=%d out of range.  Using all the training data.\n", k);
        k=-1;
        //err=PARAMETER_OUT_OF_RANGE;
      }
    }
  }

  template <class real, class cls_t>
  agf_classifier<real, cls_t>::~agf_classifier() {
  }

  template <class real, class cls_t>
  cls_t agf_classifier<real, cls_t>::classify(real *x, real *p, real *praw) {
    cls_t c;
    agf_diag_param diag;

    if (k==-1) {
      c=agf_classify(this->train, this->D1, this->cls, this->ntrain,
		      this->ncls, x, var, W, p, &diag);
    } else {
      c=agf_classify(this->train, this->D1, this->cls, this->ntrain, 
		      this->ncls, x, var, k, W, p, &diag);
      if (ntrial_k==0) min_f=diag.f;
      if (diag.f < min_f) min_f=diag.f;
      else if (diag.f > max_f) max_f=diag.f;
      total_f+=diag.f;
      ntrial_k++;
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

    return c;
  }

  template <class real, class cls_t>
  knn_classifier<real, cls_t>::knn_classifier(const char *fbase, const char *options) {
    linked_list<char *> args;
    char **argv;
    char *ocpy;
    long argc;
    int argc2;
    int wsflag=0;
    agf_command_opts opt_args;
    int err;

    logfs=stderr;
    this->type='K';

    this->name=new char[strlen(fbase)+1];
    strcpy(this->name, fbase);
    this->options=new char[strlen(options)+1];
    strcpy(this->options, options);

    err=agf_read_train(fbase, this->train, this->cls, this->ntrain, this->D);
    if (err!=0) exit(err);

    //count number of classes
    for (nel_ta i=0; i<this->ntrain; i++) {
      if (this->cls[i]>this->ncls) this->ncls=this->cls[i]+1;
    }

    //parse the options:
    ocpy=new char [strlen(options)+1];
    strcpy(ocpy, options);

    if (isspace(options[0])) {
      wsflag=1;
    } else {
      wsflag=0;
      args.add(ocpy);
    }

    for (int i=1; ocpy[i]!='\0'; i++) {
      if (wsflag && isspace(ocpy[i])==0) {
        args.add(ocpy+i);
	wsflag=0;
      } else if (wsflag==0 && isspace(ocpy[i])) {
        ocpy[i]='\0';
	wsflag=1;
      }
    }

    argv=args.make_array(argc);
    argc2=argc;

    opt_args.k=K_DEFAULT_KNN;
    err=agf_parse_command_opts(argc2, argv, "k:m:", &opt_args);
    if (err!=0) exit(err);

    k=opt_args.k;

    delete [] argv;
    delete [] ocpy;

    //check the range of k:
    if (k < 1 || k >= this->ntrain) {
      if (k != -1) {
        fprintf(stderr, "knn: Parameter k=%d out of range.  Using default.\n", k);
        k=K_DEFAULT_KNN;
        //err=PARAMETER_OUT_OF_RANGE;
      }
    }
  }

  template <class real, class cls_t>
  knn_classifier<real, cls_t>::~knn_classifier() {
  }

  template <class real, class cls_t>
  cls_t knn_classifier<real, cls_t>::classify(real *x, real *p, real *praw) {
    cls_t c;
    c=knn(global_metric2, this->train, this->D1, this->ntrain, this->cls, this->ncls, 
		    x, k, p);
    return c;
  }

  template <class real, class cls_t>
  general_classifier<real, cls_t>::general_classifier(const char *com, const char *file, cls_t nc, int mf, int kf) { 
    this->type='G';
    //store the command in the "options" field:
    this->options=new char[strlen(com)+1];
    strcpy(this->options, com);
    this->name=new char[strlen(file)+1];
    strcpy(this->name, file);
    this->ncls=nc;
    this->Mflag=mf;
    Kflag=kf;
  }

  template <class real, class cls_t>
  general_classifier<real, cls_t>::~general_classifier() { 
    //delete [] command;
  }

  template <class real, class cls_t>
  cls_t general_classifier<real, cls_t>::classify(real *x, real *p, real *praw) {
    cls_t c;
    batch_classify(&x, &c, &p, 1, this->D1);
    return c;
  }

  template <class real, class cls_t>
  void general_classifier<real, cls_t>::batch_classify(real **x, cls_t *cls, real **p, nel_ta n, dim_ta nvar) {
    FILE *fs;
    unsigned long uid;
    cls_t ncls;
    int err;

    char infile[25];
    char outfile[25];
    char *syscall;

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
    if (this->Mflag==0) fprintf(fs, "%d\n", nvar);
    for (nel_ta i=0; i<n; i++) {
      cls[i]=0;			//clear class data
      //we write to the input file:
      if (this->Mflag) {
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
    syscall=new char[strlen(this->options)+strlen(infile)+strlen(this->name)+strlen(outfile)+4];
    sprintf(syscall, "%s %s %s %s", this->options, infile, this->name, outfile);

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

    if (ncls!=this->ncls) {
      fprintf(stderr, "general_classifier: found %d classes in output file, expecting %d\n", ncls, this->ncls);
      exit(PARAMETER_OUT_OF_RANGE);
    }

    if (err!=n) {
      fprintf(stderr, "general_classifier: error reading output file, %s, exiting\n", outfile);
      exit(FILE_READ_ERROR);
    }

    delete [] syscall;
      
  }

  template <class real, class cls_t>
  dim_ta general_classifier<real, cls_t>::n_feat() {
    return -1;
  }

  template class direct_classifier<real_a, cls_ta>;
  template class agf_classifier<real_a, cls_ta>;
  template class knn_classifier<real_a, cls_ta>;
  template class general_classifier<real_a, cls_ta>;

} //end namespace libagf

