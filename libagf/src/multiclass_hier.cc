#include <stdio.h>
#include <getopt.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <assert.h>

#include <typeinfo>

#include "randomize.h"
#include "peteys_tmpl_lib.h"
#include "read_ascii_all.h"

#include "../../libpetey/linked.cc"

#include "agf_lib.h"

using namespace std;
using namespace libpetey;

namespace libagf {

  template <class real, class cls_t>
  multiclass_hier<real, cls_t>::multiclass_hier() {
    classifier=NULL;
    children=NULL;
    this->ncls=0;
  }

  //high level initialization for classification:
  template <class real, class cls_t>
  multiclass_hier<real, cls_t>::multiclass_hier(const char *file, int type, char *prefix, const char *com, int mf, int kf, int sigcode, int Zflag) {
    FILE *fs;
    int err;

    fs=fopen(file, "r");
    if (fs==NULL) {
      fprintf(stderr, "multiclass_hier: Unable to open control file, %s\n", file);
      throw UNABLE_TO_OPEN_FILE_FOR_READING;
    }
    err=init(fs, type, prefix, com, mf, kf, sigcode, Zflag);
    if (err!=0) throw err;
  }

  template <class real, class cls_t>
  multiclass_hier<real, cls_t>::multiclass_hier(FILE *fs, int type, char *prefix, const char *com, int mf, int kf, int sigcode, int Zflag) {
    int err=init(fs, type, prefix, com, mf, kf, sigcode, Zflag);
    if (err!=0) throw err;
  }

  //high level initializion for training purposes:
  template <class real, class cls_t>
  multiclass_hier<real, cls_t>::multiclass_hier(const char *file, int argc, char **argv, int maxstacksize) {
    FILE *fs;
    int err;
    fs=fopen(file, "r");
    if (fs==NULL) {
      fprintf(stderr, "multiclass_hier: Unable to open control file, %s\n", file);
      throw UNABLE_TO_OPEN_FILE_FOR_READING;
    }
    err=init(fs, argc, argv, maxstacksize);
    if (err!=0) throw err;
  }
    
  template <class real, class cls_t>
  multiclass_hier<real, cls_t>::multiclass_hier(FILE *fs, int argc, char **argv, int maxstacksize) {
    int err;
    err=init(fs, argc, argv, maxstacksize);
    if (err!=0) throw err;
  }
    
  template <class real, class cls_t>
  int multiclass_hier<real, cls_t>::init(FILE *fs, int type, char *prefix, const char *com, int mf, int kf, int sigcode, int Zflag) {
    multi_parse_param param;
    int err;

    //transfer method arguments to parameter structure:
    param.commandname=NULL;
    if (com!=NULL) {
      param.commandname=new char [strlen(com)+1];
      strcpy(param.commandname, com);
    }
    param.trainflag=0;		//we are not training a model here
    param.Mflag=mf;		//LIBSVM file format
    param.Kflag=kf;		//keep temporary files
    param.cw=1;			//constraint weight
    param.type=type;		//how to solve for the conditional prob.
    param.sigcode=sigcode;	//sigmoid function
    //need to set this with parameter later:
    param.prefix=prefix;	//path to data files
    param.Zflag=Zflag;		//use in house SVM codes

    param.infs=fs;		//input stream for control file
    param.lineno=0;		//keep track of line number

    err=init(param);
    delete [] param.commandname;

    return err;
  }

  //initialize for training purposes:
  template <class real, class cls_t>
  int multiclass_hier<real, cls_t>::init(FILE *fs, int argc, char **argv, int maxstacksize) {
    multi_parse_param param;
    int optlen;
    int err=0;

    //transfer method arguments to parameter structure:
    param.trainflag=1;
    param.prefix=NULL;

    param.infs=fs;
    param.lineno=0;

    //concatinate options:
    optlen=argc+1;
    for (int i=0; i<argc; i++) optlen+=strlen(argv[i]);
    param.maxnstack=maxstacksize;
    param.optstack=new char *[param.maxnstack];
    param.optstack[0]=new char [optlen];
    strcpy(param.optstack[0], "");
    for (int i=0; i<argc; i++) {
      strcat(param.optstack[0], argv[i]);
      strcat(param.optstack[0], " ");
    }
    for (int i=1; i<param.maxnstack; i++) param.optstack[i]=NULL;
    param.stackptr=1;

    //set non-training options:
    param.type=0;
    param.Mflag=0;

    init(param);

    //clean up:
    for (int i=0; i<param.maxnstack; i++) {
      if (param.optstack[i]!=NULL) delete [] param.optstack[i];
    }
    delete [] param.optstack;
    return err;
  }

  //low-level initialization (from a control structure):  
  template <class real, class cls_t>
  int multiclass_hier<real, cls_t>::init(multi_parse_param &param) {
    int err=0;
    char *fname0;			//name of file containing model
    char *fname;			//name of file containing model
    int flag;
    int c1;
    char c2;
    cls_t *cls;
    cls_t npart;
    cls_t topcls;			//possible class label at top of hierarchy
    char *options=NULL;		//options for direct classification

    fname0=parse_multi_start(&param, flag, options);
    if (param.prefix!=NULL) {
      //should really start using some better string libraries...
      fname=new char[strlen(param.prefix)+strlen(fname0)+1];
      sprintf(fname, "%s%s", param.prefix, fname0);
      delete fname0;
    } else {
      fname=fname0;
    }

    //three possibilities:
    assert(flag!=2);
    if (flag==1) {
      //multi-class classifier (non-hierarchical):
      classifier=new multiclass<real, cls_t>();
      ((multiclass<real, cls_t> *) classifier)->init(param);

      //we have to scan for the opening brackets:
      do {
        c1=fgetc(param.infs);
        if (c1==EOF) {
          fprintf(stderr, "multiclass_hier: %d, unexpected end of file(1).\n", param.lineno);
          throw FILE_READ_ERROR;
        }
        c2=(char) c1;
        //printf("6.(%c)\n", c2);
        if (c2=='{') {
          //we've moved to the next level:
          break;
        }
        if (isspace(c2)==0) {
          fprintf(stderr, "multiclass_hier: %d, syntax error in control file(2) at \"%c\".\n", param.lineno, c2);
          throw FILE_READ_ERROR;
        }
        if (c2=='\n') param.lineno++;
      } while (1);
    } else {
      //otherwise we read in the two-class model:
      //printf("multiclass_hier: initializing two class classifier, %s\n", fname);
      if (param.trainflag) {
        //if we are training a model, check for default options or
        //"previous option" option:
	//(*** options stack should be updated with modern containers 
	//and C++ string classes... ***)
        if (param.optstack[param.stackptr]!=NULL &&
                        strcmp(fname, ".")==0) {
          delete [] fname;
          fname=new char[strlen(param.optstack[param.stackptr])+1];
          strcpy(fname, param.optstack[param.stackptr]);
        } else {
          if (strlen(fname)==0) {
            delete [] fname;
            fname=new char[strlen(param.optstack[0])+1];
            strcpy(fname, param.optstack[0]);
          } else if (strcmp(fname, ".")==0) {
            delete [] fname;
            fname=new char[strlen(param.optstack[param.stackptr-1])+1];
            strcpy(fname, param.optstack[param.stackptr-1]);
          }
          if (param.optstack[param.stackptr]!=NULL) {
            delete [] param.optstack[param.stackptr];
          }
          param.optstack[param.stackptr]=new char[strlen(fname)+1];
          //option stack has to be independent of parsing operations:
          strcpy(param.optstack[param.stackptr], fname);
        }
        if (param.stackptr<param.maxnstack) {
          param.stackptr++;
        } else {
          fprintf(stderr, "multiclass_hier: option stack exausted (%d levels)\n", param.stackptr);
          throw PARAMETER_OUT_OF_RANGE;
        }
	if (flag=='A' || flag=='K' || flag=='G') {
	  //defer initializing until we know how many classes it has:
          classifier=NULL;
	} else {
          classifier=new binaryclassifier<real, cls_t>(fname);
	}
      } else {
        //performing classifications not training a model:
	//"direct" classifiers:
        if (flag=='A') {
          //AGF:
          classifier=new agf_classifier<real, cls_t>(fname, options);
	} else if (flag=='K') {
          //KNN:
          classifier=new knn_classifier<real, cls_t>(fname, options);
	} else if (flag=='G') {
          //external classifier (need to count number of classes first):
          classifier=NULL;
	} else if (param.Zflag) {
          //"in-house" SVM codes:
          classifier=new svm2class<real, cls_t>(fname);
	} else if (param.commandname==NULL) {
          //default borders classifier:
          //classifier=new borders_classifier<real, cls_t>(fname, param.sigcode);
          classifier=new borders_calibrated<real, cls_t>(fname);
        } else {
          //external binary classifier:
          classifier=new general2class<real, cls_t>(fname, 
		param.commandname, param.Mflag, param.Kflag);
        }
      }
    }
    if (options!=NULL) delete options;

    linked_list<classifier_obj<real, cls_t> *> brood;

    //go to the next level in the hierarchy:
    do {
      char *fname2;
      c1=scan_nowhitespace(param.infs, param.lineno);
      if (c1==EOF) {
        fprintf(stderr, "multiclass_hier: %d, unexpected end of file(2). (missing closing brace?)\n", param.lineno);
        throw FILE_READ_ERROR;
      }
      c2=(char) c1;
      if (c2=='}') break;
      fseek(param.infs, -1, SEEK_CUR);
      fname2=scan_class_label(param.infs, param.lineno);
      if (fname2==NULL) {
        multiclass_hier<real, cls_t> *ch;
        ch=new multiclass_hier<real, cls_t>();
        ch->init(param);
        brood.add(ch);
      } else {
        oneclass<real, cls_t> *ch;
        sscanf(fname2, "%d", &topcls);
        //printf("_hier: oneclass with class, %d\n", topcls);
        //"one-class" classifier:
        ch=new oneclass<real, cls_t>(topcls);
        brood.add(ch);
        delete [] fname2;
      }
    } while(1);

    long cnt;
    children=brood.make_array(cnt);
    npart=cnt;

    //printf("multiclass_hier: found %d children\n", npart);

    //now that we've counted the number of children:
    if (classifier==NULL) {
      if (param.trainflag) {
        classifier=new direct_classifier<real, cls_t>(npart, fname, flag, param.Mflag);
      } else {
        classifier=new general_classifier<real, cls_t>(fname, options, npart,
			param.Mflag, param.Kflag);
      }
    } else if (npart!=classifier->n_class()) {
      fprintf(stderr, "multiclass_hier: expected %d children, found %d\n", classifier->n_class(), npart);
      throw SAMPLE_COUNT_MISMATCH;
    }

    delete [] fname;

    //check to make sure everything's Kosher:
    this->ncls=-1;
    n_class();
    //printf("_hier->ncls=%d\n", this->ncls);
    cls=new cls_t[this->ncls];
    class_list(cls);
    heapsort_inplace(cls, this->ncls);

    //printf("multiclass_hier class list: %d\n", cls[0]);
    for (int i=1; i<this->ncls; i++) {
      //printf("multiclass_hier class list: %d\n", cls[i]);
      if (cls[i-1]==cls[i]) {
        fprintf(stderr, "multiclass_hier: error, duplicate classes\n");
        throw OTHER_ERROR;
      }
    }

    if (max_depth()==1) nonh_flag=1; else nonh_flag=0;

    //delete everything in the stack above the current pointer:
    if (param.trainflag) {
      if (param.optstack[param.stackptr+1]!=NULL) {
        delete [] param.optstack[param.stackptr+1];
        param.optstack[param.stackptr+1]=NULL;
      }
      param.stackptr--;
    }

    delete [] cls;

    return err;
  }

  template <class real, class cls_t>
  multiclass_hier<real, cls_t>::~multiclass_hier() {
    if (children!=NULL) {
      for (int i=0; i<classifier->n_class(); i++) delete children[i];
      delete [] children;
    }
    if (classifier!=NULL) delete classifier;
  }

  template <class real, class cls_t>
  cls_t multiclass_hier<real, cls_t>::classify(real *x, real &p, real *praw) {
    cls_t cls1, cls2;
    real pdum;

    cls1=classifier->classify(x, pdum, praw);
    cls2=children[cls1]->classify(x, p, praw);
    //printf("cls1=%d; cls2=%d; pdum=%g; p=%g; psum=%g\n", cls1, cls2, pdum, p, p*pdum);
    p=p*pdum;

    return cls2;
  }

  template <class real, class cls_t>
  cls_t multiclass_hier<real, cls_t>::classify(real *x, real *p, real *praw) {
    cls_t cls1, cls2;

    if (nonh_flag) {
      real pdum;
      cls1=classifier->classify(x, p, praw);
      cls2=children[cls1]->classify(x, pdum, praw);
    } else {
      cls2=this->classifier_obj<real, cls_t>::classify(x, p, praw);
    }

    return cls2;
  }

  template <class real, class cls_t>
  int multiclass_hier<real, cls_t>::ltran_model(real **mat, real *b, dim_ta d1, dim_ta d2) {
    int err=0;

    err=classifier->ltran_model(mat, b, d1, d2);
    if (err!=0) {
      fprintf(stderr, "multiclass_hier::ltran_model: an error occurred transforming partition classifier\n");
      return err;
    }
    for (int i=0; i<classifier->n_class(); i++) {
      err=children[i]->ltran_model(mat, b, d1, d2);
      if (err!=0) {
        fprintf(stderr, "multiclass_hier::ltran_model: an error occurred transforming child #%d\n", i);
        return err;
      }
    }
    if (this->mat==NULL) {
      this->D1=0;
      this->D=n_feat();
    } else {
      assert(n_feat()==d2);
    }
    return err;
  }

  template <class real, class cls_t>
  cls_t multiclass_hier<real, cls_t>::n_class() {
    cls_t nchild;

    if (this->ncls<=0) {
      this->ncls=0;
      nchild=classifier->n_class();
      for (cls_t i=0; i<nchild; i++) this->ncls+=children[i]->n_class();
    }

    return this->ncls;
  }

  template <class real, class cls_t>
  dim_ta multiclass_hier<real, cls_t>::n_feat() {
    cls_t nchild;
    cls_t D2;

    if (this->D1<=0) {
      //fprintf(stderr, "multiclass_hier: checking dimensions of children\n");
      nchild=classifier->n_class();
      this->D1=classifier->n_feat();
      for (cls_t i=0; i<nchild; i++) {
        if (children[i]->max_depth() == 0) continue;
        D2=children[i]->n_feat();
        if (D2!=this->D1) {
          fprintf(stderr, "multiclass_hier: number of features in classifier does not mathc that in child %d", i);
          fprintf(stderr, "                 %d vs. %d\n", this->D1, D2);
          throw DIMENSION_MISMATCH;
        }
      }
    }

    return this->D1;
  }

  template <class real, class cls_t>
  int multiclass_hier<real, cls_t>::max_depth(int cur) {
    int maxdepth;
    int depth;
    int npart;

    npart=classifier->n_class();
    maxdepth=0;
    cur++;
    for (int i=0; i<npart; i++) {
      depth=children[i]->max_depth(cur);
      if (depth>maxdepth) maxdepth=depth;
    }

    return maxdepth+1;

  }

  template <class real, class cls_t>
  cls_t multiclass_hier<real, cls_t>::class_list(cls_t *cls) {
    cls_t nchild;
    cls_t nc_child;
    cls_t ncls1=0;

    nchild=classifier->n_class();
    //printf("_hier::class_list: %d children\n", nchild);
    for (cls_t i=0; i<nchild; i++) {
      nc_child=children[i]->class_list(cls+ncls1);
      //printf("_hier::class_list: child %d has %d classes\n", i, nc_child);
      ncls1+=nc_child;
    }
    assert(this->ncls==ncls1);

    return this->ncls;
  }

  template <class real, class cls_t>
  void multiclass_hier<real, cls_t>::batch_classify(real **x, cls_t *cls, real *p, nel_ta n, dim_ta nvar) {
    cls_t cls1[n];
    cls_t ncls1;		//number of children
    cls_t *nlab;		//number of instances of each class
    nel_ta *ind;		//one index for each child
    cls_t **cls2;		//class results from each child
    real **p2;			//probabilities from each child
    real ***x2;			//feature data to pass to each child

    if (n==0) return;

    //this would be really easy in IDL...
    //printf("multiclass_hier: batch_classify for partition\n");
    classifier->batch_classify(x, cls1, p, n, nvar);
    ncls1=classifier->n_class();
    nlab=new cls_t[ncls1];
    for (cls_t i=0; i<ncls1; i++) nlab[i]=0;
    for (nel_ta i=0; i<n; i++) {
      if (cls1[i]>=0 && cls1[i]<=ncls1-1) {
        nlab[cls1[i]]++;
      } else {
        fprintf(stderr, "multiclass_hier->batch_classify: warning, class value (%d) out of range [0, %d), skipping\n", cls1[i], ncls1);
      }
    }

    //allocate results arrays:    
    cls2=new cls_t *[ncls1];		//this one for class data
    cls2[0]=new cls_t[n];

    p2=new real *[ncls1];		//this one for probabilities
    p2[0]=new real[n];

    //we need a place to put the different feature sets:
    //printf("Found:\n%d %d labels\n", nlab[0], 0);
    x2=new real **[ncls1];
    x2[0]=new real*[n];
    for (nel_ta i=1; i<ncls1; i++) {
      //printf("%d %d labels\n", nlab[i], i);
      cls2[i]=cls2[i-1]+nlab[i-1];
      p2[i]=p2[i-1]+nlab[i-1];
      x2[i]=x2[i-1]+nlab[i-1];
    }

    //arrange the pointers:
    ind=new nel_ta[ncls1];
    for (dim_ta i=0; i<ncls1; i++) ind[i]=0;
    for (nel_ta i=0; i<n; i++) {
      cls_t cl=cls1[i];
      if (cl<0 || cl>=ncls1) continue;
      x2[cl][ind[cl]]=x[i];
      ind[cl]++;
    }

    //perform classifications with each of the children:
    for (cls_t i=0; i<ncls1; i++) {
      //printf("multiclass_hier: batch_classify for child %d\n", i);
      children[i]->batch_classify(x2[i], cls2[i], p2[i], nlab[i], nvar);
    }

    //resolve the results (classes in the right place, probabilities multiplied through):
    for (dim_ta i=0; i<ncls1; i++) ind[i]=0;
    for (nel_ta i=0; i<n; i++) {
      cls_t cl=cls1[i];
      if (cl<0 || cl>=ncls1) continue;
      cls[i]=cls2[cl][ind[cl]];
      p[i]=p[i]*p2[cl][ind[cl]];
      ind[cl]++;
    }

    //clean up:
    delete [] ind;
    delete [] nlab;
  
    delete [] x2[0];
    delete [] x2;
    delete [] cls2[0];
    delete [] cls2;
    delete [] p2[0];
    delete [] p2;
  }

  template <class real, class cls_t>
  void multiclass_hier<real, cls_t>::batch_classify(real **x, cls_t *cls, real **p, nel_ta n, dim_ta nvar) {
    if (nonh_flag) {
      cls_t list[this->ncls];
      class_list(list);
      classifier->batch_classify(x, cls, p, n, nvar);
      for (nel_ta i=0; i<n; i++) cls[i]=list[cls[i]];
    } else {
      this->classifier_obj<real, cls_t>::batch_classify(x, cls, p, n, nvar);
    }
  }

  template <class real, class cls_t>
  void multiclass_hier<real, cls_t>::print(FILE *fs, char *fbase, int depth) {
    cls_t npart;
    char *fbase2=NULL;
    if (fbase!=NULL) fbase2=new char[strlen(fbase)+4];
    
    classifier->print(fs, fbase, depth);
    npart=classifier->n_class();
    fprintf(fs, "\n");
    for (int i=0; i<depth; i++) fprintf(fs, "  ");
    fprintf(fs, "{\n");
    for (cls_t i=0; i<npart; i++) {
      if (fbase!=NULL) sprintf(fbase2, "%s.%2.2d", fbase, i);
      children[i]->print(fs, fbase2, depth+1);
    }
    for (int i=0; i<depth; i++) fprintf(fs, "  ");
    fprintf(fs, "}\n");

    if (fbase2!=NULL) delete [] fbase2;
  }

  template <class real, class cls_t>
  int multiclass_hier<real, cls_t>::generate_commands(FILE *fs,
		char *train, char *fbase, char *command, char *partcom,
		char *concom, int Kflag, unsigned long session_id) {
    multi_train_param param;
    cls_t *clist;

    param.commandfs=fs;
    param.Kflag=Kflag;
    param.train=new char[strlen(train)+1];
    strcpy(param.train, train);
    if (command!=NULL) {
      param.commandname=new char[strlen(command)+1];
      strcpy(param.commandname, command);
    } else {
      param.commandname=new char[strlen(AGF_COMMAND_PREFIX)+
                strlen(AGF_BINARY_CLASSIFIER)+strlen(AGF_OPT_VER)+20];
      sprintf(param.commandname, "%s%s%s", AGF_COMMAND_PREFIX,
                AGF_BINARY_CLASSIFIER, AGF_OPT_VER);
    }

    //set optional arguments:
    param.partcom=NULL;
    param.concom=NULL;
    param.precom=NULL;		//not used
    if (partcom!=NULL) {
      param.partcom=new char [strlen(partcom)+1];
      strcpy(param.partcom, partcom);
      if (concom!=NULL) {
        param.concom=new char [strlen(concom)+1];
        strcpy(param.concom, concom);
      }
      if (session_id==0) param.session_id=seed_from_clock();
		else param.session_id=session_id;
    }

    clist=new cls_t[this->ncls+1];

    commands(param, &clist, fbase);

    delete [] param.train;
    delete [] param.commandname;
    delete [] param.partcom;
    delete [] param.concom;
    delete [] clist;
  }

  template <class real, class cls_t>
  int multiclass_hier<real, cls_t>::commands(multi_train_param &param, cls_t **clist, char *fbase) {
    cls_t **clist2;		//list of class labels for each child
    cls_t nchild=classifier->n_class();
    char *fbase2;

    clist2=new cls_t*[nchild+1];
    clist2[0]=*clist;
    fbase2=new char[strlen(fbase)+4];

    for (cls_t i=0; i<nchild; i++) {
      sprintf(fbase2, "%s.%2.2d", fbase, i);
      children[i]->commands(param, clist2+i, fbase2);
      clist2[i+1]=clist2[i]+children[i]->n_class();
    }

    classifier->commands(param, clist2, fbase);

    delete [] clist2;
    delete [] fbase2;

    return this->ncls;
  }

  //set id's of binary classifiers for collating raw probabilities:
  template <class real, class cls_t>
  void multiclass_hier<real, cls_t>::set_id(cls_t *id) {
    cls_ta nchild=classifier->n_class();
    if (nchild>2) {
      cls_t *idlist=new cls_t[nchild-1];
      for (cls_t i=0; i<nchild-1; i++) {
        children[i]->set_id(id);
        idlist[i]=*id;
        (*id)++;
      }
      children[nchild-1]->set_id(id);
      classifier->set_id(idlist);
      delete [] idlist;
    } else {
      //just like traversing a binary tree:
      children[0]->set_id(id);
      classifier->set_id(id);
      children[1]->set_id(id);
    }
  }

  template <typename real, typename cls_t>
  int multiclass_hier<real, cls_t>::load(FILE *fs, int ct) {
    char **sub;
    int nsub;
    char *type;
    type=fget_line(fs, 1);	//first line describes type
    if (strcmp(type, "1v1")!=0 && strcmp(type, "1vR")!=0 && strcmp(type, "ADJ")!=0) {
      delete [] type;
      return PARAMETER_OUT_OF_RANGE;
    }
    char *line=fget_line(fs);		//second line has number of classes
    sscanf(line, "%d", &this->ncls);
    delete [] line;
    //third line contains class labels:
    line=fget_line(fs, 1);
    sub=split_string_destructive(line, nsub);
    if (nsub!=this->ncls) throw DIMENSION_MISMATCH;
    children=new classifier_obj<real, cls_t> *[this->ncls];
    for (int i=0; i<this->ncls; i++) {
      children[i]=new oneclass<real, cls_t>(atoi(sub[i]));
    }
    delete [] line;
    delete [] sub;
    classifier=new multiclass<real, cls_t>(ct);
    rewind(fs);
    classifier->load(fs);
    nonh_flag=1;
    delete [] type;
    return 0;
  }

  template <typename real, typename cls_t>
  int multiclass_hier<real, cls_t>::load(FILE *fs) {
    return load(fs, -1);
  }

  template <typename real, typename cls_t>
  int multiclass_hier<real, cls_t>::save(FILE *fs) {
    if (nonh_flag==0) {
      fprintf(stderr, "multiclass_hier::save: cannot save; not the right type\n");
      throw PARAMETER_OUT_OF_RANGE;
    }
    classifier->save(fs);
    fseek(fs, 9, SEEK_SET);		//***ugly hack
    					//(could stick labels at very end...)
    for (int i=0; i<this->ncls; i++) children[i]->save(fs);
    return 0;
  }

  template <typename real, typename cls_t>
  void multiclass_hier<real, cls_t>::train(real **train, cls_t *cls, nel_ta ntrain, int type, real *param) {
    cls_t *map;				//for partitioning the classes
    cls_t maxcls=0;			//largest value for class label
    cls_t label[this->ncls];
    cls_t nchild=classifier->n_class();		//number of children
    cls_t cls2[ntrain];

    //partition the classes (this part is quite general--should make a
    //generalized routine that calls a different routine at the end...):
    this->class_list(label);		//have to add this extra bullshit
    for (cls_t i=0; i<this->ncls; i++) if (label[i]>=maxcls) maxcls=label[i]+1;
    map=new cls_t[maxcls];
    for (cls_t i=0; i<maxcls; i++) map[i]=-1;

    cls_t k=0;
    for (cls_t i=0; i<nchild; i++) {
      cls_t nchildcls;
      nchildcls=children[i]->class_list(label);
      for (cls_t j=0; j<nchildcls; j++) map[j+k]=i;
      k+=nchildcls;
    }
    for (nel_ta i=0; i<ntrain; i++) {
      if (cls[i]<0 || cls[i]>maxcls) {
        cls2[i]=-1;
      } else {
        cls2[i]=map[cls[i]];
      }
    }
    classifier->train(train, cls2, ntrain, type, param);
    for (cls_t i=0; i<nchild; i++) {
      children[i]->train(train, cls, ntrain, type, param);
    }
    delete [] map;
  }

  template class multiclass_hier<real_a, cls_ta>;

}

