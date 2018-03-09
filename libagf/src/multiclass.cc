#include <stdio.h>
#include <getopt.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <assert.h>

#include <vector>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>

#include "constrained.h"
#include "full_util.h"
#include "gsl_util.h"
#include "read_ascii_all.h"
#include "peteys_tmpl_lib.h"

#include "agf_lib.h"

using namespace std;
using namespace libpetey;

namespace libagf {

  template <typename real, typename cls_t>
  multiclass<real, cls_t>::multiclass(int ct) {
    type=ct;
    set_solve_type(type);

    //set model variables to NULL:
    twoclass=NULL;

    //set number parameters to 0:
    nmodel=0;
    this->ncls=0;

    //"polarity":
    pol=NULL;

  }

  template <typename real, typename cls_t>
  multiclass<real, cls_t>::multiclass(const char *file, int clstyp, const char *com, int mf, int kf, int sigcode, int Zflag) {
    int err;
    multi_parse_param param;

    param.infs=fopen(file, "r");
    if (param.infs==NULL) {
      fprintf(stderr, "multiclass: Unable to open control file, %s\n", file);
      exit(UNABLE_TO_OPEN_FILE_FOR_READING);
    }
    param.trainflag=0;
    param.lineno=0;

    fclose(param.infs);
    //multi_partition_strict(fname, part, nmodel);

    param.commandname=NULL;
    if (com!=NULL) {
      param.commandname=new char [strlen(com)+1];
      strcpy(param.commandname, com);
    }
    param.Mflag=mf;
    param.Kflag=kf;
    param.Zflag=Zflag;
    param.sigcode=sigcode;  
    err=init(param);
    if (err!=0) exit(err);

    //clean up:
    delete [] param.commandname;

    //classification method:
    type=clstyp;
    set_solve_type(type);

    //"polarity":
    pol=NULL;
  }

  template <typename real, typename cls_t>
  int multiclass<real, cls_t>::init(multi_parse_param &param) {
    //char *name[MAXNPART];
    int npart;
    char **name;
    //cls_t *part[MAXNPART*2];
    cls_t **coding_matrix;
    int ncls;
    int c1;
    char c2;
    int err=0;

    //nmodel=parse_multi_partitions(&param, name, part, MAXNPART);
    parse_multi_partitions(&param, name, coding_matrix, npart, ncls);
    fseek(param.infs, -1, SEEK_CUR);

    //pass to another initialization routine (duh...):
    init(name, coding_matrix, npart, ncls, param.prefix, param.type, 
		    param.trainflag, param.commandname, 
		    param.Mflag, param.Kflag, param.Zflag);

    //clean up:
    for (int i=0; i<nmodel; i++) {
      delete [] name[i];
    }
    delete [] name;
    delete [] coding_matrix[0];
    delete [] coding_matrix;

    return err;
  }

  template <typename real, typename cls_t>
  int multiclass<real, cls_t>::init(char **fname, int **coding_matrix, int npart, cls_t ncls, char *prefix, int method,
		int tflag, char *com, int Mflag, int Kflag, int Zflag) {
    nmodel=npart;

    //classification type and constraint weight:
    type=method;
    set_solve_type(type);

    //pass in the number of classes:
    this->ncls=ncls;

    //add file name prefix if applicable:
    if (prefix!=NULL) {
      char *fname1;
      for (int i=0; i<nmodel; i++) {
        fname1=new char[strlen(prefix)+strlen(fname[i])+1];
	sprintf(fname1, "%s%s", prefix, fname[i]);
	delete fname[i];
	fname[i]=fname1;
      }
    }

    //initialize binary classifiers:
    twoclass=new binaryclassifier<real, cls_t> *[nmodel];
    for (cls_t i=0; i<nmodel; i++) {
      if (tflag) {
        twoclass[i]=new binaryclassifier<real, cls_t>(fname[i]);
      } else {
	if (com==NULL) {
          //twoclass[i]=new borders_classifier<real, cls_t>(fname[i], sigcode);
          //if (binparam==NULL) 
	  //just bite the bullet and use a couple of global variables??
	  // -- don't think the algorithm's thread safe anyways ...
          twoclass[i]=binclass_init<real, cls_t>(fname[i], Zflag);
	  //		else twoclass[i]=new binclass(fname[i], binparam);
        } else {
          twoclass[i]=new general2class<real, cls_t>(fname[i], 
			com, Mflag, Kflag);
        }
      }
    }

    //create the mapping:
    code=zero_matrix<real>(nmodel, this->ncls);	//only one format

    for (int i=0; i<nmodel; i++) {
      for (cls_t j=0; j<this->ncls; j++) {
        code[i][j]=coding_matrix[i][j];
      }
    }

    strictflag=1;
    for (int i=0; i<npart; i++) {
      for (cls_t j=0; j<ncls; j++) {
        if (code[i][j]==0) {
          strictflag=0;
	  break;
	}
      }
      if (strictflag==0) break;
    }

    return 0;
  }

  template <typename real, typename cls_t>
  multiclass<real, cls_t>::~multiclass() {
    //delete the binary classifiers and the decision matrix:
    for (int i=0; i<nmodel; i++) delete twoclass[i];
    delete [] twoclass;

    if (pol!=NULL) delete [] pol;

    delete_matrix(code);
  }

  template <typename real, typename cls_t>
  void multiclass<real, cls_t>::set_solve_type(int ct) {
    switch (ct) {
      case (0):
        solve_class=&solve_class_constrained2<real, real>;
	break;
      case (1):
	solve_class=&solve_class_scratch<real, real>;
        break;
      case (2):
	solve_class=&solve_class_vote_pdf<real, real>;
        break;
      case (3):
	solve_class=&solve_class_vote<real, real>;
        break;
      case (4):
	solve_class=&solve_class_norm1<real, real>;
        break;
      case (5):
	solve_class=&solve_class_norm2<real, real>;
        break;
      case (6):
	//solve_class=&solve_class_renorm<real, real>;
	solve_class=&solve_class_renorm<real>;
        break;
      case (7):
        solve_class=&solve_class_constrained1<real, real>;
	break;
      case (8):
	solve_class=&solve_class_vote_pdf2<real, real>;
	break;
      case (9):
	solve_class=&solve_class_1vR<real, real>;
	break;
      case (10):
	solve_class=&solve_class_Zadrozny<real, real>;
	break;
      case (11):
	solve_class=&solve_class_interior<real, real>;
	break;
      default:
        solve_class=&solve_class_constrained2<real, real>;
        break;
    }
  }

  //load a linear transformation to apply to the test points:
  template <typename real, typename cls_t>
  int multiclass<real, cls_t>::ltran_model(real **mat, real *b, dim_ta d1, dim_ta d2) {
    int err=0;
    //just pass it to the binary classifiers:
    for (int i=0; i<nmodel; i++) {
      err=twoclass[i]->ltran_model(mat, b, d1, d2);
      if (err!=0) {
        fprintf(stderr, "multiclass::ltran: an error occured transforming partition #%d\n", i);
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

  template <typename real, typename cls_t>
  cls_t multiclass<real, cls_t>::classify(real *x, real *p, real *praw) {
    real r[nmodel];
    real pt=0;

    for (int i=0; i<nmodel; i++) {
      r[i]=twoclass[i]->R(x, praw);
      //printf("%12.6g", r[i]);
    }
    //printf("\n");

    (*solve_class)(code, nmodel, this->ncls, r, p);
    //for (cls_t i=0; i<this->ncls; i++) pt+=p[i];
    //printf("pt=%g\n", pt);
    return choose_class(p, this->ncls);
  }

  template <typename real, typename cls_t>
  dim_ta multiclass<real, cls_t>::n_feat() {
    cls_t D2;

    if (this->D1<=0) {
      this->D1=twoclass[0]->n_feat();
      for (cls_t i=1; i<nmodel; i++) {
        D2=twoclass[i]->n_feat();
        if (D2!=this->D1) {
          fprintf(stderr, "multiclass: number of features in classifier %d does not match that in child %d", 0, i);
          fprintf(stderr, "                 %d vs. %d\n", this->D1, D2);
          exit(DIMENSION_MISMATCH);
        }
      }
      if (this->mat==NULL) this->D=this->D1;
    }

    return this->D1;
  }

  template <typename real, typename cls_t>
  void multiclass<real, cls_t>::batch_classify(real **x, cls_t *cls, real **p1, nel_ta n, dim_ta nvar) {
    real **r1;
    real r2[nmodel];

    //printf("multiclass: performing classifications with %d test vectors on %d models\n", n, nmodel);
    r1=allocate_matrix<real>(nmodel, n);
    for (cls_t i=0; i<nmodel; i++) {
      //printf("multiclass: batchR model %d\n", i);
      twoclass[i]->batchR(x, r1[i], n, nvar);
    }
    for (nel_ta i=0; i<n; i++) {
      for (int j=0; j<nmodel; j++) r2[j]=r1[j][i];
      (*solve_class)(code, nmodel, this->ncls, r2, p1[i]);
      cls[i]=choose_class(p1[i], this->ncls);
    }

    delete_matrix(r1);

  }

  //for accessing "raw probabilities" for use in continuum predictions:
  template <typename real, typename cls_t>
  void multiclass<real, cls_t>::set_id(cls_t *id) {
    cls_t niter1;
    //this does actually makes sense: check multiclass_hier.cc for the
    //explanation...
    if (nmodel<this->ncls-1) niter1=nmodel; else niter1=this->ncls-1;
    for (cls_t i=0; i<niter1; i++) {
      twoclass[i]->set_id(id+i);
      //remove side-effect:
      id[i]--;
    }
    for (cls_t i=niter1; i<nmodel; i++) {
      twoclass[i]->set_id(id+this->ncls-1);
    }

    //for (cls_t i=0; i<nmodel; i++) twoclass[i]->set_id(id+i);
    //return nmodel;
  }

  //since they do about the same thing, we put them next to each other...
  template <typename real, typename cls_t>
  cls_t multiclass<real, cls_t>::collect_binary_classifiers(binaryclassifier<real, cls_t> **list) {
    for (cls_t i=0; i<nmodel; i++) {
      //this makes no sense:
      twoclass[i]->collect_binary_classifiers(list+i);
      //list[i]=twoclass[i];
    }
    return nmodel;
  }

  template <typename real, typename cls_t>
  void multiclass<real, cls_t>::print(FILE *fs, char *fbase, int depth) {
    char *fbase2=NULL;
    if (fbase!=NULL) fbase2=new char[strlen(fbase)+5];

    for (int i=0; i<nmodel; i++) {
      if (fbase!=NULL) sprintf(fbase2, "%s-%2.2d", fbase, i);
      twoclass[i]->print(fs, fbase2, depth);
      for (int j=0; j<this->ncls; j++) {
        if (code[i][j]<0) fprintf(fs, " %d", j);
      }
      fprintf(fs, " /");
      for (int j=0; j<this->ncls; j++) {
        if (code[i][j]>0) fprintf(fs, " %d", j);
      }
      fprintf(fs, ";");
      if (i!=nmodel-1) fprintf(fs, "\n");
    }

    if (fbase2!=NULL) delete [] fbase2;
  }

  template <typename real, typename cls_t>
  int multiclass<real, cls_t>::commands(multi_train_param &param,
                cls_t **clist, char *fbase) {
    char *fbase2=NULL;
    cls_t *clist2;
    cls_t *clist3[3];
    cls_t nct=0;		//number of classes total
    cls_t nc1, nc2;		//number in first, second partition resp.

    //how many classes total?
    for (cls_t i=0; clist[0]+i != clist[this->ncls]; i++) nct++;
    clist2=new cls_t[nct+1];		//need that extra element hanging off
					//the end as a flag for stopping iteration

    if (fbase!=NULL) fbase2=new char[strlen(fbase)+4];

    for (int i=0; i<nmodel; i++) {
      //generate file name:
      if (fbase!=NULL) sprintf(fbase2, "%s-%2.2d", fbase, i);
      //print command name, options, training data, output model:
      //fprintf(param.commandfs, "%s %s %s %s", param.commandname, name[i],
	//	param.train, fbase2);

      //print class partitions:
      //gather the class labels in each partition:
      nc1=0;
      clist3[0]=clist2;
      for (cls_t j=0; j<this->ncls; j++) {
        if (code[i][j]<0) {
          for (cls_t k=0; clist[j]+k!=clist[j+1]; k++) {
            clist2[nc1]=clist[j][k];
            nc1++;
            //fprintf(param.commandfs, " %d", clist[j][k]);
          }
        }
      }
      clist3[1]=clist3[0]+nc1;
      //fprintf(param.commandfs, " /");
      nc2=0;
      for (cls_t j=0; j<this->ncls; j++) {
        if (code[i][j]>0) {
          for (cls_t k=0; clist[j]+k!=clist[j+1]; k++) {
            clist2[nc1+nc2]=clist[j][k];
            nc2++;
            //fprintf(param.commandfs, " %d", clist[j][k]);
          }
        }
      }
      clist3[2]=clist3[1]+nc2;

      //pass the whole business one level up:
      twoclass[i]->commands(param, clist3, fbase2);
    }
    //clean up:
    if (fbase2!=NULL) delete [] fbase2;
    delete [] clist2;

    return this->ncls;
  }

  template <typename real, typename cls_t>
  int multiclass<real, cls_t>::get_code(int **code2, char **model) {
    //int dum[2];		//no idea why this doesn't work...
    int *dum=new int[2];
    for (int i=0; i<nmodel; i++) {
      for (int j=0; j<this->ncls; j++) code2[i][j]=code[i][j];
      twoclass[i]->get_code(&dum, model+i);
    }
    delete [] dum;
    return nmodel;
  }
		    
  template <typename real, typename cls_t>
  int multiclass<real, cls_t>::detect_type() {
    int spec_type;		//0=1 v. rest; 1=1 v. 1; 2=adj.
    int **proto;			//prototype to compare against
    int ind[nmodel];
    long dum[nmodel];
    vector<int> sd[nmodel];
    vector<int> p2[nmodel];

    for (int i=0; i<nmodel; i++) {
      sd[i].resize(this->ncls);
      for (int j=0; j<this->ncls; j++) sd[i][j]=code[i][j];
    }

    spec_type=-1;			//which type is it?
    pol=new int[nmodel];
    for (int tc=0; tc<3; tc++) {
      int nm;			//required number of binary classifiers
      switch (tc) {
        case(0):
          //check for 1 v rest:
          proto=one_against_all<int>(this->ncls);
	  nm=this->ncls;
	  break;
        case(1):
          //check for 1 v 1:
          proto=one_against_one<int>(this->ncls);
	  nm=this->ncls*(this->ncls-1)/2;
	  break;
        case(2):
	  //check for adjacent partitioning:
          proto=partition_adjacent<int>(this->ncls);
	  nm=this->ncls-1;
	  break;
      }
      if (nmodel==nm) {
        for (int i=0; i<nm; i++) {
          p2[i].resize(this->ncls);
          for (int j=0; j<this->ncls; j++) p2[i][j]=proto[i][j];
        }
	//fuck it, just use a brute-force method:
	for (int i=0; i<nm; i++) ind[i]=-1;
        for (int i=0; i<nm; i++) {
          for (int j=0; j<nm; j++) {
            if (sd[i]!=p2[j]) {
              vector<int> tmp(this->ncls);
	      //should have copy of negative coding matrix for "efficiency"
	      //ha ha...
              for (int k=0; k<this->ncls; k++) tmp[k]=-p2[j][k];
              if (sd[i]==tmp) {
                pol[i]=-1;
                ind[i]=j;
	      }
            } else {
              pol[i]=1;
              ind[i]=j;
	    }
	  }
	}
	spec_type=tc;
	for (int i=0; i<nm; i++) {
          printf("%d ", ind[i]);
          if (ind[i] == -1) {
            spec_type=-1;
	    break;
	  }
	}
	printf("\n");
      }
      if (spec_type!=-1) break;
      delete [] proto[0];
      delete [] proto;
    }

    if (spec_type==-1) {
      delete [] pol;
      pol=NULL;
    } else {
      //rearrange the coding matrix and binary classifiers so they're in
      //the "right" order:
      real ** newcode=allocate_matrix<real, int>(nmodel, this->ncls);
      binaryclassifier<real, cls_t> **twoclass2=new binaryclassifier<real, cls_t>*[nmodel];
      for (int i=0; i<nmodel; i++) {
        for (int j=0; j<this->ncls; j++) {
          newcode[i][j]=code[ind[i]][j];
        }
	twoclass2[i]=twoclass[ind[i]];
      }
      delete [] code[0];
      delete [] code;
      code=newcode;

      delete [] twoclass;
      twoclass=twoclass2;
      type=10+spec_type;
    }

    return spec_type;
  }

  template <typename real, typename cls_t>
  int multiclass<real, cls_t>::load(FILE *fs) {
    int err=0;
    char **sub;
    int nsub;
    char *typestr=fget_line(fs, 1);
    char *line=fget_line(fs);
    sscanf(line, "%d", &this->ncls);
    delete [] line;
    //type specific stuff:
    if (strcmp(typestr, "1vR")==0) {
      nmodel=this->ncls;
      code=one_against_all<real>(this->ncls);
      if (type<0) {
        type=10;		//deprecated, really...
        solve_class=&solve_class_1vR<real>;
      }
      strictflag=1;
    } else if (strcmp(typestr, "1v1")==0) {
      nmodel=this->ncls*(this->ncls-1)/2;
      code=one_against_one<real>(this->ncls);
      if (type<0) {
        type=11;
        solve_class=&solve_class_norm2<real>;
      }
      strictflag=0;
    } else if (strcmp(typestr, "ADJ")==0) {
      nmodel=this->ncls-1;
      code=partition_adjacent<real>(this->ncls);
      if (type<0) {
        type=12;
	solve_class=&solve_class_constrained2<real>;
      }
      strictflag=1;
    } else {
      fprintf(stderr, "multiclass::load: type, %s, not recognized\n", typestr);
      throw PARAMETER_OUT_OF_RANGE;
    }
    line=fget_line(fs);
    //don't need labels:
    delete [] line;
    //"polarity":
    line=fget_line(fs);
    sub=split_string_destructive(line, nsub);
    if (nsub<nmodel) {
      fprintf(stderr, "multiclass::load: Not enough \"polarities\" found in initialization file\n");
      fprintf(stderr, "  %d vs. %d\n", this->ncls, nsub);
      throw SAMPLE_COUNT_MISMATCH;
    }
    pol=new int[nmodel];
    for (int i=0; i<nmodel; i++) {
      pol[i]=atoi(sub[i]);
      for (int j=0; j<this->ncls; j++) code[i][j]=pol[i]*code[i][j];
    }
    delete [] line;
    delete [] sub;

    //read in binary classifiers:
    twoclass=new binaryclassifier<real, cls_t>*[nmodel];
    for (int i=0; i<nmodel; i++) {
      twoclass[i]=new borders_classifier<real, cls_t>();
      err=twoclass[i]->load(fs);
      if (err!=0) {
        fprintf(stderr, "multiclass::load: error loading border vectors\n");
        throw err;
      }
    }

    return err;
  }

  template <typename real, typename cls_t>
  int multiclass<real, cls_t>::save(FILE *fs) {
    int err=0;
    detect_type();
    switch (type) {
      case(10):
        fprintf(fs, "1vR\n");
	break;
      case(11):
	fprintf(fs, "1v1\n");
	break;
      case(12):
	fprintf(fs, "ADJ\n");
	break;
      default:
	fprintf(stderr, "multiclass::save: not a recognized type\n");
	throw PARAMETER_OUT_OF_RANGE;
    }
    fprintf(fs, "%4d\n", this->ncls);
    //need place filler for labels (how are we going to do this?):
    for (int i=0; i<this->ncls; i++) fprintf(fs, "%4d", i);
    fprintf(fs, "\n");
    for (int i=0; i<nmodel; i++) fprintf(fs, "%d ", pol[i]);
    fprintf(fs, "\n");
    for (int i=0; i<nmodel; i++) {
      err=twoclass[i]->save(fs);
      if (err!=0) throw err;
    }
    return err;
  }

  template <typename real, typename cls_t>
  void multiclass<real, cls_t>::train(real **train, cls_t *cls, nel_ta ntrain, int type, real *param) {
    cls_t cls2[ntrain];

    cls_t k=0;
    for (cls_t i=0; i<nmodel; i++) {
      for (nel_ta j=0; j<ntrain; j++) {
        real mapel;
        if (cls[j]<0 || cls[j]>=nmodel) {
          cls2[j]=-1;
	} else {
          mapel=code[i][cls[j]];
	  if (mapel>0) {
            cls2[j]=1;
          } else if (mapel<0) {
            cls2[j]=0;
          } else {
            cls2[j]=-1;
	  }
	}
      }
      twoclass[i]->train(train, cls2, ntrain, type, param);
    }
  }

  template class multiclass<real_a, cls_ta>;

}

