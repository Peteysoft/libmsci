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

  template <class real, class cls_t>
  multiclass<real, cls_t>::multiclass(int ct) {
    type=ct;

    //set model variables to NULL:
    twoclass=NULL;
    u=NULL;
    vt=NULL;
    s=NULL;
    map=NULL;

    //set number parameters to 0:
    nmodel=0;
    this->ncls=0;

    //"polarity":
    pol=NULL;

  }

  template <class real, class cls_t>
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

    //"polarity":
    pol=NULL;
  }

  template <class real, class cls_t>
  int multiclass<real, cls_t>::init(multi_parse_param &param) {
    char *name[MAXNPART];
    cls_t *part[MAXNPART*2];
    int c1;
    char c2;
    int err=0;

    nmodel=parse_multi_partitions(&param, name, part, MAXNPART);
    err=multi_partition_strict(name, part, nmodel);
    if (err<0) throw PARAMETER_OUT_OF_RANGE;
    if (err>0) strictflag=0; else strictflag=1;

    //classification type and constraint weight:
    type=param.type;

    //pass to another initialization routine (duh...):
    init(name, part, nmodel, param.prefix, param.trainflag, param.commandname, 
		    param.Mflag, param.Kflag, param.sigcode, param.Zflag);

    //clean up:
    for (int i=0; i<nmodel; i++) {
      delete [] name[i];
      delete [] part[2*i];
      delete [] part[2*i+1];
    }

    fseek(param.infs, -1, SEEK_CUR);

    return err;

  }

  template <class real, class cls_t>
  int multiclass<real, cls_t>::init(char **fname, cls_t **part, int npart, char *prefix, 
		int tflag, char *com, int Mflag, int Kflag, int sigcode, int Zflag) {
    nmodel=npart;

    switch (type) {
      case (0):
        solve_class=&solve_class_constrained2<real>;
	break;
      case (1):
	solve_class=&solve_class_scratch<real>;
        break;
      case (2):
	solve_class=&solve_class_vote_pdf<real>;
        break;
      case (3):
	solve_class=&solve_class_vote<real>;
        break;
      case (4):
	solve_class=&solve_class_norm1<real>;
        break;
      case (5):
	solve_class=&solve_class_norm2<real>;
        break;
      case (6):
	solve_class=&solve_class_renorm<real>;
        break;
      case (7):
        solve_class=&solve_class_constrained1<real>;
	break;
      case (8):
	solve_class=&solve_class_vote_pdf2<real>;
	break;
      case (9):
	solve_class=&solve_class_1vR<real>;
	break;
	break;
      default:
        solve_class=&solve_class_constrained2<real>;
        break;
    }

    //figure out how many classes:
    this->ncls=0;
    for (int i=0; i<nmodel; i++) {
      for (int j=0; part[i*2][j]>=0; j++) {
        if (part[i*2][j]>=this->ncls) this->ncls=part[i*2][j]+1;
      }
      for (int j=0; part[i*2+1][j]>=0; j++) {
        if (part[i*2+1][j]>=this->ncls) this->ncls=part[i*2+1][j]+1;
      }
    }

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
        if (Zflag) {
          twoclass[i]=new svm2class<real, cls_t>(fname[i]);
	} else if (com==NULL) {
          //twoclass[i]=new borders_classifier<real, cls_t>(fname[i], sigcode);
          twoclass[i]=new borders_calibrated<real, cls_t>(fname[i]);
        } else {
          twoclass[i]=new general2class<real, cls_t>(fname[i], 
			com, Mflag, Kflag);
        }
      }
    }

    //create the mapping:
    map=gsl_matrix_alloc(nmodel+1, this->ncls);
    code=zero_matrix<real>(nmodel, this->ncls);	//two formats
    gsl_matrix_set_zero(map);
    //the sum of the conditional probabilities should always equal 1:
    for (int i=0; i<this->ncls; i++) gsl_matrix_set(map, nmodel, i, 1);

    //sum of the conditional probabilities on one side of the partition
    //is equal to the conditional probability returned from the 2-class
    //classification result:
    for (int i=0; i<nmodel*2; i++) {
      double sgn=2*(i%2)-1;
      for (int j=0; part[i][j]>=0; j++) {
        gsl_matrix_set(map, i/2, part[i][j], sgn);
	code[i/2][part[i][j]]=sgn;
      }
    }

    //find the inverse of this matrix when we need it (we don't yet):
    u=NULL;
    vt=NULL;
    s=NULL;

    return 0;

  }

  template <class real, class cls_t>
  multiclass<real, cls_t>::~multiclass() {
    //delete the binary classifiers and the decision matrix:
    for (int i=0; i<nmodel; i++) delete twoclass[i];
    delete [] twoclass;
    gsl_matrix_free(map);

    //delete the decomposition of the decision matrix:
    gsl_vector_free(s);
    gsl_matrix_free(vt);
    gsl_matrix_free(u);

    if (pol!=NULL) delete [] pol;

    delete_matrix(code);
  }

  //find the singular value decomposition of the coding matrix:
  template <class real, class cls_t>
  int multiclass<real, cls_t>::code_svd() {
    gsl_vector *work;
    int err;

    if (u!=NULL) return 0;

    //now we find the inverse of this matrix:
    u=gsl_matrix_alloc(nmodel+1, this->ncls);
    gsl_matrix_memcpy(u, map);

    //print_gsl_matrix(stdout, u);
    vt=gsl_matrix_alloc(this->ncls, this->ncls);
    s=gsl_vector_alloc(this->ncls);
    work=gsl_vector_alloc(this->ncls);
    err=gsl_linalg_SV_decomp(u, vt, s, work);
    //gsl_linalg_SV_decomp_jacobi(u, vt, s);

    /*
    printf("U:\n");
    print_gsl_matrix(stdout, u);
    printf("S:\n");
    for (int i=0; i<s->size; i++) printf("%10.5g ", gsl_vector_get(s, i));
    printf("\nV^T:\n");
    print_gsl_matrix(stdout, vt);
    */

    gsl_vector_free(work);

    return err;

  }

  //load a linear transformation to apply to the test points:
  template <class real, class cls_t>
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

  template <class real, class cls_t>
  cls_t multiclass<real, cls_t>::classify(real *x, real *p, real *praw) {
    real r[nmodel];
    real pt=0;

    if (praw!=NULL) {
      //printf("multiclass raw pdfs: ");
      for (int i=0; i<nmodel; i++) {
        r[i]=twoclass[i]->R(x, praw);
      }
    } else {
      for (int i=0; i<nmodel; i++) r[i]=twoclass[i]->R(x);
    }
    (*solve_class)(code, nmodel, this->ncls, r, p);
    //for (cls_t i=0; i<this->ncls; i++) pt+=p[i];
    //printf("pt=%g\n", pt);
    return choose_class(p, this->ncls);
  }

  template <class real, class cls_t>
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

  template <class real, class cls_t>
  void multiclass<real, cls_t>::batch_classify(real **x, cls_t *cls, real **p1, nel_ta n, dim_ta nvar) {
    real **r;

    //printf("multiclass: performing classifications with %d test vectors on %d models\n", n, nmodel);
    r=new real*[nmodel];
    r[0]=new real[nmodel*n];
    for (cls_t i=0; i<nmodel; i++) {
      r[i]=r[0]+i*n;
      //printf("multiclass: batchR model %d\n", i);
      twoclass[i]->batchR(x, r[i], n, nvar);
    }
    for (nel_ta i=0; i<n; i++) {
      (*solve_class)(code, nmodel, this->ncls, r[i], p1[i]);
      cls[i]=choose_class(p1[i], this->ncls);
    }

    delete [] r[0];
    delete [] r;

  }

  //for accessing "raw probabilities" for use in continuum predictions:
  template <class real, class cls_t>
  void multiclass<real, cls_t>::set_id(cls_t *id) {
    if (nmodel != this->ncls-1) {
      fprintf(stderr, "id model only valid if # binary classifiers = ncls-1\n");
      exit(PARAMETER_OUT_OF_RANGE);
    }
    for (cls_t i=0; i<nmodel; i++) {
      twoclass[i]->set_id(id+i);
      //remove side-effect:
      id[i]--;
    }
  }

  template <class real, class cls_t>
  void multiclass<real, cls_t>::print(FILE *fs, char *fbase, int depth) {
    char *fbase2=NULL;
    if (fbase!=NULL) fbase2=new char[strlen(fbase)+5];

    for (int i=0; i<nmodel; i++) {
      if (fbase!=NULL) sprintf(fbase2, "%s-%2.2d", fbase, i);
      twoclass[i]->print(fs, fbase2, depth);
      for (int j=0; j<this->ncls; j++) {
        if (gsl_matrix_get(map, i, j)<0) fprintf(fs, " %d", j);
      }
      fprintf(fs, " /");
      for (int j=0; j<this->ncls; j++) {
        if (gsl_matrix_get(map, i, j)>0) fprintf(fs, " %d", j);
      }
      fprintf(fs, ";");
      if (i!=nmodel-1) fprintf(fs, "\n");
    }

    if (fbase2!=NULL) delete [] fbase2;
  }

  template <class real, class cls_t>
  int multiclass<real, cls_t>::commands(multi_train_param &param,
                cls_t **clist, char *fbase) {
    double coef;
    char *fbase2;
    cls_t *clist2;
    cls_t *clist3[3];
    cls_t nct=0;		//number of classes total
    cls_t nc1, nc2;		//number in first, second partition resp.

    //how many classes total?
    for (cls_t i=0; clist[0]+i != clist[this->ncls]; i++) nct++;
    clist2=new cls_t[nct+1];		//need that extra element hanging off
					//the end as a flag for stopping iteration

    fbase2=new char[strlen(fbase)+4];

    for (int i=0; i<nmodel; i++) {
      //generate file name:
      sprintf(fbase2, "%s-%2.2d", fbase, i);
      //print command name, options, training data, output model:
      //fprintf(param.commandfs, "%s %s %s %s", param.commandname, name[i],
	//	param.train, fbase2);

      //print class partitions:
      //gather the class labels in each partition:
      nc1=0;
      clist3[0]=clist2;
      for (cls_t j=0; j<this->ncls; j++) {
        coef=gsl_matrix_get(map, i, j);
        if (coef<0) {
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
        coef=gsl_matrix_get(map, i, j);
        if (coef>0) {
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
    delete [] fbase2;
    delete [] clist2;

    return this->ncls;
  }

  //converts matrix to a sorted array of STL vectors:
  template <typename scalar>
  void matrix2sorted(scalar **mat,		//original matrix
			int m,
			int n, 
			vector<scalar> *sd,	//vector of sorted arrays
			long *sind) {		//sorting indices
    vector<scalar> sd2[m];

    for (int i=0; i<m; i++) {
      sd[i].resize(n);		//I thought these fucking things were supposed
      sd2[i].resize(n);		//to resize themselves??
      for (int j=0; j<n; j++) sd[i][j]=mat[i][j];
    }

    heapsort(sd, sind, m);

    for (int i=0; i<m; i++) sd2[i]=sd[sind[i]];
    for (int i=0; i<m; i++) sd[i]=sd2[i];
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
      for (int j=0; j<this->ncls; j++) sd[i][j]=gsl_matrix_get(map, i, j);
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
      gsl_matrix *newmap=gsl_matrix_alloc(nmodel+1, this->ncls);
      binaryclassifier<real, cls_t> **twoclass2=new binaryclassifier<real, cls_t>*[nmodel];
      for (int i=0; i<nmodel; i++) {
        for (int j=0; j<this->ncls; j++) {
          gsl_matrix_set(newmap, i, j, gsl_matrix_get(map, ind[i], j));
        }
	twoclass2[i]=twoclass[ind[i]];
      }
      for (int j=0; j<this->ncls; j++) {
        gsl_matrix_set(newmap, nmodel, j, gsl_matrix_get(map, nmodel, j));
      }
      gsl_matrix_free(map);
      delete [] twoclass;
      map=newmap;
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
      if (type<0) type=10;
      strictflag=1;
    } else if (strcmp(typestr, "1v1")==0) {
      nmodel=this->ncls*(this->ncls-1)/2;
      code=one_against_one<real>(this->ncls);
      if (type<0) type=11;
      strictflag=0;
    } else if (strcmp(typestr, "ADJ")==0) {
      nmodel=this->ncls-1;
      code=partition_adjacent<real>(this->ncls);
      if (type<0) type=12;
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
      twoclass[i]=new borders_calibrated<real, cls_t>();
      err=twoclass[i]->load(fs);
      if (err!=0) {
        fprintf(stderr, "multiclass::load: error loading border vectors\n");
        throw err;
      }
    }

    //convert integer coding matrix to floating point, GSL compatible one:
    map=gsl_matrix_alloc(nmodel+1, this->ncls);
    for (int i=0; i<nmodel; i++) {
      for (int j=0; j<this->ncls; j++) gsl_matrix_set(map, i, j, code[i][j]);
    }
    for (int j=0; j<this->ncls; j++) gsl_matrix_set(map, nmodel, j, 1);

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
    cls_t *map2;			//for partitioning the classes
    cls_t cls2[ntrain];

    cls_t k=0;
    for (cls_t i=0; i<nmodel; i++) {
      for (nel_ta j=0; j<ntrain; j++) {
        real mapel;
        if (cls[j]<0 || cls[j]>=nmodel) {
          cls2[j]=-1;
	} else {
          mapel=gsl_matrix_get(map, i, cls[j]);
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

