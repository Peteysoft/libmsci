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
    constraint_weight=1;

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
  multiclass<real, cls_t>::multiclass(const char *file, int clstyp, const char *com, int mf, int kf, int sigcode) {
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

    //weight for normalization constraint:
    param.cw=1;

    param.commandname=NULL;
    if (com!=NULL) {
      param.commandname=new char [strlen(com)+1];
      strcpy(param.commandname, com);
    }
    param.Mflag=mf;
    param.Kflag=kf;
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
    constraint_weight=param.cw;

    init(name, part, nmodel, param.trainflag, param.commandname, 
		    param.Mflag, param.Kflag, param.sigcode);

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
  int multiclass<real, cls_t>::init(char **fname, cls_t **part, int npart, 
		int tflag, char *com, int Mflag, int Kflag, int sigcode) {
    nmodel=npart;
    this->ncls=0;
    for (int i=0; i<nmodel; i++) {
      for (int j=0; part[i*2][j]>=0; j++) {
        if (part[i*2][j]>=this->ncls) this->ncls=part[i*2][j]+1;
      }
      for (int j=0; part[i*2+1][j]>=0; j++) {
        if (part[i*2+1][j]>=this->ncls) this->ncls=part[i*2+1][j]+1;
      }
    }

    twoclass=new binaryclassifier<real, cls_t> *[nmodel];
    for (cls_t i=0; i<nmodel; i++) {
      if (tflag) {
        twoclass[i]=new binaryclassifier<real, cls_t>(fname[i]);
      } else {
        if (com==NULL) {
          twoclass[i]=new agf2class<real, cls_t>(fname[i], sigcode);
        } else {
          twoclass[i]=new general2class<real, cls_t>(fname[i], 
			com, Mflag, Kflag);
        }
      }
    }

    //create the mapping:
    map=gsl_matrix_alloc(nmodel+1, this->ncls);
    gsl_matrix_set_zero(map);
    //the sum of the conditional probabilities should always equal 1:
    for (int i=0; i<this->ncls; i++) gsl_matrix_set(map, nmodel, i, constraint_weight);

    //sum of the conditional probabilities on one side of the partition
    //is equal to the conditional probability returned from the 2-class
    //classification result:
    for (int i=0; i<nmodel*2; i++) {
      double sgn=2*(i%2)-1;
      for (int j=0; part[i][j]>=0; j++) gsl_matrix_set(map, i/2, part[i][j], sgn);
    }

    //find the inverse of this matrix when we need it (we don't yet):
    u=NULL;
    vt=NULL;
    s=NULL;

    //initialize constraints when we need them:
    cnorm=NULL;
    cthresh=NULL;

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

    if (cnorm!=NULL) {
      gsl_matrix_free(cnorm);
      gsl_vector_free(cthresh);
    }

    if (pol!=NULL) delete [] pol;
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

  //initialize the constraint coefficents and thresholds:
  template <class real, class cls_t>
  void multiclass<real, cls_t>::init_constraint() {
    if (cnorm==NULL) {
      gsl_vector_view lastrow;
      cnorm=gsl_matrix_alloc(this->ncls, this->ncls-1);
      gsl_matrix_set_identity(cnorm);
      cthresh=gsl_vector_alloc(this->ncls);
      gsl_vector_set_zero(cthresh);
      lastrow=gsl_matrix_row(cnorm, this->ncls-1);
      gsl_vector_set_all(&lastrow.vector, -1);
      gsl_vector_set(cthresh, this->ncls-1, -1);
    }
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

  //vote based on class values from raw probabilities:
  template <class real, class cls_t>
  cls_t multiclass<real, cls_t>::vote_label(gsl_vector *b, real *tly) {
    real val;
    cls_t cls;
    //we just inline it:
    for (int i=0; i<this->ncls; i++) tly[i]=0;
    for (int i=0; i<nmodel; i++) {
      val=gsl_vector_get(b, i);
      for (int j=0; j<this->ncls; j++) {
        tly[j]+=gsl_matrix_get(map, i, j)*val/fabs(val);
      }
    }
    //correction isn't perfect, but should move tallies closer to 
    //conditional probability:
    for (int j=0; j<this->ncls; j++) tly[j]=(tly[j]+1)/(nmodel+1);
    cls=choose_class(tly, this->ncls);
    return cls;
  }

  //vote based on probabilities from the binary classifier:
  template <class real, class cls_t>
  cls_t multiclass<real, cls_t>::vote_pdf(gsl_vector *b, real *tly) {
    real val;
    cls_t cls;
    //we just inline it:
    for (int i=0; i<this->ncls; i++) tly[i]=0;
    for (int i=0; i<nmodel; i++) {
      val=gsl_vector_get(b, i);
      for (int j=0; j<this->ncls; j++) {
        tly[j]+=gsl_matrix_get(map, i, j)*val;
      }
    }
    cls=choose_class(tly, this->ncls);
    return cls;
  }

  //vote based on probabilities from the binary classifier:
  //corrected and re-normalized:
  template <class real, class cls_t>
  cls_t multiclass<real, cls_t>::vote_pdf2(gsl_vector *b, real *tly) {
    real val;
    cls_t cls;
    real pt=0;			//total of computed cond. prob.
    cls_t ind[this->ncls];	//flag the good values
    cls_t ng=0;			//number of good values
    cls_t ng2;
    float k;			//additive normalization constant

    //we just inline it:
    for (int i=0; i<this->ncls; i++) tly[i]=0;
    for (int i=0; i<nmodel; i++) {
      val=gsl_vector_get(b, i);
      for (int j=0; j<this->ncls; j++) {
        tly[j]+=gsl_matrix_get(map, i, j)*val;
      }
    }
    //must be done in two steps (I think...):
    //correction step:
    for (int i=0; i<this->ncls; i++) {
      tly[i]=(tly[i]+1)/(nmodel+1);
      if (tly[i]<0) {
        tly[i]=0;
      } else {
        pt+=tly[i];
	ind[ng]=i;
	ng++;
      }
    }
    //re-normalization step (farm out to another unit):
    //(use this version since we assume:
    //- A^T*A=nI where A is the coding matrix
    //- top row of A is all 1's which is not explicitly included in the control
    //  file but is included in above calculation
    //- we solve: p = A^T*r where the first value in r is a free parameter 
    //  which we vary for normalization if the other constraints in p are 
    //  violated)
    p_renorm3(tly, this->ncls);

    cls=choose_class(tly, this->ncls);
    return cls;
  }

  //gets results from all the binary models:
  template <class real, class cls_t>
  void multiclass<real, cls_t>::raw_classify(real *x, gsl_vector *b) {
    //printf("multiclass raw pdfs: ");
    for (int i=0; i<nmodel; i++) {
      gsl_vector_set(b, i, twoclass[i]->R(x));
    }
    //printf("\n");
    gsl_vector_set(b, nmodel, constraint_weight);
  }

  //classification from pseudo-inverse (linear least squares via SVD):
  template <class real, class cls_t>
  cls_t multiclass<real, cls_t>::classify_basic(gsl_vector *b, real *p) {
    cls_t cls;
    gsl_vector *p1=gsl_vector_alloc(this->ncls);

    if (strictflag==0) classify_scratch(b, p);

    code_svd();

    gsl_linalg_SV_solve(u, vt, s, b, p1);
    for (cls_t i=0; i<this->ncls; i++) p[i]=gsl_vector_get(p1, i);
    cls=choose_class(p, this->ncls);
    gsl_vector_free(p1);

    return cls;
  }

  //solve matrix equation from scratch so we can apply regularizations to it:
  template <class real, class cls_t>
  cls_t multiclass<real, cls_t>::classify_scratch(gsl_vector *b, real *p) {
    cls_t cls;
    gsl_vector *s1;
    gsl_matrix *u1;
    gsl_matrix *vt1;
    gsl_vector *work;
    gsl_vector *b1;

    gsl_vector *p1=gsl_vector_alloc(this->ncls);

    vt1=gsl_matrix_alloc(this->ncls, this->ncls);
    s1=gsl_vector_alloc(this->ncls);
    work=gsl_vector_alloc(this->ncls);

    //copy mapping:
    u1=gsl_matrix_alloc(nmodel+1, this->ncls);
    gsl_matrix_memcpy(u1, map);

    //normalize on the basis of the raw probabilities:
    b1=gsl_vector_alloc(nmodel+1);
    gsl_vector_set(b1, nmodel, 1);
    for (int i=0; i<nmodel; i++) {
      real val=gsl_vector_get(b, i);
      for (int j=0; j<this->ncls; j++) {
        real map_el=gsl_matrix_get(u1, i, j);
	if (map_el==0) {
          //gsl_matrix_set(u1, i, j, 1);
          gsl_matrix_set(u1, i, j, val);
	} else {
          //gsl_matrix_set(u1, i, j, map_el/val);
          gsl_matrix_set(u1, i, j, 1);
	}
      }
      gsl_vector_set(b1, i, 1);
    }

    gsl_linalg_SV_decomp(u1, vt1, s1, work);

    gsl_linalg_SV_solve(u1, vt1, s1, b, p1);
    for (cls_t i=0; i<this->ncls; i++) p[i]=gsl_vector_get(p1, i);
    cls=choose_class(p, this->ncls);

    //clean up:
    gsl_vector_free(p1);
    gsl_vector_free(work);
    gsl_vector_free(s1);
    gsl_matrix_free(u1);
    gsl_matrix_free(vt1);
    gsl_vector_free(b1);

    return cls;
  }

  //"constrained" version (not completely rigorous, but works quite well):
  int solve_cond_prob(gsl_matrix *, gsl_vector *, gsl_vector *);

  template <class real, class cls_t>
  cls_t multiclass<real, cls_t>::classify_special(gsl_vector *b, real *p) {
    real pt=0;
    gsl_matrix *map1;

    if (strictflag) {
      map1=map;
    } else {
      map1=gsl_matrix_alloc(nmodel+1, this->ncls);
      gsl_matrix_memcpy(map1, map);
      for (cls_t i=0; i<nmodel; i++) {
        cls_t cnt=0;
	real r_i=gsl_vector_get(b, i);
        for (cls_t j=0; j<this->ncls; j++) {
          real map_el=gsl_matrix_get(map, i, j);
          if (map_el==0) gsl_matrix_set(map1, i, j, r_i);
	}
      }
    }

    //old method:
    if (type==7) {	//(logic isn't pretty, but simplest way to implement)
      gsl_vector *p1=gsl_vector_alloc(this->ncls);
      solve_cond_prob(map1, b, p1);
      for (cls_t i=0; i<this->ncls; i++) {
        p[i]=gsl_vector_get(p1, i);
      }
      gsl_vector_free(p1);
    //new method:
    } else {
      gsl_vector *bt=gsl_vector_alloc(nmodel);
      gsl_vector *p2=gsl_vector_alloc(this->ncls-1);
      //apply normalization constraint:  
      //to avoid any biases produced by using the same variable each time
      int ind=ranu()*this->ncls;
      gsl_matrix *at=gsl_matrix_alloc(nmodel, this->ncls-1);
      init_constraint();

      for (int i=0; i<nmodel; i++) {
        double aind=gsl_matrix_get(map1, i, ind);
        gsl_vector_set(bt, i, gsl_vector_get(b, i)-aind);
        for (int j=0; j<ind; j++) {
          gsl_matrix_set(at, i, j, gsl_matrix_get(map1, i, j)-aind);
        }
        for (int j=ind+1; j<this->ncls; j++) {
          gsl_matrix_set(at, i, j-1, gsl_matrix_get(map1, i, j)-aind);
        }
      }

      constrained(at, bt, cnorm, cthresh, p2);

      //reconstitute missing variable and extract the rest:
      p[ind]=1;
      for (int j=0; j<ind; j++) {
        p[j]=gsl_vector_get(p2, j);
        p[ind]-=p[j];
      }
      for (int j=ind+1; j<this->ncls; j++) {
        p[j]=gsl_vector_get(p2, j-1);
        p[ind]-=p[j];
      }
      gsl_vector_free(p2);
      gsl_vector_free(bt);
      gsl_matrix_free(at);
    }

    if (strictflag!=1) gsl_matrix_free(map1);

    return choose_class(p, this->ncls);

  }

  template <class real, class cls_t>
  cls_t multiclass<real, cls_t>::classify_1vR(gsl_vector *b, real *p) {
    int ind[this->ncls];		//still being worked on
    int nwork;				//how many prob. still being worked on
    int k=ranu()*this->ncls;		//omit random classifier
    int dflag;

    for (cls_t i=0; i<this->ncls; i++) p[i]=pol[i]*(1-gsl_vector_get(b, i))/2;

    for (cls_t i=0; i<k; i++) ind[i]=i;
    for (cls_t i=k+1; i<this->ncls; i++) ind[i-1]=i;
    nwork=this->ncls-1;

    do {
      real pt;
      dflag=1;
      for (cls_t i=0; i<nwork; i++) {
        if (p[ind[i]]<0) {
          p[ind[i]]=0;
	  for (cls_t j=i+1; i<nwork; j++) ind[i-1]=ind[i];
	  nwork--;
	  dflag=0;
	  break;
	}
      }
      if (nwork==0) break;
      pt=0;
      for (cls_t i=0; i<nwork; i++) pt+=p[ind[i]];
      if (pt>1) {
        for (cls_t i=0; i<nwork; i++) p[ind[i]]=p[ind[i]]-(pt-1)/nwork;
	dflag=0;
      }
    } while (dflag==0);

    p[k]=1;
    for (int i=0; i<nwork; i++) p[k]-=p[ind[i]];

    return choose_class(p, this->ncls);
  }

  template <class real, class cls_t>
  cls_t multiclass<real, cls_t>::classify_1v1(gsl_vector *b, real *p) {
    real **praw=allocate_matrix<real, int32_t>(this->ncls, this->ncls);
    int k=0;
    for (cls_t i=0; i<this->ncls; i++) {
      for (cls_t j=i+1; j<this->ncls; j++) {
        praw[i][j]=(1-pol[k]*gsl_vector_get(b, k))/2;
	k++;
      }
    }
    solve_cond_prob_1v1(praw, this->ncls, p);
    delete_matrix(praw);
    return choose_class(p, this->ncls);
  }

  template <class real, class cls_t>
  cls_t multiclass<real, cls_t>::solve_class(gsl_vector *b, real *pdf) {
    cls_t cls1, cls2;
    real pt;		//total of conditional probability estimates
    real k;		//correction value
    real tally[this->ncls];
    char flag[this->ncls];
    int nc2=0;

    cls2=-1;

    //printf("classification type=%d\n", type);
    switch(type) {
      case (0):
	cls1=classify_special(b, pdf);
	break;
      case (1):
        cls1=classify_basic(b, pdf);
        break;
      case (2):
        cls1=vote_pdf(b, pdf);
        break;
      case (3):
        cls1=vote_label(b, pdf);
        break;
      case (4):
	cls1=classify_basic(b, pdf);
        cls2=vote_pdf(b, tally);
        break;
      case (5):
        cls1=classify_basic(b, pdf);
        cls2=vote_label(b, tally);
        break;
      case (6):
        cls1=classify_scratch(b, pdf);
        break;
      case (7):
	cls1=classify_special(b, pdf);
	break;
      case (8):
        cls1=vote_pdf2(b, pdf);
	break;
      case (10):
	cls1=classify_1vR(b, pdf);
	break;
      case (11):
	cls1=classify_1v1(b, pdf);
	break;
      default:
        cls1=classify_special(b, pdf);
        break;
    }

    if (cls2 >= 0) {
      //correct the resultant conditional probabilities (=hack):
      //(we use these specific forms of renormalization because they tend to
      //maximize the "peakedness" of the distribution which seem to better
      //reflect most real distributions)
      p_renorm1(pdf, this->ncls);

      //if voting is different from matrix inversion, correct the results using a crude hack:
      if (cls1!=cls2) {
        cls_t swp=cls1;
        cls1=cls2;
        cls2=swp;
        pt=pdf[cls2];
        for (cls_t i=0; i<this->ncls; i++) if (i!=cls1) pt+=pdf[i];
        k=(pdf[cls2]-pdf[cls1])/pt;
        for (cls_t i=0; i<this->ncls; i++) if (i!=cls1) pdf[i]=(1-k)*pdf[i];
        pdf[cls1]=pdf[cls2];
        //pt=0;
        //for (cls_t i=0; i<this->ncls; i++) pt+=pdf[i];
        //printf("pt (2)=%g\n", pt);
      }
    }
    pt=0;
    for (cls_t i=0; i<this->ncls; i++) pt+=pdf[i];
    printf("pt (2)=%g\n", pt);

    return cls1;
  }

  template <class real, class cls_t>
  cls_t multiclass<real, cls_t>::classify(real *x, real *pdf, real *praw) {
    cls_t cls;
    gsl_vector *b=gsl_vector_alloc(nmodel+1);

    if (praw!=NULL) {
      //printf("multiclass raw pdfs: ");
      for (int i=0; i<nmodel; i++) {
        gsl_vector_set(b, i, twoclass[i]->R(x, praw));
      }
      gsl_vector_set(b, nmodel, constraint_weight);
    } else {
      raw_classify(x, b);
    }
    cls=solve_class(b, pdf);
    gsl_vector_free(b);
    return cls;
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
    gsl_vector *b=gsl_vector_alloc(nmodel+1);

    //printf("multiclass: performing classifications with %d test vectors on %d models\n", n, nmodel);
    r=new real*[nmodel];
    r[0]=new real[nmodel*n];
    for (cls_t i=0; i<nmodel; i++) {
      r[i]=r[0]+i*n;
      //printf("multiclass: batchR model %d\n", i);
      twoclass[i]->batchR(x, r[i], n, nvar);
    }
    for (nel_ta i=0; i<n; i++) {
      for (int j=0; j<nmodel; j++) {
        gsl_vector_set(b, j, r[j][i]);
      }
      gsl_vector_set(b, nmodel, constraint_weight);
      cls[i]=solve_class(b, p1[i]);
    }

    delete [] r[0];
    delete [] r;
    gsl_vector_free(b);

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
      for (int j=0; j<n; j++) sd[i][j]=mat[i][j];
    }

    heapsort(sd, sind, m);

    for (int i=0; i<m; i++) sd2[i]=sd[sind[i]];
    for (int i=0; i<m; i++) sd[i]=sd2[i];
  }

  template <typename real, typename cls_t>
  int multiclass<real, cls_t>::detect_type() {
    int spec_type;		//0=1 v. rest; 1=1 v. 1; 2=adj.
    int **code;				//coding matrix (as integers)
    int **proto;			//prototype to compare against
    long sind[this->nmodel];
    long dum[nmodel];
    vector<int> sd[nmodel];
    vector<int> p2[nmodel];

    code=allocate_matrix<int, int32_t>(nmodel, this->ncls);
    for (int i=0; i<nmodel; i++) {
      for (int j=0; j<this->ncls; j++) code[i][j]=gsl_matrix_get(map, i, j);
    }

    matrix2sorted(code, nmodel, this->ncls, sd, sind);

    spec_type=-1;			//which type is it?
    pol=new int[nmodel];
    for (int tc=0; tc<3; tc++) {
      int nm;			//required number of binary classifiers
      switch (tc) {
        case(0):
          //check for 1 v rest:
          proto=one_against_all(this->ncls);
	  nm=this->ncls;
	  break;
        case(1):
          //check for 1 v 1:
          proto=one_against_one(this->ncls);
	  nm=this->ncls*(this->ncls-1)/2;
	  break;
        case(2):
	  //check for adjacent partitioning:
          proto=partition_adjacent(this->ncls);
	  nm=this->ncls-1;
	  break;
      }
      if (nmodel==nm) {
        matrix2sorted(proto, this->nmodel, this->ncls, p2, dum);
        spec_type=tc;
        for (int i=0; i<this->ncls; i++) {
          if (code[i]!=proto[i]) {
            for (int j=0; j<this->ncls; j++) proto[i][j]=-proto[i][j];
            if (code[i]!=proto[i]) {
              spec_type=-1;
              break;
            } else {
              pol[i]=-1;
	    }
          } else {
            pol[i]=1;
	  }
	}
      }
      delete [] proto[0];
      delete [] proto;
      if (spec_type==-1) {
        delete [] pol;
	pol=NULL;
        break;
      }
    }
    if (spec_type!=-1) {
      //rearrange the coding matrix and binary classifiers so they're in
      //the "right" order:
      gsl_matrix *newmap=gsl_matrix_alloc(nmodel+1, this->ncls);
      binaryclassifier<real, cls_t> **twoclass2=new binaryclassifier<real, cls_t>*[nmodel];
      for (int i=0; i<nmodel; i++) {
        for (int j=0; j<this->ncls; j++) {
          gsl_matrix_set(newmap, dum[i], j, gsl_matrix_get(map, sind[i], j));
        }
	twoclass2[dum[i]]=twoclass2[sind[i]];
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

    //clean up:
    delete [] code[0];
    delete [] code;
      
    return spec_type;
  }

  template <typename real, typename cls_t>
  int multiclass<real, cls_t>::load(FILE *fs) {
    int err=0;
    int **code;				//coding matrix (as integer)
    char **sub;
    int nsub;
    char *typestr=fget_line(fs, 1);
    char *line=fget_line(fs);
    sscanf(line, "%d", &this->ncls);
    delete [] line;
    //type specific stuff:
    if (strcmp(typestr, "1vR")==0) {
      nmodel=this->ncls;
      code=one_against_all(this->ncls);
      if (type<0) type=10;
    } else if (strcmp(typestr, "1v1")==0) {
      nmodel=this->ncls*(this->ncls-1)/2;
      code=one_against_one(this->ncls);
      if (type<0) type=11;
    } else if (strcmp(typestr, "ADJ")==0) {
      nmodel=this->ncls-1;
      code=partition_adjacent(this->ncls);
      if (type<0) type=12;
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
      for (int j=0; j<this->ncls; j++) code[i][j]=-pol[i]*code[i][j];
    }
    delete [] line;
    delete [] sub;

    //read in binary classifiers:
    twoclass=new binaryclassifier<real, cls_t>*[nmodel];
    for (int i=0; i<nmodel; i++) {
      twoclass[i]=new agf2class<real, cls_t>();
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
    constraint_weight=1;

    print_gsl_matrix(stdout, map);
    
    return err;
  }

  template <typename real, typename cls_t>
  int multiclass<real, cls_t>::save(FILE *fs) {
    int err=0;
    detect_type();
    switch (type) {
      case(10):
        fprintf(fs, "1vR");
	break;
      case(11):
	fprintf(fs, "1v1");
	break;
      case(12):
	fprintf(fs, "ADJ");
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

  template class multiclass<real_a, cls_ta>;

}

