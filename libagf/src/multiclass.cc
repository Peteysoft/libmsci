#include <stdio.h>
#include <getopt.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <assert.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>

#include "constrained.h"
#include "full_util.h"
#include "gsl_util.h"
//#include "peteys_tmpl_lib.h"

#include "agf_lib.h"

using namespace std;
using namespace libpetey;

namespace libagf {

  template <class real, class cls_t>
  multiclass<real, cls_t>::multiclass() {
    constraint_weight=1;
    type=2;

    //set model variables to NULL:
    twoclass=NULL;
    u=NULL;
    vt=NULL;
    s=NULL;
    map=NULL;
    imap=NULL;

    //set number parameters to 0:
    nmodel=0;
    this->ncls=0;
  }

  template <class real, class cls_t>
  multiclass<real, cls_t>::multiclass(const char *file, int clstyp, real cw, const char *com, int mf, int kf, int sigcode) {
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
    constraint_weight=cw;

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

    type=clstyp;

  }

  template <class real, class cls_t>
  int multiclass<real, cls_t>::init(multi_parse_param &param) {
    char *name[MAXNPART];
    cls_t *part[MAXNPART*2];
    int c1;
    char c2;
    int flag;
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

    //check for a mapping:
    flag=0;
    do {
      c1=fgetc(param.infs);
      if (c1==EOF) break;
      c2=(char) c1;
      //printf("1.(%c)\n", c2);
      if (c2=='\n') param.lineno++;
      //filename has to start with a letter:
      if (c2=='[') {
        flag=1;
        break;
      }
    } while (isspace(c1));

    if (flag) {
      int32_t m1, n1;
      imap=scan_matrix<real, int32_t>(param.infs, m1, n1);
      if (imap==NULL) {
        fprintf(stderr, "multiclass::init: error parsing mapping; skipping\n");
	err=FILE_READ_ERROR;
      }
      if (m1 != this->ncls || n1 != nmodel) {
        fprintf(stderr, "multiclass::init: map has wrong dimensions: [%d*%d]; [%d*%d] expected\n", m1, n1, this->ncls, nmodel);
	delete_matrix(imap);
	err=SAMPLE_COUNT_MISMATCH;
      }

      //find closing bracket:
      flag=0;
      do {
        c1=fgetc(param.infs);
        if (c1==EOF) break;
        c2=(char) c1;
        //printf("1.(%c)\n", c2);
        if (c2=='\n') param.lineno++;
        //filename has to start with a letter:
        if (c2==']') {
          flag=1;
          break;
        }
      } while (isspace(c1));

      if (flag!=1) {
        fprintf(stderr, "multiclass::init: syntax error on line %d;\n", param.lineno);
        fprintf(stderr, "     no closing bracket for mapping; exiting...\n");
	exit(SYNTAX_ERROR);
      }
    } else {
      fseek(param.infs, -1, SEEK_CUR);
    }

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

    if (imap!=NULL) {
      delete [] imap[0];
      delete [] imap;
    }
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
    printf("multiclass: initializing constraints\n");
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

  //train mapping with actual data:
  template <class real, class cls_t>
  int multiclass<real, cls_t>::train_map(real **train, cls_t *cls, nel_ta n) {
    gsl_matrix *a;
    gsl_vector *b;
    gsl_vector *p;
    gsl_vector *x;

    gsl_matrix *vt2;
    gsl_vector *s2;
    gsl_vector *work;
    cls_t weight[this->ncls];
    real_a sump[nmodel];

    for (cls_t i=0; i<this->ncls; i++) weight[i]=0;
    for (cls_t i=0; i<nmodel; i++) sump[i]=0;
    for (nel_ta i=0; i<n; i++) weight[cls[i]]++;

    a=gsl_matrix_alloc(n*(this->ncls+1)+this->ncls, nmodel*this->ncls);
    b=gsl_vector_alloc(n*(this->ncls+1)+this->ncls);
    p=gsl_vector_alloc(nmodel+1);
    for (nel_ta i=0; i<n; i++) {
      raw_classify(train[i], p);
      for (cls_ta j=0; j<this->ncls; j++) {
        for (cls_ta k=0; k<nmodel; k++) {
          gsl_matrix_set(a, i*(this->ncls+1)+j, j*nmodel+k, 
			  gsl_vector_get(p, k)/weight[j]);
	  gsl_matrix_set(a, i*(this->ncls+1)+this->ncls, j*nmodel+k, 
			  gsl_vector_get(p, k));
        }
        if (cls[i]==j) {
          gsl_vector_set(b, i*(this->ncls+1)+j, 1./weight[j]);
        } else {
          gsl_vector_set(b, i*(this->ncls+1)+j, 0);
        }
	gsl_vector_set(b, i*(this->ncls+1)+this->ncls, 1.);
      }
      for (cls_ta k=0; k<nmodel; k++) {
        sump[k]+=gsl_vector_get(p, k);
      }
    }

    //constraint: <p_i>=n_i/n
    for (int i=0; i<this->ncls; i++) {
      for (int j=0; j<nmodel; j++) {
        gsl_matrix_set(a, n*(this->ncls+1)+i, i*nmodel+j, sump[j]/weight[i]);
      }
      gsl_vector_set(b, n*(this->ncls+1)+i, 1.);
    }

    print_gsl_matrix(stdout, a);

    //solve the linear system:
    x=gsl_vector_alloc(nmodel*this->ncls);
    vt2=gsl_matrix_alloc(nmodel*this->ncls, nmodel*this->ncls);
    s2=gsl_vector_alloc(nmodel*this->ncls);
    work=gsl_vector_alloc(nmodel*this->ncls);

    gsl_linalg_SV_decomp(a, vt2, s2, work);
    gsl_linalg_SV_solve(a, vt2, s2, b, x);

    imap=allocate_matrix<real, int>(this->ncls, nmodel);
    for (cls_t i=0; i<this->ncls; i++) {
      for (cls_t j=0; j<nmodel; j++) {
        imap[i][j]=gsl_vector_get(x, i*nmodel+j);
      }
    }
    //need some error checking here...
    return 0;
  }

  template <class real, class cls_t>
  int multiclass<real, cls_t>::write_map(FILE *fs) {
    return write_matrix<real, int32_t>(fs, imap, nmodel, this->ncls);
  }

  template <class real, class cls_t>
  int multiclass<real, cls_t>::load_map(FILE *fs) {
    int32_t m1, n1;
    if (imap!=NULL) delete_matrix(imap);
    imap=read_matrix<real, int32_t>(fs, m1, n1);
    if (imap==NULL) return FILE_READ_ERROR;
    if (m1!=nmodel || n1!=this->ncls) {
      fprintf(stderr, "multiclass::load_map: failed to load map;\n");
      fprintf(stderr, "     dimensions incorrect: [%d * %d]; expected [%d * %d]\n",
		      m1, n1, nmodel, this->ncls);
      return SAMPLE_COUNT_MISMATCH;
    }
    return 0;
  }

  //load a linear transformation to apply to the test points:
  template <class real, class cls_t>
  int multiclass<real, cls_t>::ltran(real **mat, real *b, dim_ta d1, dim_ta d2, int flag) {
    int err=0;
    //just pass it to the binary classifiers:
    for (int i=0; i<nmodel; i++) {
      err=twoclass[i]->ltran(mat, b, d1, d2, flag);
      if (err!=0) {
        fprintf(stderr, "multiclass::ltran: an error occured transforming partition #%d\n", i);
        return err;
      }
    }
    this->D=-1;
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

  //classify from a mapping (inverse) constructed earlier:
  template <class real, class cls_t>
  cls_t multiclass<real, cls_t>::classify_map(gsl_vector *b, real *p) {
    for (cls_t i=0; i<this->ncls; i++) {
      p[i]=0;
      for (cls_t j=0; j<nmodel; j++) p[i]+=gsl_vector_get(b, j)*imap[i][j];
    }
    return choose_class(p, this->ncls);
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
    gsl_vector *bt=gsl_vector_alloc(nmodel);
    //solve by old method:
    gsl_vector *p1=gsl_vector_alloc(this->ncls);
    //solve by new method:
    gsl_vector *p2=gsl_vector_alloc(this->ncls-1);

    init_constraint();

    //print_gsl_matrix(stdout, map);
    //printf("\n");

    if (strictflag) {
      solve_cond_prob(map, b, p1);
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
      solve_cond_prob(map1, b, p1);
      gsl_matrix_free(map1);
    }

    //apply normalization constraint:  
    //to avoid any biases produced by using the same variable each time
    int ind=ranu()*this->ncls;
    gsl_matrix *at=gsl_matrix_alloc(nmodel, this->ncls-1);

    printf("multiclass::classify_special: applying first constraint\n");
    for (int i=0; i<nmodel; i++) {
      printf("multiclass::classify_special: i=%d\n", i);
      double aind=gsl_matrix_get(map1, i, ind);
      gsl_vector_set(bt, i, gsl_vector_get(b, i)-aind);
      for (int j=0; j<ind; j++) {
        gsl_matrix_set(at, i, j, gsl_matrix_get(map1, i, j)-aind);
      }
      for (int j=ind+1; j<this->ncls; j++) {
        gsl_matrix_set(at, i, j-1, gsl_matrix_get(map1, i, j)-aind);
      }
    }

    printf("multiclass::classify_special: calling constrained subroutine\n");
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

    for (cls_t i=0; i<this->ncls; i++) {
      printf("%g ", gsl_vector_get(p1, i));
      p[i]=gsl_vector_get(p1, i);
      //if (p[i]<0) printf("p[%d]=%g out-of-bounds\n", i, p[i]);
      //pt+=p[i];
    }
    printf("\n");
    //printf("pt=%g\n", pt);
    gsl_vector_free(p1);
    gsl_vector_free(p2);
    if (strictflag!=1) gsl_matrix_free(map1);
    gsl_matrix_free(at);

    return choose_class(p, this->ncls);

  }

  template <class real, class cls_t>
  cls_t multiclass<real, cls_t>::solve_class(gsl_vector *b, real *pdf) {
    cls_t cls1, cls2;
    real pt;		//total of conditional probability estimates
    real tally[this->ncls];

    cls2=-1;

    //printf("classification type=%d\n", type);
    switch(type) {
      case (-1):
        cls1=classify_map(b, pdf);
        break;
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
      default:
        cls1=classify_special(b, pdf);
        break;
    }

    if (cls2 >= 0) {        
      //renormalize the conditional prob.:
      pt=0;
      for (cls_t i=0; i<this->ncls; i++) {
        if (pdf[i]<0) pdf[i]=0;
        pt+=pdf[i];
      }
      //printf("pt=%g\n", pt);

      //correct the resultant conditional probabilities (=hack):
      for (int i=0; i<this->ncls; i++) pdf[i]=pdf[i]/pt;

      //if voting is different from matrix inversion, correct the results using a crude hack:
      if (cls1!=cls2) {
        real k;
        cls_t swp=cls1;
        cls1=cls2;
        cls2=swp;
        pt=pdf[cls2];
        for (cls_t i=0; i<this->ncls; i++) if (i!=cls1) pt+=pdf[i];
        k=(pdf[cls2]-pdf[cls1])/pt;
        for (cls_t i=0; i<this->ncls; i++) if (i!=cls1) pdf[i]=(1-k)*pdf[i];
        pdf[cls1]=pdf[cls2];
        pt=0;
        for (cls_t i=0; i<this->ncls; i++) pt+=pdf[i];
        //printf("pt (2)=%g\n", pt);
      }
    }

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
    cls_t nchild;
    cls_t D1, D2;

    if (this->D<=0) {
      D1=twoclass[0]->n_feat();
      //printf("multiclass: partition %d has %d features\n", 0, D1);
      for (cls_t i=1; i<nmodel; i++) {
        D2=twoclass[i]->n_feat();
        //printf("multiclass: partition %d has %d features\n", i, D2);
        if (D2!=D1) {
          fprintf(stderr, "multiclass: number of features in classifier %d does not match that in child %d", 0, i);
          fprintf(stderr, "                 %d vs. %d\n", D1, D2);
          exit(DIMENSION_MISMATCH);
        }
      }
      this->D=D1;
    }

    return this->D;
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

    if (imap!=NULL) {
      //screw it, lets just write the map directly into the control file:
      fprintf(fs, "\n[");
      print_matrix(fs, imap, this->ncls, nmodel);
      fprintf(fs, "]");
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

  template class multiclass<real_a, cls_ta>;

}

