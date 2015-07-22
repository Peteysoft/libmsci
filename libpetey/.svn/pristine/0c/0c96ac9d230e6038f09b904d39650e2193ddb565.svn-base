#include <string.h>

#include "error_codes.h"

#include "sc_type.h"

using namespace libpetey;
using namespace sparse;

namespace sparse_calc {

//default is, it's illegal:
sc_t * sc_t::mult(sc_t *v) {
  return NULL;
}

sc_t * sc_t::add(sc_t *v) {
  return NULL;
}

sc_t * sc_t::sub(sc_t *v) {
  return NULL;
}

sc_t * sc_t::sub(sc_t *v1, sc_t *v2) {
  return NULL;
}

sc_t * sc_t::cprod(sc_t *v) {
  return NULL;
}

sc_t * sc_t::norm() {
  return NULL;
}

//shouldn't happen:
void sc_t::neg() {
  assert(0);
}

//shouldn't happen:
int sc_t::load(FILE *fs) {
  assert(0);
}

//shouldn't happen:
int sc_t::save(FILE *fs) {
  assert(0);
}

//shouldn't happen:
int sc_t::read(FILE *fs) {
  assert(0);
}

//shouldn't happen:
int sc_t::print(FILE *fs) {
  assert(0);
}

//*********************  string type  *************************
//
sc_str_t::sc_str_t(char *str) {
  s=new char[strlen(str)+1];
  strcpy(s, str);
}

//*********************  scalar type  *************************
//
template <class scalar_t>
sc_scal_t<scalar_t>::sc_scal_t(scalar_t val) {
  value=val;
}

template <class scalar_t>
sc_t * sc_scal_t<scalar_t>::mult(sc_t *cand) {
  sc_t *result;
  int t=sc_type_of(cand);

  switch (t) {
    //multiply with another scalar:
    case (SC_SCALAR_T):
      result=new sc_scal_t<scalar_t>(((sc_scal_t<scalar_t> *) cand)->value*value);
      break;
    //multiply with vector:
    case(SC_VECTOR_T):
    //multiply with matrix:
    case(SC_FULL_T):
    case(SC_SPARSE_T):
    case(SC_SPARSE_ARRAY_T):
    //multiply with list type:
    case(SC_LIST_T):
      result=cand->mult(this);
      break;
    default:
      result=NULL;
  }

  return result;
}

template <class scalar_t>
sc_t * sc_scal_t<scalar_t>::add(sc_t *cand) {
  sc_t *result;
  int t=sc_type_of(cand);

  switch (t) {
    //add to another scalar:
    case (SC_SCALAR_T):
      result=new sc_scal_t<scalar_t>(((sc_scal_t<scalar_t> *) cand)->value+value);
      break;
    //add to a vector:
    case(SC_VECTOR_T):
    //add to a matrix:
    case(SC_FULL_T):
    case(SC_SPARSE_T):
    case(SC_SPARSE_ARRAY_T):
    //add to list type:
    case(SC_LIST_T):
      result=cand->add(this);
      break;
    default:
      result=NULL;
  }

  return result;
}

template <class scalar_t>
void sc_scal_t<scalar_t>::neg() {
  value=-value;
}

template <class scalar_t>
void sc_scal_t<scalar_t>::norm() {
  return new sc_scal_t<scalar_t>(fabs(value));
}

template <class scalar_t>
void sc_scal_t<scalar_t>::sub(sc_t *sub) {
  return NULL; 
}

//*********************  vector type  *************************

template <class scalar_t>
sc_vec_t<scalar_t>::sc_vec_t() {
  n=0;
  data=NULL;
}

template <class scalar_t>
sc_vec_t<scalar_t>::~sc_vec_t() {
  delete [] data;
}

template <class scalar_t>
sc_vec_t<scalar_t>::sc_vec_t(scalar_t start, scalar_t finish) {
  n=(finish-start);
  data=new scalar_t[n];
  for (scalar_t i=0; i<n; i++) data[i]=start+i;
}

template <class scalar_t>
sc_vec_t<scalar_t>::sc_vec_t(sc_int_t n1) {
  n=n1;
  data=new scalar_t[n];
}

template <class scalar_t>
int sc_vec_t<scalar_t>::load(char *fname) {
  FILE *fs=fopen(fname, "r");
  if (fs==NULL) {
    fprintf(stderr, "Error loading vector %s\n", fname);
    return UNABLE_TO_OPEN_FILE_FOR_READING;
  }
  fseek(fs, 0, SEEK_END);
  n=ftell(fs)/sizeof(scalar_t);
  fseek(fs, 0, SEEK_SET);
  data=new scalar_t[n];
  fread(data, n, sizeof(scalar_t), fs);
  fclose(fs);

  return 0;
}

template <class scalar_t>
int sc_vec_t<scalar_t>::save(char *fname) {
  int nwrit;
  FILE *fs=fopen(fname, "w");
  if (fs==NULL) {
    fprintf(stderr, "Error saving vector %s\n", fname);
    return UNABLE_TO_OPEN_FILE_FOR_WRITING;
  }
  nwrit=fwrite(data, n, sizeof(scalar_t), fs);
  if (nwrit!=n) {
    fprintf(stderr, "Error loading vector %s\n", fname);
    return FILE_WRITE_ERROR;
  }
  fclose(fs);
  return 0;
}

template <class scalar_t>
int sc_vec_t<scalar_t>::read(FILE *fs) {
  int err;
  scalar_t val;
  linked_list<scalar_t> list;
  long n1;

  while (feof(fs)==0) {
    err=fscanf(fs, "%g", val);
    if (err!=1) break;
    list.add(val);
  }

  data=list.make_array(n1);
  n=n1;
}

template <class scalar_t>
int sc_vec_t<scalar_t>::print(FILE *fs) {
  for (sc_int_t i=0; i<n; i++) fprintf(fs, "%g\n", data[i]);
}

template <class scalar_t>
sc_t * sc_vec_t<scalar_t>::mult(sc_t *cand) {
  sc_t *result;
  int t=sc_type_of(cand);

  switch (t) {
    //multiply with a scalar:
    case (SC_SCALAR_T):
      sc_vec_t<scalar_t> *r1;
      scalar_t value=((sc_scal_t<scalar_t> *) cand)->value;
      r1=new sc_vec_t<scalar_t>(n);
      for (int i=0; i<n; i++) r1->data[i]=data[i]*value;
      result=r1;
      break;
    //multiply with vector:
    case(SC_VECTOR_T):
      sc_vec_t<scalar_t> *c1;
      scalar_t r1;
      c1=(sc_vec_t<scalar_t> *) cand;
      if (c1->n != n) {
        fprintf(stderr, "Vector multiplication: vectors are different sizes (%d vs. %d)\n", n, c1->n);
        result=NULL;
        break;
      }
      r1=0;
      for (sc_int_t i=0; i<n; i++) r1+=data[i]*c1->data[i];
      result=new sc_scal_t<scalar_t>(r1);
      break;
    //multiply with matrix:
    case(SC_FULL_T):
    case(SC_SPARSE_T):
    case(SC_SPARSE_ARRAY_T):
      sc_int_t m1, n1;
      sc_vec_t<scalar_t> *r1;
      sc_mat_t<scalar_t> *c1=(sc_mat_t<scalar_t> *) cand;
      cand->mat->dimensions(m1, n1);
      if (m1 != n) {
        fprintf(stderr, "Vector-matrix product: vector dimension doesn't match matrix inner dimension (%d vs. %d)\n", n, m1);
        result=NULL;
        break;
      }
      r1=new sc_vec_t<scalar_t>(n1);
      c1->mat->left_mult(data, r1->data);
      result=r1;

    //multiply with list type:
    case(SC_LIST_T):
      break;
    default:
      result=NULL;
  }

  return result;
}

template <class scalar_t>
sc_t * sc_vec_t<scalar_t>::add(sc_t *cand) {
  sc_t *result;
  int t=sc_type_of(cand);

  switch (t) {
    //multiply with a scalar:
    case (SC_SCALAR_T):
      sc_vec_t<scalar_t> *r1;
      scalar_t value=((sc_scal_t<scalar_t> *) cand)->value;
      r1=new sc_vec_t<scalar_t>(n);
      for (int i=0; i<n; i++) r1->data[i]=data[i]+value;
      result=r1;
      break;
    //multiply with vector:
    case(SC_VECTOR_T):
      sc_vec_t<scalar_t> *c1;
      c1=(sc_vec_t<scalar_t> *) cand;
      if (c1->n != n) {
        fprintf(stderr, "Vector inner product: vectors are different sizes (%d vs. %d)\n", n, c1->n);
        result=NULL;
      } else {
        sc_vec_t<scalar_t> *r1;
        r1=new sc_vec_t<scalar_t>(n);
        for (sc_int_t i=0; i<n; i++) r1->data[i]=data[i]+c1->data[i];
        result=r1;
      }
      break;
    //multiply with matrix:
    case(SC_FULL_T):
    case(SC_SPARSE_T):
    case(SC_SPARSE_ARRAY_T):
    case(SC_LIST_T):
      cand->add(this);
      break;
    default:
      result=NULL;
  }

  return result;
}

template <class scalar_t>
void sc_vec_t<scalar_t>::neg() {
  for (sc_int_t i; i<n; i++) data[i]=-data[i];
}

template <class scalar_t>
sc_t * sc_vec_t<scalar_t>::sub (sc_t *sub) {
  integer ind;
  sc_t *result;
  int t=sc_type_of(sub);

  switch (t) {
    //subscript with scalar:
    case (SC_SCALAR_T):
      ind=((sc_scal_t<scalar_t> *) sub)->value;
      if (ind<0 || ind>=n) {
        fprintf(stderr, "Scalar subscript (%d) of vector out-of-range\n", ind);
        result=NULL;
        break;
      }
      result=new sc_scal_t<scalar_t>(data[ind]);
      break;
    //subscript with vector:
    case(SC_VECTOR_T):
      sc_vec_t<scalar_t> *s1=(sc_vec_t<scalar_t> *) sub;
      sc_vec_t<scalar_t> *r1=new sc_vec_t<scalar_t>(s1->n);
      for (sc_int_t i=0; i<s1->n; i++) {
        ind=s1->data[i];
        if (ind<0 || ind>=n) {
          fprintf(stderr, "Vector subscript (%d) of vector out-of-range\n", ind);
          result=NULL;
          delete [] r1;
          break;
        }
        if (r1==NULL) break;
        r1->data[i]=data[ind];
      }
    //type mismatch:
    case(SC_FULL_T):
    case(SC_SPARSE_T):
    case(SC_SPARSE_ARRAY_T):
    case(SC_LIST_T):
    default:
      fprintf(stderr, "Type mismatch: only scalar or vector subscripts allowed\n");
      result=NULL;
  }

  return result;
}

template <class scalar_t>
sc_t * sc_vec_t<scalar_t>::cprod(sc_t *cand) {
  integer m, n;
  sc_mat_t<scalar_t> *result=NULL;
  int t=sc_type_of(cand);

  if (t==SC_SPARSE_ARRAY_T) {
    sparse_array_t *c1=(sparse_array_t *) mat;
    sc_vec_t<scalar_t> *cor=(sc_vec_t<scalar_t> *) cand;
    cor->dimensions(m, n);
    if (n != c1->m) {
      fprintf(stderr, "Cumulative multiplication: inner dimensions must match (%d vs. %d)\n", c1->n, m);
      return NULL;
    }
    result=new sc_mat_t<scalar_t>();
    result->mat=mat->cmult(cor->data);
  } else {
    fprintf(stderr, "Type mismatch in cumulative product operator (#).");
  }

  return result;
}

template <class scalar_t>
sc_t * sc_vec_t<scalar_t>::norm() {
  sc_scal_t<scalar_t> *r1=new sc_scal_t<scalar_t>(0);
  for (integer i=0; i<n; i++) r1->val+=data[i]*data[i];
  return r1;
}

//*********************  matrix type  *************************

template <class scalar_t>
sc_t * sc_mat_t<scalar_t>::sc_mat_t() {
  mat=NULL;
}

template <class scalar_t>
sc_t * sc_mat_t<scalar_t>::sc_mat_t(int type) {
  switch (type) {
    case(SPARSE_T):
      mat=new sparse<scalar_t, sc_int_t>();
      break;
    case(SPARSE_ARRAY_T):
      mat=new sparse_array<scalar_t, sc_int_t>();
      break;
    case(FULL_T):
      mat=new full_matrix<scalar_t, sc_int_t>();
      break;
  }
}


template <class scalar_t>
sc_t * sc_mat_t<scalar_t>::~sc_mat_t() {
  delete mat;
}

template <class scalar_t>
sc_t * sc_mat_t<scalar_t>::load(char *fname) {
  FILE *fs=fopen(fname, "r");
  if (fs==NULL) {
    fprintf(stderr, "Error loading matrix %s\n", fname);
    return UNABLE_TO_OPEN_FILE_FOR_READING;
  }
  mat->read(fs);
  fclose(fs);
  return 0;
}

template <class scalar_t>
sc_t * sc_mat_t<scalar_t>::save(char *fname) {
  FILE *fs=fopen(fname, "w");
  if (fs==NULL) {
    fprintf(stderr, "Error saving matrix %s\n", fname);
    return UNABLE_TO_OPEN_FILE_FOR_WRITING;
  }
  mat->write(fs);
  fclose(fs);
  return 0;
}

template <class scalar_t>
sc_t * sc_mat_t<scalar_t>::read(FILE *fs) {
  mat->scan(fs);
}

template <class scalar_t>
sc_t * sc_mat_t<scalar_t>::print(FILE *fs) {
  mat->print(fs);
}


template <class scalar_t>
sc_t * sc_mat_t<scalar_t>::mult(sc_t *cand) {
  sc_t *result;
  int t=sc_type_of(cand);

  switch (t) {
    //multiply with a scalar:
    case (SC_SCALAR_T):
      sc_mat_t<scalar_t> *r1=new sc_mat_t<scalar_t>();
      r1->mat=copy_matrix(this);
      r1->mat->scal_mult(((sc_scal_t<scalar_t> *) cand)->value);
      result=r1;
      break;
    //multiply with a vector:
    case(SC_VECTOR_T):
      sc_vect_t<scalar_t> * c1=(sc_vect_t<scalar_t> *) cand;
      sc_int_t m, n;
      mat->dimensions(m, n);
      if (c1->n != n) {
        fprintf(stderr, "Vector multiplication: inner dimensions must match (%d vs. %d)\n", n, c1->n);
        result=NULL;
      } else {
        r1->data=mat->vect_mult(c1->data);
        r1->n=m;
      }
      result=r1;
      break;
    //multiply with matrix:
    case(SC_FULL_T):
    case(SC_SPARSE_T):
    case(SC_SPARSE_ARRAY_T):
      sc_mat_t<scalar_t> *c1=(sc_mat_t<scalar_t> *) cand;
      sc_mat_t<scalar_t> *r1=new sc_mat_t<scalar_t>();
      r1->mat=mat->mult(c1->mat);
      result=r1;
      break;
    //multiply with list type:
    case(SC_LIST_T):
      sc_list_t *c1=(sc_list_t<scalar_t> *) cand;
      result=c1->left_mult(this);
      break;
    default:
      result=NULL;
  }

  return result;
}

template <class scalar_t>
sc_t * sc_mat_t<scalar_t>::add(sc_t *cand) {
  sc_mat_t<scalar_t> *result;
  sc_int_t m, n;

  //gather some information:
  int t=sc_type_of(cand);
  mat->dimensions(m, n);

  //result is always a matrix:
  result=new sc_mat_t();

  switch (t) {
    //add to a scalar:
    case (SC_SCALAR_T):
      full_matrix<sc_int_t, scalar_t> *c1=new full_matrix<sc_int_t, scalar_t>(m, n);
      c1->ones();
      c1->scal_mult(((sc_scal_t<scalar_t> *) cand)->value);
      result->mat=mat->add(c1);
      break;
    //add to a vector:
    case(SC_VECTOR_T):
      sc_vect_t<scalar_t> * c1=(sc_vect_t<scalar_t> *) cand;
      sparse<sc_int_t, scalar_t> *c2=new sparse<sc_int_t, scalar_t>(m, n);
      c2->identity();
      for (sc_int_t i=0; i<m && i<n && i<c1->n; i++) {
        c2->cel(c1->data[i], i, i);
      }
      result->data=mat->add(c2);
      break;
    //add to matrix:
    case(SC_FULL_T):
    case(SC_SPARSE_T):
    case(SC_SPARSE_ARRAY_T):
      sc_mat_t<scalar_t> *c1=(sc_mat_t<scalar_t> *) cand;
      result->mat=mat->add(c1->mat);
      break;
    //add to list type:
    case(SC_LIST_T):
      cand->add(this);
      break;
    default:
      delete [] result;
      result=NULL;
  }

  return result;
}

template <class scalar_t>
sc_t * sc_mat_t<scalar_t>::sub (sc_t *sub) {
  integer m, n;
  sc_t *result;
  sc_int_t ind;
  int t=sc_type_of(sub);

  mat->dimensions(m, n);

  switch (t) {
    //extract a row:
    case (SC_SCALAR_T):
      sc_vec_t<scalar_t> *r1;
      ind=((sc_scal_t<scalar_t> *) sub)->value;
      if (ind<0 || ind>=sub->m) {
        fprintf(stderr, "Scalar subscript (%d) of matrix out-of-range\n", ind);
        result=NULL;
        break;
      }
      r1=new sc_vec_t<scalar_t>(n);
      mat->get_row(ind, r1->data);
      result=r1;
      break;
    //subscript with vector:
    case(SC_VECTOR_T):
      fprintf(stderr, "Vector subscripting of rows not yet supported.  Sorry\n");
      result=NULL;
    //type mismatch:
    case(SC_FULL_T):
    case(SC_SPARSE_T):
    case(SC_SPARSE_ARRAY_T):
    case(SC_LIST_T):
    default:
      fprintf(stderr, "Type mismatch: only scalar or vector subscripts allowed\n");
      result=NULL;
  }

  return result;
}

template <class scalar_t>
sc_t * sc_mat_t<scalar_t>::sub (sc_t *sub1, sc_t *sub2) {
  integer m, n;
  sc_t *result;
  sc_int_t ind;
  int t1=sc_type_of(sub1);
  int t2=sc_type_of(sub2);

  mat->dimensions(m, n);

  switch (t) {
    //extract a row:
    case (SC_SCALAR_T):
      sc_vec_t<scalar_t> *r1;
      ind=((sc_scal_t<scalar_t> *) sub)->value;
      if (ind<0 || ind>=sub->m) {
        fprintf(stderr, "Scalar subscript (%d) of matrix out-of-range\n", ind);
        result=NULL;
        break;
      }
      r1=new sc_vec_t<scalar_t>(n);
      mat->get_row(ind, r1->data);
      result=r1;
      break;
    //subscript with vector:
    case(SC_VECTOR_T):
      fprintf(stderr, "Vector subscripting of rows not yet supported.  Sorry\n");
      result=NULL;
    //type mismatch:
    case(SC_FULL_T):
    case(SC_SPARSE_T):
    case(SC_SPARSE_ARRAY_T):
    case(SC_LIST_T):
    default:
      fprintf(stderr, "Type mismatch: only scalar or vector subscripts allowed\n");
      result=NULL;
  }

  return result;
}

template <class scalar_t>
sc_t * sc_mat_t<scalar_t>::cprod(sc_t *cand) {
  sc_int_t m, n;
  sc_mat_t<scalar_t> *result=NULL;
  int t=sc_type_of(cand);

  if (typeid(*mat)!=typeid(sparse_array_t)) {
    fprintf(stderr, "Type mismatch");
    return NULL;
  }

  if (t==SC_VEC_T) {
    sparse_array_t *cor=(sparse_array_t *) mat;
    sc_vec_t<scalar_t> *c1=(sc_vec_t<scalar_t> *) cand;
    cor->dimensions(m, n);
    if (n != c1->n) {
      fprintf(stderr, "Cumulative multiplication: inner dimensions must match (%d vs. %d)\n", n, c1->n);
      return NULL;
    }
    result=new sc_mat_t<scalar_t>();
    result->mat=cor->cmult(c1->data);
  } else {
    //need to define for matrix (full and sparse) multiplicands...
    fprintf(stderr, "Type mismatch");
  }

  return result;
}

template <class scalar_t>
sc_t * sc_mat_t<scalar_t>::norm() {
  sc_scal_t<scalar_t> *r1=new sc_scal_t<scalar_t>(mat->norm());
  return r1;
}
    
//*********************  list type  *************************

template <typename scalar_t>
sc_t * sc_list_t<scalar_t>::sc_list_t() {
  list=NULL;
}

template <typename scalar_t>
sc_t * sc_list_t<scalar_t>::sc_list_t(int n) {
  list=new vector_e<sc_t *>(n);
}

template <typename scalar_t>
sc_list_t<scalar_t>::~sc_list_t() {
  int n=list->size();
  for (int i=0; i<n; i++) delete list[i];
}

template <typename scalar_t>
int sc_list_t<scalar_t>::load(char *fname) {
  int n=0;
  for (int i=0; i<n; i++) delete list[i];
  
}

template <typename scalar_t>
sc_t * sc_list_t<scalar_t>::mult(sc_t *cand) {
  sc_int_t n=list->size();
  sc_t *cur, *mnew;

  if (n==0) {
    fprintf(stderr, "Empty list\n");
    return NULL;
  }

  cur=cand->mult(list[n-1]);
  for (sc_int_t i=n-2; i>=0; i--) {
    if (cur==NULL) return NULL;
    mnew=cur->mult(list[i]);
    delete cur;
    cur=mnew;
  }

  return cur;

}

template <typename scalar_t>
sc_t * sc_list_t<scalar_t>::left_mult(sc_t *cor) {
  sc_int_t n=list->size();
  sc_t *cur, *mnew;

  if (n==0) {
    fprintf(stderr, "Empty list\n");
    return NULL;
  }

  cur=list[0]->mult(cor);
  for (sc_int_t i=1; i<n; i++) {
    if (cur==NULL) return NULL;
    mnew=list[i]->mult(cur);
    delete cur;
    cur=mnew;
  }

  return cur;

}

template <typename scalar_t>
sc_t * sc_list_t<scalar_t>::resolve() {
  sc_int_t n=list->size();
  sc_t *cur, *mnew;

  if (n==0) {
    fprintf(stderr, "Empty list\n");
    return NULL;
  }

  cur=list[0];
  for (sc_int_t i=1; i<n; i++) {
    if (cur==NULL) return NULL;
    mnew=list[i]->mult(cur);
    delete cur;
    cur=mnew;
  }

  return cur;

}

template <typename scalar_t>
sc_t * sc_list_t<scalar_t>::add(sc_t *b) {
  sc_int_t n=list->size();
  sc_t *cur, *mnew;

  if (n==0) {
    fprintf(stderr, "Empty list\n");
    return NULL;
  }

  cur=list[0];
  for (sc_int_t i=1; i<n; i++) {
    mnew=cur->mult(list[i]);
    delete cur;
    if (mnew==NULL) return NULL;
    cur=mnew;
  }

  cur->add(cand);

}

template <typename scalar_t>
sc_t * sc_list_t<scalar_t>::sub(sc_t *sub) {
  sc_t *res=resolve();
  sc_t *result=res->sub(sub);
  delete res;
  return result;
}

template <class scalar_t>
sc_t * sc_mat_t<scalar_t>::norm() {
  sc_t *res=resolve();
  sc_scal_t<scalar_t> *r1=new sc_scal_t<scalar_t>(res->norm());
  return r1;
}
    
