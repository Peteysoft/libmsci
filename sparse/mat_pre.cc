#include "mat_pre.h"

#include "error_codes.h"

#include "full_matrix.h"
#include "full_util.h"
#include "sparse_array.h"

//initialization routines:
template <class index_t, class scalar>
mat_pre<index_t, scalar>::mat_pre() {
  m=0;
  n=0;
}

//initialization routines:
template <class index_t, class scalar>
mat_pre<index_t, scalar>::mat_pre(index_t m1, index_t n1) {
  m=m1;
  n=n1;
}

template <class index_t, class scalar>
mat_pre<index_t, scalar>::mat_pre(matrix_base<index_t, scalar> *other) {
  dum->dimensions(m, n);
};

//destructor:
template <class index_t, class scalar>
mat_pre<index_t, scalar>::~mat_pre() {
}

template <class index_t, class scalar>
scalar mat_pre<index_t, scalar>::operator() (index_t i, index_t j) {
  return 0;
}

template <class index_t, class scalar>
scalar *mat_pre<index_t, scalar>::operator() (index_t i) {
  scalar *row;
  row=new scalar[n];
  return row;
}

template <class index_t, class scalar>
long mat_pre<index_t, scalar>::cel (scalar val, index_t i, index_t j) {
  return i*m+j;
}

//multiply with full matrix:
template <class index_t, class scalar>
matrix_base<index_t, scalar> * mat_pre<index_t, scalar>::mat_mult(full_matrix<index_t, scalar> *cand) {
  mat_pre<index_t, scalar> *result;

  if (n!=cand->m) {
    fprintf(stderr, "Matrix inner dimensions must match ([%dx%d]*[%dx%d]\n", 
		    m, n, cand->m, cand->n);
    return NULL;
  }
  result=new mat_pre<index_t, scalar>(m, cand->n);
  return result;
} 

//multiply with sparse matrix:
template <class index_t, class scalar>
matrix_base<index_t, scalar> * mat_pre<index_t, scalar>::mat_mult(sparse<index_t, scalar> *cand) {
  mat_pre<index_t, scalar> *result;
  index_t m1, n1;

  cand->dimensions(m1, n1);

  if (n!=m1) {
    fprintf(stderr, "Matrix inner dimensions must match ([%dx%d]*[%dx%d]\n", 
		    m, n, m1, n1);
    return NULL;
  }
  result=new mat_pre<index_t, scalar>(m, n1);
  return result;
} 

//multiply with sparse array:
template <class index_t, class scalar>
matrix_base<index_t, scalar> * mat_pre<index_t, scalar>::mat_mult(sparse_array<index_t, scalar> *cand) {
  index_t m1, n1;
  mat_pre<index_t, scalar> *result;
  cand->dimensions(m1, n1);

  if (n!=m1) {
    fprintf(stderr, "Matrix inner dimensions must match ([%dx%d]*[%dx%d]\n", 
  		    m, n, m1, n1);
    return NULL;
  }
  result=new mat_pre<index_t, scalar>(m, n1);
  return result;
} 

//multiply with sparse array:
template <class index_t, class scalar>
matrix_base<index_t, scalar> * mat_pre<index_t, scalar>::mat_mult(mat_pre<index_t, scalar> *cand) {
  index_t m1, n1;
  mat_pre<index_t, scalar> *result;
  cand->dimensions(m1, n1);

  if (n!=m1) {
    fprintf(stderr, "Matrix inner dimensions must match ([%dx%d]*[%dx%d]\n", 
  		    m, n, m1, n1);
    return NULL;
  }
  result=new mat_pre<index_t, scalar>(m, n1);
  return result;
} 

//matrix addition:
template <class index_t, class scalar>
matrix_base<index_t, scalar> * mat_pre<index_t, scalar>::add(full_matrix<index_t, scalar> *b) {
  mat_pre<index_t, scalar> *result;
  index_t m1, n1;
  b->dimensions(m1, n1);
  if (n!=n1 || m1!=m) {
    fprintf(stderr, "Matrix dimensions must match ([%dx%d]+[%dx%d]\n", 
		    m, n, b->m, b->n);
    return NULL;
  }
  result=new mat_pre<index_t, scalar>();
  return result;
}

template <class index_t, class scalar>
matrix_base<index_t, scalar> * mat_pre<index_t, scalar>::add(sparse<index_t, scalar> *b) {
  mat_pre<index_t, scalar> *result;
  index_t m1, n1;
  b->dimensions(m1, n1);
  if (n!=n1 || m1!=m) {
    fprintf(stderr, "Matrix dimensions must match ([%dx%d]+[%dx%d]\n", 
		    m, n, b->m, b->n);
    return NULL;
  }
  result=new mat_pre<index_t, scalar>();
  return result;
}

template <class index_t, class scalar>
matrix_base<index_t, scalar> * mat_pre<index_t, scalar>::add(sparse_array<index_t, scalar> *b) {
  mat_pre<index_t, scalar> *result;
  index_t m1, n1;
  b->dimensions(m1, n1);
  if (n!=n1 || m1!=m) {
    fprintf(stderr, "Matrix dimensions must match ([%dx%d]+[%dx%d]\n", 
		    m, n, b->m, b->n);
    return NULL;
  }
  result=new mat_pre<index_t, scalar>();
  return result;
}

//multiply with a vector:
template <class index_t, class scalar>
scalar * mat_pre<index_t, scalar>::vect_mult(scalar *cand) {
  scalar *result;
  result=new scalar[n];
  return result;
} 

template <class index_t, class scalar>
scalar * mat_pre<index_t, scalar>::left_mult(scalar *cor) {
  scalar *result;
  result=new scalar[m];
  return result;
} 

//multiply with a vector:
template <class index_t, class scalar>
void mat_pre<index_t, scalar>::vect_mult(scalar *cand, scalar *result) {
} 

template <class index_t, class scalar>
void mat_pre<index_t, scalar>::left_mult(scalar *cor, scalar *result) {
} 

//copy constructor:
template <class index_t, class scalar>
mat_pre<index_t, scalar>::mat_pre(mat_pre<index_t, scalar> &other) {
  m=other.m;
  n=other.n;
}

//type conversion:
template <class index_t, class scalar>
matrix_base<index_t, scalar> & mat_pre<index_t, scalar>::operator = 
		(full_matrix<index_t, scalar> &other) {
  m=other.m;
  n=other.n;
  return *this;
}

//copy constructor
template <class index_t, class scalar>
mat_pre<index_t, scalar>::mat_pre(sparse<index_t, scalar> &other) {
  other.dimensions(m, n);
}

//type conversion:
template <class index_t, class scalar>
matrix_base<index_t, scalar> & mat_pre<index_t, scalar>::operator = 
		(sparse<index_t, scalar> &other) {
  other.dimensions(m, n);
  return *this;
}

template <class index_t, class scalar>
mat_pre<index_t, scalar>::mat_pre(sparse_array<index_t, scalar> &other) {
  other.dimensions(m, n);
}

template <class index_t, class scalar>
matrix_base<index_t, scalar> & mat_pre<index_t, scalar>::operator = 
		(sparse_array<index_t, scalar> &other) {
  other.dimensions(m, n);
  return *this;
}

template <class index_t, class scalar>
void mat_pre<index_t, scalar>::dimensions(index_t &mout, index_t &nout) const {
  mout=m;
  nout=n;
}

//should probably return a new matrix:
template <class index_t, class scalar>
void mat_pre<index_t, scalar>::transpose() {
}

//should probably return a new matrix:
template <class index_t, class scalar>
void mat_pre<index_t, scalar>::scal_mult(scalar cand) {
}
	
template <class index_t, class scalar>
int mat_pre<index_t, scalar>::scan(FILE *fs) {
  char line[1000];

  if (fgets(line, MAXLL, fptr)==NULL) return FILE_READ_ERROR;
  sscanf(line, "%d %d\n", &m, &n);
 
  return 0;
}

template <class index_t, class scalar>
size_t mat_pre<index_t, scalar>::read(FILE *fs, int mattype) {
  size_t nread=0;
  if (mattype==FULL_T) {
    nread=fread(&n, sizeof(m), 1, fs);
    fseek(fptr, 0, SEEK_END);
    m=(ftell(fptr)-sizeof(n))/sizeof(real)/n;
    fseek(fptr, sizeof(n), SEEK_SET);
  } else if (mattype==SPARSE_ARRAY_T || mattype==SPARSE_T) {
    nread=fread(&m, sizeof(m), 1, fptr);
    nread+=fread(&n, sizeof(n), 1, fptr);
  } else {
    m=0; n=0;
  }

  return nread;
}

template <class index_t, class scalar>
size_t mat_pre<index_t, scalar>::write(FILE *fs) {
  return 0;
}

template <class index_t, class scalar>
void mat_pre<index_t, scalar>::print(FILE *fs) {
  fprintf("%d %d\n", m, n);
}

template <class index_t, class scalar>
matrix_base<index_t, scalar> * mat_pre<index_t, scalar>::clone() {
  mat_pre<index_t, scalar> *result;
  result=new mat_pre<index_t, scalar>(m, n);
  return result;

}

template <class index_t, class scalar>
scalar full_matrix<index_t, scalar>::norm() {
  return 0;
}

//template class full_matrix<int16_t, float>;
template class mat_pre<ind_t, float>;
template class mat_pre<ind_t, double>;

