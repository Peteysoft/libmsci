#include "error_codes.h"

#include "full_matrix.h"
#include "full_util.h"
#include "sparse_array.h"

namespace libpetey {
namespace libsparse {

//initialization routines:
template <class index_t, class scalar>
full_matrix<index_t, scalar>::full_matrix() {
  m=0;
  n=0;
  data=NULL;
}

template <class index_t, class scalar>
full_matrix<index_t, scalar>::full_matrix(matrix_base<index_t, scalar> *other) {
  if (typeid(*other) == typeid(sparse<index_t, scalar>)) {
    sparse<index_t, scalar> * dum;
    //printf("Assigning sparse matrix to full matrix\n");
    dum=(sparse<index_t, scalar> *) other;
    dum->dimensions(m, n);
    data=(scalar **) (*dum);
    //print_matrix(stdout, data, m, n);
  } else if (typeid(*other) == typeid(sparse_array<index_t, scalar>)) {
    sparse_array<index_t, scalar> * dum;
    printf("Converting sparse array to full:");
    other->dimensions(m, n);
    dum=(sparse_array<index_t, scalar>*) other;
    data=(scalar **) (*dum);
  } else if (typeid(*other) == typeid(full_matrix<index_t, scalar>)) {
    full_matrix<index_t, scalar> * dum;
    dum=(full_matrix<index_t, scalar> *) other;
    m=dum->m;
    n=dum->n;
    data=libpetey::copy_matrix<scalar, index_t>(dum->data, m, n);
  } else {
    fprintf(stderr, "Failed to initialize full matrix from base class\n");
  }

};

/*
template <class index_t, class scalar>
full_matrix<index_t, scalar> & full_matrix<index_t, scalar>::operator = (matrix_base<index_t, scalar> &other) {
  if (typeid(other) == typeid(sparse<index_t, scalar>)) {
    return *this=*((sparse<index_t, scalar> *) &other);
  } else if (typeid(other) == typeid(sparse_array<index_t, scalar>)) {
    return *this=*(sparse_array<index_t, scalar> *) &other;
  } else if (typeid(other) == typeid(full_matrix<index_t, scalar>)) {
    return *this=*(full_matrix<index_t, scalar> *) &other;
  }

};
*/

//initialization routines:
template <class index_t, class scalar>
full_matrix<index_t, scalar>::full_matrix(index_t min, index_t nin) {
  m=min;
  n=nin;
  data=zero_matrix<scalar, index_t>(m, n);
}

template <class index_t, class scalar>
void full_matrix<index_t, scalar>::ones() {
  for (index_t i=0; i<m; i++) {
    for (index_t j=0; j<n; j++) {
      data[i][j]=1;
    }
  }
}

template <class index_t, class scalar>
void full_matrix<index_t, scalar>::identity() {
  zero_matrix<scalar, index_t>(data, m, n);
  for (index_t i=0; i<m && i<n; i++) data[i][i]=1;
}

template <class index_t, class scalar>
full_matrix<index_t, scalar>::full_matrix(scalar *dat, index_t min, index_t nin) {
  m=min;
  n=nin;
  data=allocate_matrix<scalar, index_t>(m, n);
  for (index_t i=0; i<n*m && i<m; i++) data[0][i]=dat[i];
}

template <class index_t, class scalar>
full_matrix<index_t, scalar>::full_matrix(scalar **dat, index_t min, index_t nin) {
  m=min;
  n=nin;
  data=libpetey::copy_matrix<scalar, index_t>(dat, m, n);
}

//destructor:
template <class index_t, class scalar>
full_matrix<index_t, scalar>::~full_matrix() {
  delete_matrix(data);
}

template <class index_t, class scalar>
scalar full_matrix<index_t, scalar>::operator() (index_t i, index_t j) {
  return data[i][j];
}

template <class index_t, class scalar>
void full_matrix<index_t, scalar>::get_row(index_t i, scalar *row) {
  //pretty simple:
  for (index_t j=0; j<n; j++) row[j]=data[i][j];
}

template <class index_t, class scalar>
scalar *full_matrix<index_t, scalar>::operator() (index_t i) {
  scalar *row;
  row=new scalar[n];
  get_row(i, row);
  return row;
}

template <class index_t, class scalar>
long full_matrix<index_t, scalar>::cel (scalar val, index_t i, index_t j) {
  data[i][j]=val;
  return i*m+j;
}

//multiply with another full matrix:
template <class index_t, class scalar>
matrix_base<index_t, scalar> * full_matrix<index_t, scalar>::mat_mult(matrix_base<index_t, scalar> *cand) {
  matrix_base<index_t, scalar> *result;
  index_t m1, n1;
  cand->dimensions(m1, n1);
  if (n!=m1) {
    fprintf(stderr, "Matrix inner dimensions must match ([%dx%d]*[%dx%d]\n", 
		    m, n, m1, n1);
    result=NULL;
  } else if (typeid(*cand) == typeid(full_matrix<index_t, scalar>)) {
    result=new full_matrix<index_t, scalar>(m, n1);
    //ugly way of doing it:
    matrix_mult<scalar, index_t>(data, 
		((full_matrix<index_t, scalar> *) cand)->data, 
		((full_matrix<index_t, scalar> *) result)->data, 
		m, n, n1);
  } else if (typeid(*cand) == typeid(sparse<index_t, scalar>)) {
    result=new full_matrix<index_t, scalar>(m, n1);
    ((sparse<index_t, scalar> *) cand)->left_m_mult(this->data, 
		((full_matrix<index_t, scalar> *) result)->data, m);
  } else if (typeid(*cand) == typeid(sparse_array<index_t, scalar>)) {
    if (m != n) {
      full_matrix<index_t, scalar> cand2;
      cand2=*(full_matrix<index_t, scalar> *) cand;
      result=mat_mult(&cand2);
    } else {
      result=new full_matrix<index_t, scalar>(m, n1);
      ((sparse_array<index_t, scalar> *) cand)->left_m_mult(this->data, 
			((full_matrix<index_t, scalar> *) result)->data);
    }
  } else if (typeid(*cand) == typeid(matrix_base<index_t, scalar>)) {
    result=new matrix_base<index_t, scalar>(m, n1);
    //do nothing
  } else {
    fprintf(stderr, "full_matrix::mat_mult - error, matrix type not recognized\n");
    result=NULL;
  } 
  return result;
} 

//matrix addition:
template <class index_t, class scalar>
matrix_base<index_t, scalar> * full_matrix<index_t, scalar>::add(matrix_base<index_t, scalar> *b) {
  matrix_base<index_t, scalar> *result;
  index_t m1, n1;
  b->dimensions(m1, n1);

  if (n!=n1 || m1!=m) {
    fprintf(stderr, "Matrix dimensions must match ([%dx%d]+[%dx%d]\n", 
		    m, n, m1, n1);
    result=NULL;
  } else if (typeid(*b) == typeid(full_matrix<index_t, scalar>)) {
    result=new full_matrix<index_t, scalar>();
    ((full_matrix<index_t, scalar> *) result)->data =
		libpetey::copy_matrix<scalar, index_t>(data, m, n);
    matrix_add<scalar, index_t>(
		((full_matrix<index_t, scalar> *) result)->data, 
		((full_matrix<index_t, scalar> *) b)->data, m, n);
    ((full_matrix<index_t, scalar> *) result)->m=m;
    ((full_matrix<index_t, scalar> *) result)->n=n;
  } else if (typeid(*b) == typeid(sparse<index_t, scalar>)) {
    result=b->add(this);
  } else if (typeid(*b) == typeid(sparse_array<index_t, scalar>)) {
    result=b->add(this);
  } else if (typeid(*b) == typeid(matrix_base<index_t, scalar>)) {
    result=new matrix_base<index_t, scalar>(m, n1);
    //do nothing
  } else {
    fprintf(stderr, "full_matrix::add - error, matrix type not recognized\n");
    result=NULL;
  } 
  return result;
}

//multiply with a vector:
template <class index_t, class scalar>
scalar * full_matrix<index_t, scalar>::vect_mult(scalar *cand) {
  return vector_mult<scalar, index_t>(data, cand, m, n);
} 

template <class index_t, class scalar>
scalar * full_matrix<index_t, scalar>::left_mult(scalar *cor) {
  return left_vec_mult<scalar, index_t>(cor, data, m, n);
} 

//multiply with a vector:
template <class index_t, class scalar>
void full_matrix<index_t, scalar>::vect_mult(scalar *cand, scalar *result) {
  vector_mult<scalar, index_t>(data, cand, result, m, n);
} 

template <class index_t, class scalar>
void full_matrix<index_t, scalar>::left_mult(scalar *cor, scalar *result) {
  left_vec_mult<scalar, index_t>(cor, data, result, m, n);
} 

//copy constructor:
template <class index_t, class scalar>
full_matrix<index_t, scalar>::full_matrix(full_matrix<index_t, scalar> &other) {
  m=other.m;
  n=other.n;
  data=libpetey::copy_matrix<scalar, index_t>(other.data, m, n);
}

//type conversion:
template <class index_t, class scalar>
matrix_base<index_t, scalar> & full_matrix<index_t, scalar>::operator = 
		(full_matrix<index_t, scalar> &other) {
  if (data != NULL) delete_matrix(data);
  m=other.m;
  n=other.n;
  data=libpetey::copy_matrix<scalar, index_t>(other.data, m, n);
  return *this;
}

//copy constructor
template <class index_t, class scalar>
full_matrix<index_t, scalar>::full_matrix(sparse<index_t, scalar> &other) {
  other.dimensions(m, n);
  data=(scalar **) other;
}

//type conversion:
template <class index_t, class scalar>
matrix_base<index_t, scalar> & full_matrix<index_t, scalar>::operator = 
		(sparse<index_t, scalar> &other) {
  if (data != NULL) delete_matrix(data);
  other.dimensions(m, n);
  data=(scalar **) other;
  return *this;
}

template <class index_t, class scalar>
full_matrix<index_t, scalar>::full_matrix(sparse_array<index_t, scalar> &other) {
  other.dimensions(m, n);
  data=(scalar **) other;
}

template <class index_t, class scalar>
matrix_base<index_t, scalar> & full_matrix<index_t, scalar>::operator = 
		(sparse_array<index_t, scalar> &other) {
  if (data != NULL) delete_matrix(data);
  other.dimensions(m, n);
  printf("Converting sparse array to full:");
  data=(scalar **) other;
  return *this;
}

template <class index_t, class scalar>
void full_matrix<index_t, scalar>::dimensions(index_t &mout, index_t &nout) const {
  mout=m;
  nout=n;
}

//should probably return a new matrix:
template <class index_t, class scalar>
void full_matrix<index_t, scalar>::transpose() {
  scalar **nm;
  scalar sw1;

  if (m==n) {
    matrix_transpose<scalar, index_t>(data, m);
  } else {
    nm=matrix_transpose<scalar, index_t>(data, m, n);
    delete_matrix(data);
    data=nm;
    sw1=m;
    m=n;
    n=sw1;
  }

}

//should probably return a new matrix:
template <class index_t, class scalar>
void full_matrix<index_t, scalar>::scal_mult(scalar cand) {
  for (long i=0; i<m*n; i++) data[0][i]=cand*data[0][i];
}
	
template <class index_t, class scalar>
int full_matrix<index_t, scalar>::scan(FILE *fs) {
  data=scan_matrix<scalar, index_t>(fs, m, n);
  //should do a bit better error handling here...
  if (data!=NULL) return 0;
  return FILE_READ_ERROR;
}

template <class index_t, class scalar>
size_t full_matrix<index_t, scalar>::read(FILE *fs) {
  if (data!=NULL) delete_matrix(data);
  data=read_matrix<scalar, index_t>(fs, m, n);
  if (data==NULL) return 0;
  return m*n+2;
}

template <class index_t, class scalar>
size_t full_matrix<index_t, scalar>::write(FILE *fs) {
  return write_matrix(fs, data, m, n);
}

template <class index_t, class scalar>
void full_matrix<index_t, scalar>::print(FILE *fs) {
  print_matrix(fs, data, m, n);
}

template <class index_t, class scalar>
matrix_base<index_t, scalar> * full_matrix<index_t, scalar>::clone() {
  full_matrix<index_t, scalar> *result;
  result=new full_matrix<index_t, scalar>(this);
  return result;

}

template <class index_t, class scalar>
scalar full_matrix<index_t, scalar>::norm() {
  return matrix_norm(data, m, n);
}

/*
template <class index_t, class data_t>
full_matrix<index_t, data_t>::operator sparse<index_t, data_t>& () {
  sparse<index_t, data_t> *result;
  result=new sparse<index_t, data_t>(*this);
  return *result;
}

template <class index_t, class data_t>
full_matrix<index_t, data_t>::operator sparse_array<index_t, data_t>& () {
  sparse_array<index_t, data_t> *result;
  result=new sparse_array<index_t, data_t>(*this);
  return *result;
}

template <class index_t, class data_t>
full_matrix<index_t, data_t>::operator full_matrix<index_t, data_t>& () {
  full_matrix<index_t, data_t> *result;
  result=new full_matrix(*this);
  return *result;
}
*/

//template class full_matrix<int16_t, float>;
template class full_matrix<ind_t, float>;
template class full_matrix<ind_t, double>;

} //end namespace libsparse
} //end namespace libpetey

