#ifndef FULL_MATRIX_H
#define FULL_MATRIX_H

#include "matrix_base.h"

namespace libpetey {
namespace libsparse {

template <class index_t, class scalar>
class full_matrix:public matrix_base<index_t, scalar> {
  protected:
    index_t m;
    index_t n;
    scalar **data;
  public:
    friend class sparse<index_t, scalar>;
    //destroys a large piece of the functionality of the sparse_array class:
    friend class sparse_array<index_t, scalar>;

    full_matrix();
    full_matrix(index_t min, index_t nin);
    full_matrix(scalar **dat, index_t min, index_t nin);
    full_matrix(scalar *dat, index_t min, index_t nin);

    full_matrix(matrix_base<index_t, scalar> *other);

    full_matrix(full_matrix<index_t, scalar> &other);
    full_matrix(sparse<index_t, scalar> &other);
    full_matrix(sparse_array<index_t, scalar> &other);

    void identity();
    void ones();

    virtual ~full_matrix();

    //access elements:
    virtual scalar operator ( ) (index_t i, index_t j);
    virtual long cel(scalar val, index_t i, index_t j);

    //access rows:
    virtual scalar *operator ( ) (index_t i);
    virtual void get_row(index_t i, scalar *row);

    //matrix multiplication:
    virtual matrix_base<index_t, scalar> * mat_mult(matrix_base<index_t, scalar> *cand);

    //matrix addition:
    virtual matrix_base<index_t, scalar> * add(matrix_base<index_t, scalar> * b);

    //vector multiplication:
    virtual scalar * vect_mult(scalar *cand);
    virtual scalar * left_mult(scalar *cor);

    virtual void vect_mult(scalar *cand, scalar *result);
    virtual void left_mult(scalar *cor, scalar *result);

    //scalar multiplication:
    virtual void scal_mult(scalar cand);

    //take transpose:
    virtual void transpose();

    //informational:
    virtual void dimensions(index_t &mout, index_t &nout) const;

    //norm:
    virtual scalar norm();

    //IO:
    virtual size_t read(FILE *fs);
    virtual size_t write(FILE *fs);
    virtual void print(FILE *fs);
    virtual int scan(FILE *fs);

    //type conversion:
    virtual matrix_base<index_t, scalar> & operator = (full_matrix<index_t, scalar> &other);
    virtual matrix_base<index_t, scalar> & operator = (sparse<index_t, scalar> &other);
    virtual matrix_base<index_t, scalar> & operator = (sparse_array<index_t, scalar> &other);

    virtual matrix_base<index_t, scalar> * clone();

    //virtual full_matrix<index_t, scalar> & operator = (matrix_base<index_t, scalar> &other);

/*
    virtual operator sparse<index_t, scalar>& ();
    virtual operator sparse_array<index_t, scalar>& ();
    virtual operator full_matrix<index_t, scalar>& ();
*/

};

} //end namespace libsparse
} //end namespace libpetey

#endif

