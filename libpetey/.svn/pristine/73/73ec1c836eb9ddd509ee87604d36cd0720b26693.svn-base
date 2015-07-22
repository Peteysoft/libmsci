#ifndef MAT_PRE_H
#define MAT_PRE_H

#include <stdio.h>
#include <typeinfo>

#include "matrix_base.h"
#include "sparse_calc_defs.h"

template <class index_t, class scalar>
class mat_pre: public matrix_base<index_t, scalar> {
  protected:
    index_t m;
    index_t n;

  public:
    mat_pre();
    virtual ~mat_pre();
    virtual ~mat_pre(index_t m1, index_t n1);

    //access elements:
    virtual scalar operator ( ) (index_t i, index_t j);
    virtual long cel (scalar val, index_t i, index_t j);

    //access rows:
    virtual scalar *operator ( ) (index_t i);
    virtual void get_row(index_t i, scalar *row);

    //matrix multiplication:
    virtual matrix_base<index_t, scalar> * mat_mult(full_matrix<index_t, scalar> *cand);
    virtual matrix_base<index_t, scalar> * mat_mult(sparse<index_t, scalar> *cand);
    virtual matrix_base<index_t, scalar> * mat_mult(sparse_array<index_t, scalar> *cand);
    virtual matrix_base<index_t, scalar> * mat_mult(mat_pre<index_t, scalar> *cand);

    //matrix addition:
    virtual matrix_base<index_t, scalar> * add(full_matrix<index_t, scalar> * b);
    virtual matrix_base<index_t, scalar> * add(sparse<index_t, scalar> * b);
    virtual matrix_base<index_t, scalar> * add(sparse_array<index_t, scalar> * b);
    virtual matrix_base<index_t, scalar> * add(mat_pre<index_t, scalar> * b);

    //vector multiplication:
    virtual scalar * vect_mult(scalar *cand);
    virtual scalar * left_mult(scalar *cor);

    //the above are good for the sparse calculator, but not for solvers
    virtual void vect_mult(scalar *cand, scalar *result);
    virtual void left_mult(scalar *cor, scalar *result);

    //scalar multiplication:
    virtual void scal_mult(scalar cand);

    //take transpose: 
    virtual void transpose();

    //informational:
    virtual void dimensions(index_t &m, index_t &n) const;

    //norm:
    virtual scalar norm();

    //IO:
    virtual size_t read(FILE *fs, int mattype=FULL_T);
    virtual size_t write(FILE *fs);
    virtual void print(FILE *fs);
    virtual int scan(FILE *fs);

    //type conversion:
    virtual matrix_base<index_t, scalar> & operator = (full_matrix<index_t, scalar> &other);
    virtual matrix_base<index_t, scalar> & operator = (sparse<index_t, scalar> &other);
    virtual matrix_base<index_t, scalar> & operator = (sparse_array<index_t, scalar> &other);
    virtual matrix_base<index_t, scalar> & operator = (mat_pre<index_t, scalar> &other);

    //need this so we can use the "matrix_base" class in the solver routines:
    virtual matrix_base<index_t, scalar> * clone();

    //virtual operator sparse<index_t, scalar>& ()=0;
    //virtual operator sparse_array<index_t, scalar>& ()=0;
    //virtual operator full_matrix<index_t, scalar>& ()=0;

};

#endif

