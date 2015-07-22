#ifndef MATRIX_BASE_H
#define MATRIX_BASE_H

#include <stdio.h>
#include <typeinfo>

//"null" matrix: i.e. pre-compute the dimensions for error-checking purposes

namespace libpetey {
  namespace libsparse {

    //forward declarations:
    template <class index_t, class scalar> class matrix_base;
    template <class index_t, class scalar> class sparse;
    template <class index_t, class scalar> class full_matrix;
    template <class index_t, class scalar, class sparse_t=sparse<index_t, scalar> > class sparse_array;

    template <class index_t, class scalar>
    class matrix_base {
      private:
        index_t m;
        index_t n;

      public:
        matrix_base();
        matrix_base(index_t m1, index_t n1);
        matrix_base(matrix_base<index_t, scalar> &other);
        matrix_base(matrix_base<index_t, scalar> *other);
        virtual ~matrix_base();

        //access elements:
        virtual scalar operator ( ) (index_t i, index_t j);
        virtual long cel (scalar val, index_t i, index_t j);

        //access rows:
        virtual scalar *operator ( ) (index_t i);
        virtual void get_row(index_t i, scalar *row);

        //matrix multiplication:
        virtual matrix_base<index_t, scalar> * mat_mult(matrix_base<index_t, scalar> *cand);

        //matrix addition:
        virtual matrix_base<index_t, scalar> * add(matrix_base<index_t, scalar> *other);

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
        virtual size_t read(FILE *fs);
        size_t read(FILE *fs, int mtype);
        virtual size_t write(FILE *fs);
        virtual void print(FILE *fs);
        virtual int scan(FILE *fs);

        //type conversion:
        virtual matrix_base<index_t, scalar> & operator = (full_matrix<index_t, scalar> &other);
        virtual matrix_base<index_t, scalar> & operator = (sparse<index_t, scalar> &other);
        virtual matrix_base<index_t, scalar> & operator = (sparse_array<index_t, scalar> &other);

        matrix_base<index_t, scalar> & operator = (matrix_base<index_t, scalar> &other);

        //need this so we can use the "matrix_base" class in the solver routines:
        virtual matrix_base<index_t, scalar> * clone();

        //virtual operator sparse<index_t, scalar>& ()=0;
        //virtual operator sparse_array<index_t, scalar>& ()=0;
        //virtual operator full_matrix<index_t, scalar>& ()=0;

    };

    //to avoid creating a "clone" method--the designers of C++ really should
    //add polymorphic copying to the standard...
    template <class index_t, class scalar>
    matrix_base<index_t, scalar> * copy_matrix(matrix_base<index_t, scalar> *other);

  } //end namespace libsparse
} //end namespace libpetey

#endif

