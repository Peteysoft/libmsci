#ifndef SPARSE_H
#define SPARSE_H 1

#include <stdio.h>

#include <stdint.h>

#include "bit_array.h"
#include "matrix_base.h"

namespace libpetey {
  namespace libsparse {

    //forward declaration of sparse_el:
    template <class index_t, class data_t> class sparse_el;

    typedef int32_t ind_t;

    #define EPS 5.96e-8		//this is actually pretty big... (?)

    template <class index_t, class data_t>
    class sparse:public matrix_base<index_t, data_t> {
      //friend void conj_grad_norm(sparse<ind_type, float> *, data_t *, data_t *, data_t, long);
      //friend double * irreg2reg(double **xmat, long m, double *y, long n, double **grids, long *ngrid);
      private:
        //since all the "constructors" do the same thing, we have one private
        //routine to do them all:
        void copy(const sparse<index_t, data_t> &other);
        void convert(sparse_array<index_t, data_t> &other, data_t neps=EPS);
    
      protected:
        sparse_el<index_t, data_t> *matrix;

        //dimensions:
        index_t m;
        index_t n;

        //number of nonzero elements:		(these types should depend on index_t
        long nel;				//i.e. index_t = int32_t, then use int64_t)

        //index of last element searched:
        long last_search;

        //size of array
        long array_size;

        //update flag:
        //0=matrix needs updating
        //1=matrix is fully updated
        //2=matrix is sorted but has zero elements
        char update_flag;

        //"zero" tolerance:
        float eps;

        FILE *sparse_log;

      public:
        //initializers:
        sparse();
        sparse(data_t neps);
        sparse(index_t min, index_t nin, data_t neps=EPS);
        sparse(data_t **non, index_t min, index_t nin, data_t neps=EPS);
        virtual ~sparse();

        //create the identity matrix with existing dimensions:
        void identity();
        //create the identity matrix with specified dimensions:
        void identity(index_t mn, index_t nn);

        sparse(matrix_base<index_t, data_t> *other);
        sparse(const sparse<index_t, data_t> &other);
        sparse(sparse_array<index_t, data_t> &other);
        sparse(full_matrix<index_t, data_t> &other, data_t neps=EPS);

        void from_full(data_t **non, index_t min, index_t nin, data_t neps=EPS);

        //storage manipulation:
        void update();		//update matrix
        void remove_zeros();	//remove insignificant elements
        //remove elemnts of less than specified magnitude:
        void remove_zeros(data_t minmag);
        //clear matrix, but do not change storage:
        void reset(index_t mnew=0, index_t nnew=0);
        //clear matrix, deallocate storage:
        void clear(index_t mnew=0, index_t nnew=0);
        //extend storage by specified amount:
        void extend(long size);

        //add and access matrix elements:
        //add or change an element:
        long add_el(data_t val, index_t i, index_t j);
        //change existing element or insert new:
        virtual long cel(data_t val, index_t i, index_t j);

        //return value of element:
        virtual data_t operator ( ) (index_t i, index_t j);

        //access rows:
        virtual data_t *operator ( ) (index_t i);
        virtual void get_row(index_t i, data_t *row);

        //arithmetic operations:
        //multiply with a vector:
        virtual void vect_mult(data_t *cand, data_t *result);
        //multiply with a non-sparse matrix:
        void mat_mult(data_t **cand, data_t **result, index_t np);
        //multiply with sparse:
        void mat_mult(sparse<index_t, data_t> &cand, sparse<index_t, data_t> &result);
        //multiply with transpose of sparse:
        void mat_mult_t(sparse<index_t, data_t> &cand, sparse<index_t, data_t> &result);

        //add a sparse matrix:
        void sparse_add(sparse<index_t, data_t> &b);
        //adds the sparse matrix to the full matrix:
        void full_add(data_t **b);

        //left multiply:
        void left_m_mult(data_t **cor, data_t **result, index_t np);
        virtual void left_mult(data_t *cor, data_t *result);

        //take the transpose:
        void transpose();
        void transpose(sparse<index_t, data_t> &trans);	//returns sparse
        void transpose(data_t **non);			//returns non-sparse
        //float ** transpose();				//returns non-sparse

        void full(data_t ** full);
        operator data_t ** ();		//convert to non-sparse

        //informational:
        //return matrix dimensions
        virtual void dimensions(index_t &mout, index_t &nout) const;
        //return number of non-zero elements
        long size();

        //linear transform of the form: m*A+b*I
        void sl_transform(data_t m, data_t b);

        //ratio of amount of storage to equivalent full:
        float storage_ratio();
        //approx. ratio of performance (for matrix multiply) to equivalent full:
        float performance_ratio();

        //I/O:
        virtual size_t read(FILE *fptr);		//binary read
        virtual size_t write(FILE *ptr);		//binary write
        virtual void print(FILE *ptr);		//print (ascii) to file
        virtual int scan(FILE *ptr);		//read from (ascii) file
        void print_raw(FILE *ptr);

        //the "canonical" versions:
        virtual matrix_base<index_t, data_t> * mat_mult(
		    matrix_base<index_t, data_t> *cand);
        virtual matrix_base<index_t, data_t> * add(matrix_base<index_t, data_t> *b);

        virtual data_t * vect_mult(data_t *cand);
        virtual void scal_mult(data_t m);
        virtual data_t * left_mult(data_t *cor);

        virtual matrix_base<index_t, data_t> & operator = (full_matrix<index_t, data_t> &other);
        virtual matrix_base<index_t, data_t> & operator = (sparse<index_t, data_t> &other);
        virtual matrix_base<index_t, data_t> & operator = (sparse_array<index_t, data_t> &other);

        virtual data_t norm();

        virtual matrix_base<index_t, data_t> * clone();

       //sparse<index_t, data_t> & operator = (matrix_base<index_t, data_t> &other);
 
/*  
        virtual operator sparse<index_t, data_t>& ();
        virtual operator sparse_array<index_t, data_t>& ();
        virtual operator full_matrix<index_t, data_t>& ();
*/

    };

    template<class index_t, class data_t>
    inline long sparse<index_t, data_t>::size() {
      return nel;
    }

    template<class index_t, class data_t>
    inline void sparse<index_t, data_t>::remove_zeros() {
      remove_zeros(eps);
    }

    typedef sparse<ind_t, float> sparse_matrix;

  }  //end namespace libsparse
}    //end namespace libpetey

#endif

