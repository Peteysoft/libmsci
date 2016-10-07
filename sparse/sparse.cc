//
// Copywrite 2004 Peter Mills.  All rights reserved.
//
//a simple class definition for dealing with sparse matrices
//
    
#include <math.h>
#include <complex.h>
#include <string.h>
    
#include "error_codes.h"
    
#include "../libpetey/heapsort_tmpl.cc"
#include "../libpetey/bin_search.cc"
    
#include "sparse_array.h"
#include "full_util.h"
#include "full_matrix.h"
    
#include "sparse_element.h"

using namespace std;
    
namespace libpetey {

  using namespace libsparse;
    
  template long bin_search< sparse_el<short, float> >
    		(sparse_el<short, float> *, long, sparse_el<short, float>, long);
  template long bin_search< sparse_el<long, float> >
   		(sparse_el<long, float> *, long, sparse_el<long, float>, long);
  template long bin_search< sparse_el<long, double> >
    		(sparse_el<long, double> *, long, sparse_el<long, double>, long);
  template void heapsort_inplace< sparse_el<short, float> >
    		(sparse_el<short, float> *, long);
  template void heapsort_inplace< sparse_el<long, float> >
    		(sparse_el<long, float> *, long);
  template void heapsort_inplace< sparse_el<long, double> >
    		(sparse_el<long, double> *, long);

  namespace libsparse {
    
    template <class index_t, class data_t>
    sparse<index_t, data_t>::sparse() {
      m=0;
      n=0;
      nel=0;
      array_size=1;
      last_search=-1;
      update_flag=1;
      eps=EPS;
      matrix=new sparse_el<index_t, data_t>[array_size];
    
      sparse_log=stderr;
    }
    
    //fucking stupid bullshit:
    template <class index_t, class scalar>
    sparse<index_t, scalar>::sparse(matrix_base<index_t, scalar> *other) {
      sparse_log=stderr;
      if (typeid(*other) == typeid(sparse<index_t, scalar>)) {
        sparse<index_t, scalar> *dum;
        dum=(sparse<index_t, scalar> *) other;
        copy(*dum);
      } else if (typeid(*other) == typeid(sparse_array<index_t, scalar>)) {
        convert(*(sparse_array<index_t, scalar> *) other);
      } else if (typeid(*other) == typeid(full_matrix<index_t, scalar>)) {
        full_matrix<index_t, scalar> *dum;
        dum=(full_matrix<index_t, scalar> *) other;
        matrix=new sparse_el<index_t, scalar>[1];
        array_size=1;
        from_full(dum->data, dum->m, dum->n);
      } else {
        fprintf(stderr, "Failed to initialize sparse matrix from base class\n");
      }
    
    };
    
    /*
    //fucking stupid bullshit:
    template <class index_t, class scalar>
    sparse<index_t, scalar> & sparse<index_t, scalar>::operator = (matrix_base<index_t, scalar> &other) {
      if (typeid(other) == typeid(sparse<index_t, scalar>)) {
        sparse<index_t, scalar> *dum;
        dum=(sparse<index_t, scalar> *) &other;
        return *this=*dum;
      } else if (typeid(other) == typeid(sparse_array<index_t, scalar>)) {
        return *this=*(sparse_array<index_t, scalar> *) &other;
      } else if (typeid(other) == typeid(full_matrix<index_t, scalar>)) {
        return *this=*(full_matrix<index_t, scalar> *) &other;
      }
    
    };
    */
    
    
    //this initializer specifies how small an element must be before it
    //is considered "zero"
    template <class index_t, class data_t>
    sparse<index_t, data_t>::sparse(data_t neps) {
      m=0;
      n=0;
      nel=0;
      array_size=1;
      last_search=-1;
      update_flag=1;
      eps=neps;
      matrix=new sparse_el<index_t, data_t>[array_size];
    
      sparse_log=stderr;
    }
    
    template <class index_t, class data_t>
    sparse<index_t, data_t>::~sparse() {
      delete[] matrix;
    }
    
    template<class index_t, class data_t>
    void sparse<index_t, data_t>::dimensions(index_t &mout, index_t &nout) const
    {
      mout=m;
      nout=n;
    }
    
    //turn sparse matrix into the identity matrix:
    template <class index_t, class data_t>
    void sparse<index_t, data_t>::identity() {
      index_t nd;
    
      if (m > n) nd=n; else nd=m;
    
      if (array_size < nd) {
        delete[] matrix;
        matrix=new sparse_el<index_t, data_t>[nd];
        array_size=nd;
      }
    
      for (index_t i=0; i<nd; i++) {
        matrix[i].i=i;
        matrix[i].j=i;
        matrix[i].value=1;
      }
      nel=nd;
    }
    
    template <class index_t, class data_t>
    void sparse<index_t, data_t>::identity(index_t mn, index_t nn) {
      m=mn;
      n=nn;
      identity();
    }
    
    //this initializer pre-specifies the dimensions of the matrix
    template <class index_t, class data_t>
    sparse<index_t, data_t>::sparse(index_t min, index_t nin, data_t neps) {
      m=min;
      n=nin;
      nel=0;
      array_size=1;
      last_search=-1;
      update_flag=1;
      eps=neps;
      matrix=new sparse_el<index_t, data_t>[array_size];
    
      sparse_log=stderr;
    }
    
    
    //this initializer loads up a full matrix
    template <class index_t, class data_t>
    void sparse<index_t, data_t>::from_full(data_t **non, index_t min, index_t nin, data_t neps) {
      data_t val;
    
      eps=neps;
      reset(min, nin);

      for (index_t i=0; i<m; i++) {
        for (index_t j=0; j<n; j++) {
          val=non[i][j];
          add_el(val, i, j);
        }
      }
      update();

    
    }
    
    //this initializer loads up a full matrix
    template <class index_t, class data_t>
    sparse<index_t, data_t>::sparse(data_t **non, index_t min, index_t nin, data_t neps) {
    
      array_size=1;
      eps=neps;
      matrix=new sparse_el<index_t, data_t>[array_size];
      from_full(non, min, nin, neps);
      sparse_log=stderr;
    }
    
    //this initializer loads up a full matrix
    template <class index_t, class data_t>
    sparse<index_t, data_t>::sparse(full_matrix<index_t, data_t> &non, data_t neps) {
      sparse_log=stderr;
      eps=neps;
      nel=0;
      array_size=1;
      matrix=new sparse_el<index_t, data_t>[array_size];
      from_full(non.data, non.m, non.n, neps);
    }
    
    //copy constructor:
    template <class index_t, class data_t>
    void sparse<index_t, data_t>::copy(const sparse<index_t, data_t> &other) {
      //delete[] matrix;
      m=other.m;
      n=other.n;
      nel=other.nel;
      eps=other.eps;
      if (nel>0) array_size=nel; else array_size=1;
      matrix=new sparse_el<index_t, data_t>[array_size];
      for (long i=0; i<nel; i++) matrix[i]=other.matrix[i];
    
      update_flag=other.update_flag;
      last_search=other.last_search;
    
    }
    
    //copy constructor:
    template <class index_t, class data_t>
    sparse<index_t, data_t>::sparse(const sparse<index_t, data_t> &other) {
      //delete[] matrix;
      sparse_log=stderr;
      copy(other);
    }
    
    //assignment operator:
    template <class index_t, class data_t>
    matrix_base<index_t, data_t> & sparse<index_t, data_t>::operator = (sparse<index_t, data_t> &other) {
      delete[] matrix;
    
      copy(other);
    
      return *this;
    }
    
    //clears matrix, but keeps all allocated memory intact:
    template <class index_t, class data_t>
    void sparse<index_t, data_t>::reset(index_t mnew, index_t nnew) {
      m=mnew;
      n=nnew;
      nel=0;
      last_search=-1;
      update_flag=1;
    }
    
    //clears matrix, frees memory:
    template <class index_t, class data_t>
    void sparse<index_t, data_t>::clear(index_t mnew, index_t nnew) {
      delete[] matrix;
    
      m=mnew;
      n=nnew;
      nel=0;
      array_size=1;
      last_search=-1;
      update_flag=1;
      matrix=new sparse_el<index_t, data_t>[1];
    }
    
    //extends reserved memory for non-zero elements by amount size:
    template <class index_t, class data_t>
    void sparse<index_t, data_t>::extend(long size) {
      long new_size;
      sparse_el<index_t, data_t> *new_matrix;
    
      new_size=size+nel;
      if (array_size < new_size) {
    
        array_size=new_size;
        new_matrix=new sparse_el<index_t, data_t>[array_size];
        for (long k=0; k<nel; k++) new_matrix[k]=matrix[k];
        delete[] matrix;
        matrix=new_matrix;
      }
    }
    
    //adds or changes a matrix element:
    template <class index_t, class data_t>
    long sparse<index_t, data_t>::add_el(data_t value, index_t i, index_t j) {
      sparse_el<index_t, data_t> *new_matrix;
      sparse_el<index_t, data_t> new_element(i, j, value);
    
      //if (fabs(value) < eps) return -1;
    
      //if the new element is larger than the stored dimensions,
      //extend them to match:
      if (i>=m) {
        m=i+1;
      }
      if (j>=n) {
        n=j+1;
      }
    
      if (fabs(value) < eps) return nel;
    
      //if there isn't enough space, double the number of available elements:
      if (nel >= array_size) {
    
        array_size*=2;
        new_matrix=new sparse_el<index_t, data_t>[array_size];
        for (long k=0; k<nel; k++) new_matrix[k]=matrix[k];
        delete[] matrix;
        matrix=new_matrix;
      }
    
      //simply tack the new element onto the end of the array:
      if (nel == 0) {
        matrix[nel]=new_element;
        nel++;
      } else {
        if (new_element == matrix[nel-1]) {
          matrix[nel-1].value=new_element.value;
        } else {
          matrix[nel]=new_element;
          if (update_flag != 0) if (new_element <= matrix[nel-1]) {
            update_flag=0;
    //      printf("Update flag changed after %d elements\n", nel);
          }
          nel++;
        }
      }
    
      //if (update_flag==1) if (fabs(new_element.value) < eps) update_flag=2;
    
      return nel;
    }
    
    //change an existing element; if it does not exist, add it on:
    template <class index_t, class data_t>
    long sparse<index_t, data_t>::cel(data_t value, index_t i, index_t j) {
      sparse_el<index_t, data_t> new_el(i, j, value);
      sparse_el<index_t, data_t> *new_matrix;
    
      //sort the matrix elements (row first, column next):
      update();
    
      //search for the element:
      if (nel > 0) {
        last_search=bin_search(matrix, nel, new_el, last_search);
      } else {
        last_search=-1;
      }
    
      if (last_search==-1 || matrix[last_search]!=new_el) {
        add_el(value, i, j);
      } else {
        matrix[last_search].value=value;
      }
    
      return last_search;
    }
    
    //index the matrix:
    //if the desired element does not exist or the indices are out of
    //bounds, simply return 0
    template <class index_t, class data_t>
    data_t sparse<index_t, data_t>::operator ( ) (index_t i, index_t j) {
      sparse_el<index_t, data_t> search_blank(i, j, 0);
    
      if (nel == 0) return 0;
      update();
      last_search=bin_search(matrix, nel, search_blank, last_search);
      if (last_search == -1 || last_search >= nel) return 0;
      if (search_blank != matrix[last_search]) return 0;
    
      return matrix[last_search].value;
    }
    
    //sure, we can do it just with matrix multiplication, but this is
    //much more efficient:
    template <class index_t, class data_t>
    void sparse<index_t, data_t>::get_row (index_t i, data_t *row) {
      sparse_el<index_t, data_t> search_blank(i, 0, 0);
    
      for (index_t j=0; j<n; j++) row[j]=0;
      if (nel == 0) return;
      update();
      last_search=bin_search(matrix, nel, search_blank, last_search);
      if (last_search >= nel) return;
      if (last_search==-1 || matrix[last_search].i<i) {
        last_search++;
      }
    
      for (index_t j=last_search; j<nel; j++) {
        if (matrix[j].i > i) break;
        row[matrix[j].j]=matrix[j].value;
      }
    }
    
    template <class index_t, class data_t>
    data_t *sparse<index_t, data_t>::operator ( ) (index_t i) {
      data_t *row;
      row=new data_t[n];
      get_row(i, row);
      return row;
    }
    
    //for when the matrix is already sorted: remove zero elements...
    template <class index_t, class data_t>
    void sparse<index_t, data_t>::remove_zeros (data_t minmag) {
      long index[nel];	//index of non-zero elements
      long nnon;		//number of non-zero elements
      sparse_el<index_t, data_t> *new_matrix;
      
      if (update_flag==0) {
        update();
      } else {
        update_flag=1;
        nnon=0;
        for (long i=0; i<nel; i++) {
          if (fabs(matrix[i].value) > minmag) {
            index[nnon]=i;
            nnon++;
          }
        }
        if (nnon>0) {
          new_matrix=new sparse_el<index_t, data_t> [nnon];
          array_size=nnon;
        } else {
          new_matrix=new sparse_el<index_t, data_t>[1];
          array_size=1;
        }
        for (long i=0; i<nnon; i++) new_matrix[i]=matrix[index[i]];
        delete [] matrix;
        matrix=new_matrix;
        fprintf(sparse_log, "Found %ld zero or insignificant elements\n", nel-nnon);
        nel=nnon;
      }
    }
    
    //update the matrix:
    //-sort the elements by row first, then column
    //-remove elements less than eps
    //-remove duplicate elements (caveat: the order of insertion
    // does not effect which one is removed)
    //-last point has been fixed
    template <class index_t, class data_t>
    void sparse<index_t, data_t>::update() {
      long i;
      //long keep[nel];		//indices of those elements to keep
      long index[nel];
      //long nkeep;
      long offset;
      long nins, ndup;		//number of insignificant and duplicate elements respectively
      long last_el;
      long maxind;
      sparse_el<index_t, data_t> *new_matrix;
    
      if (update_flag==2) remove_zeros();
      if (update_flag==1) return;
      if (nel == 0) {
        update_flag==1;
        return;
      }
    
      //sort the elements:
      heapsort(matrix, (long *) index, nel);
    
      //remove zero and duplicate elements:
      maxind=index[0];
      /*
      if (fabs(matrix[index[0]].value) > eps) {
        keep[0]=index[0];
        nkeep++;
      }
      */
      new_matrix=new sparse_el<index_t, data_t>[array_size];
    
      nins=0;
      ndup=0;
      for (long i=1; i<nel; i++) {
        if (matrix[index[i]] != matrix[index[i-1]]) {
          if (fabs(matrix[maxind].value) > eps) {
            new_matrix[i-nins-ndup-1]=matrix[maxind];
          } else {
            nins++;
          }
          //keep[nkeep]=index[i];
          //nkeep++;
          maxind=index[i];
        } else {
          //printf("%d %d %f\n", matrix[index[i-1]].i, matrix[index[i-1]].j, matrix[index[i-1]].value);
          //printf("%d %d %f\n", matrix[index[i]].i, matrix[index[i]].j, matrix[index[i]].value);
          //matrix element is a duplicate:
          //see if it was added later or earlier than current latest duplicate
          if (index[i] > maxind) {
            maxind=index[i];
          }
          ndup++;
        }
      }
      
      offset=nins+ndup;
      if (fabs(matrix[maxind].value) > eps) {
    //    printf("%d %d %f\n", matrix[maxind].i, matrix[maxind].j, matrix[maxind].value);
        new_matrix[nel-offset-1]=matrix[maxind];
      } else {
        nins++;
      }
    
      fprintf(sparse_log, "Matrix updated: %ld insignificant and %ld duplicate elements found\n", nins, ndup);
    
      nel=nel-offset;
    
      delete[] matrix;
      matrix=new_matrix;
    
      if (nel !=0) {
        if (2*nel <= array_size) {
          array_size=nel;
          new_matrix=new sparse_el<index_t, data_t>[array_size];
          for (long i=0; i<nel; i++) {
            new_matrix[i]=matrix[i];
          }
          delete[] matrix;
          matrix=new_matrix;
        }
      } else {
        delete[] matrix;
        array_size=1;
        matrix=new sparse_el<index_t, data_t>[array_size];
      }
    
      update_flag=1;
    }
    
    //perform a vector multiplication:
    template <class index_t, class data_t>
    void sparse<index_t, data_t>::vect_mult(data_t *cand, data_t *result) {
      for (index_t i=0; i<m; i++) result[i]=0;
      for (long i=0; i<nel; i++) {
        result[matrix[i].i]+=matrix[i].value*cand[matrix[i].j];
      }
    }
    
    //perform left vector multiplication:
    template <class index_t, class data_t>
    void sparse<index_t, data_t>::left_mult(data_t *cor, data_t *result) {
      for (index_t i=0; i<n; i++) result[i]=0;
      for (long i=0; i<nel; i++) {
        result[matrix[i].j]+=matrix[i].value*cor[matrix[i].i];
      }
    }
    
    //perform a matrix multiplication on a non-sparse matrix:
    template <class index_t, class data_t>
    void sparse<index_t, data_t>::mat_mult(data_t **cand, data_t **result,
    			index_t np) {		//outer dimensions
    
      for (index_t i=0; i<m; i++) {
        for (index_t j=0; j<np; j++) {
          result[i][j]=0;
        }
      }
    
      for (index_t j=0; j<np; j++) {
        for (long k=0; k<nel; k++) {
          result[matrix[k].i][j]+=matrix[k].value*cand[matrix[k].j][j];
        }
      }
    }
    
    //perform left matrix multiplication on a non-sparse matrix:
    template <class index_t, class data_t>
    void sparse<index_t, data_t>::left_m_mult(data_t **cor, data_t **result,
    			index_t np) {		//outer dimensions
    
      for (index_t i=0; i<np; i++) {
        for (index_t j=0; j<n; j++) {
          result[i][j]=0;
        }
      }
    
      for (index_t i=0; i<np; i++) {
        for (long k=0; k<nel; k++) {
          result[i][matrix[k].j]+=matrix[k].value*cor[i][matrix[k].i];
        }
      }
    }
    
    //perform a matrix multiplication on a sparse matrix:
    template <class index_t, class data_t>
    void sparse<index_t, data_t>::mat_mult(sparse<index_t, data_t> &cand, sparse<index_t, data_t> &result) {
      sparse<index_t, data_t> tran_cand;
    
      //transpose the multiplicand:
      cand.transpose(tran_cand);
    
      //perform the matrix multiplication:
      mat_mult_t(tran_cand, result);
    }
    
    //multiplies with the transpose of a sparse matrix:
    template <class index_t, class data_t>
    void sparse<index_t, data_t>::mat_mult_t(sparse<index_t, data_t> &cand, sparse<index_t, data_t> &result) {
      data_t new_val;
      long q, r;
      long qold;
      index_t this_i_old, cand_i_old;
      index_t this_j, cand_j;
    
    //  printf("Updating multiplier:\n");
      update();
    
      result.nel=0;
      result.last_search=-1;
      result.update_flag=2;
    
    //  printf("Updating multiplicand:\n");
      cand.update();
      result.m=m;
      result.n=cand.m;
    
    //  cand.print(stdout);
    
      q=0;			//index for multiplier
      r=0;			//index for multiplicand
      this_i_old=matrix[0].i;	//current row subscript for multiplier
      cand_i_old=cand.matrix[0].i;	//current row subscript for multiplicand
      new_val=0;
      qold=0;
      do {
        this_j=matrix[q].j;		//current column subscript for multiplier
        cand_j=cand.matrix[r].j;	//current column subscript for multiplicand
    
    //    printf("q= %d r= %d, i1= %d j1= %d, i2= %d j2= %d\n", q, r,
    //    		this_i_old, this_j, cand_i_old, cand_j);
    
        //search for duplicate columns:
        if (this_j == cand_j) {
    //      printf("Hit!\n");
          //duplicate column found, increment both indices:
          new_val+=matrix[q].value*cand.matrix[r].value;
          q++;
          r++;
        } else if (this_j < cand_j) {
          //column of multiplicand is higher, increment multiplier index:
          q++;
        } else {
          //column of multiplier is higher, increment multiplicand index:
          r++;
        }
    
        if (q >= nel) {
          //index for multiplier is at the end, add new value:
          if (fabs(new_val) > eps) result.add_el(new_val, this_i_old, cand_i_old);
          new_val=0;
          //increment muliplicand index until it reaches the end, or row index
          //changes:
          while (cand_i_old == cand.matrix[r].i) {
            r++;
            if (r >= cand.nel) break;
          }
          //if multiplicand index is at the end, we are finished,
          //otherwise, proceed to all the subsequent rows of the multiplicand
          if (r >= cand.nel) {
            break;
          } else {
            q=qold;
    	cand_i_old=cand.matrix[r].i;
          }
          //other cases are dealt with in a similar fashion:
        } else if (r >= cand.nel) {
          if (fabs(new_val) > eps) result.add_el(new_val, this_i_old, cand_i_old);
          new_val=0;
          while (this_i_old == matrix[q].i) {
    	q++;
            if (q >= nel) break;
          }
          if (q >= nel) {
            break;
          } else {
            qold=q;
            this_i_old=matrix[q].i;
            r=0;
            cand_i_old=cand.matrix[0].i;
          }
        } else if (matrix[q].i != this_i_old) {
          if (fabs(new_val) > eps) result.add_el(new_val, this_i_old, cand_i_old);
          new_val=0;
          while (cand_i_old == cand.matrix[r].i) {
    	r++;
            if (r >= cand.nel) break;
          }
          if (r >= cand.nel) {
            qold=q;
    	this_i_old=matrix[q].i;
            r=0;
            cand_i_old=cand.matrix[0].i;
          } else {
            q=qold;
    	cand_i_old=cand.matrix[r].i;
          }
        } else if (cand_i_old != cand.matrix[r].i) {
          if (fabs(new_val) > eps) result.add_el(new_val, this_i_old, cand_i_old);
          new_val=0;
          q=qold;
          cand_i_old=cand.matrix[r].i;
        }
      } while (1);
    
    }
    
    //note: here this receives the addition, whereas in the full add,
    //the full matrix (i.e, the parameter) receives the addition
    template <class index_t, class data_t>
    void sparse<index_t, data_t>::sparse_add(sparse<index_t, data_t> &b) {
      update();
      b.update();
      long i, j;
    
      i=0;
      j=0;
    
      if (nel == 0) {
        delete [] matrix;
        copy(b);
        return;
      }
    
      while (j<b.nel) {
        if (matrix[i]>b.matrix[j]) {
          add_el(b.matrix[j].value, b.matrix[j].i, b.matrix[j].j);
          j++;
        } else if (matrix[i]<b.matrix[j]) {
          i++;
        } else {
          matrix[i].value+=b.matrix[j].value;
          i++;
          j++;
        }
      }
        
    }
    
    //note here the full matrix (i.e, the parameter) receives the addition
    //whereas in the sparse add, this receives it...
    template <class index_t, class data_t>
    void sparse<index_t, data_t>::full_add(data_t **b) {
      update();
      for (long i=0; i<nel; i++) {
        b[matrix[i].i][matrix[i].j]+=matrix[i].value;
      }
    }
    
    //convert matrix to its transpose:
    template <class index_t, class data_t>
    void sparse<index_t, data_t>::transpose() {
      index_t sw;
      sw=m;
      m=n;
      n=sw;
    
      for (long i=0; i<nel; i++) {
        sw=matrix[i].i;
        matrix[i].i=matrix[i].j;
        matrix[i].j=sw;
      }
      update_flag=0;
    }
    
    
    //find the transpose and return to a sparse matrix:
    template <class index_t, class data_t>
    void sparse<index_t, data_t>::transpose(sparse<index_t, data_t> &tran) {
    
      delete[] tran.matrix;
    
      tran.m=n;
      tran.n=m;
      tran.nel=nel;
      tran.array_size=nel;
      last_search=-1;
      tran.matrix=new sparse_el<index_t, data_t>[array_size];
    
      for (long i=0; i<nel; i++) {
        tran.matrix[i].i=matrix[i].j;
        tran.matrix[i].j=matrix[i].i;
        tran.matrix[i].value=matrix[i].value;
      }
      
      tran.update_flag=0;
      
    }
    
    //find the transpose and return to a full matrix:
    template <class index_t, class data_t>
    void sparse<index_t, data_t>::transpose(data_t **tran) {
    
      for (index_t i=0; i<n; i++) {
        for (index_t j=0; j<m; j++) {
          tran[i][j]=0;
        }
      }
      for (long i=0; i<nel; i++) {
        tran[matrix[i].j][matrix[i].i]=matrix[i].value;
      }
    }
    
    template <class index_t, class data_t>
    void sparse<index_t, data_t>::sl_transform(data_t m, data_t b) {
      index_t n1;
      long ind=-1;
      sparse_el<index_t, data_t> search_blank;
    
      update();
    
      if (n>m) n1=m; else n1=n;
      for (long i=0; i<nel; i++) matrix[i].value*=m;
      for (index_t i=0; i<n1; i++) {
        search_blank.i=i;
        search_blank.j=i;
        ind=bin_search(matrix, nel, search_blank, ind);
        matrix[ind].value+=b;
      }
    }
    
    /*
    //find the transpose and return as a full matrix:
    float ** sparse<index_t, data_t>::transpose() {
      float *temp, **tran;
    
      temp=new float[m*n];
      tran=new float *[m];
    
      for (index_t i=0; i<m; i++) {
        tran[i]=&temp[i*n];
        for (index_t j=0; j<n; j++) {
          tran[i][j]=0;
        }
      }
    
      for (long i=0; i<nel; i++) {
        tran[matrix[i].j][matrix[i].i]=matrix[i].value;
      }
    
      return tran;
    }
    */
    
    //convert to a full matrix:
    template <class index_t, class data_t>
    void sparse<index_t, data_t>::full(data_t ** non) {
      for (index_t i=0; i<m; i++) {
        for (index_t j=0; j<n; j++) {
          non[i][j]=0;
        }
      }
    
      for (long i=0; i<nel; i++) {
        non[matrix[i].i][matrix[i].j]=matrix[i].value;
      }
    
    }
    
    //convert to a full matrix:
    template <class index_t, class data_t>
    sparse<index_t, data_t>::operator data_t ** () {
      data_t *temp, **non;
    
      update();
    
      temp=new data_t[m*n];
      non=new data_t *[m];
    
      for (index_t i=0; i<m; i++) {
        non[i]=&temp[i*n];
        for (index_t j=0; j<n; j++) {
          non[i][j]=0;
        }
      }
    
      for (long i=0; i<nel; i++) {
        non[matrix[i].i][matrix[i].j]=matrix[i].value;
      }
    
      return non;
    }
    
    template <class index_t, class data_t>
    size_t sparse<index_t, data_t>::read(FILE *fptr) {
      size_t nread;
      if (feof(fptr)!=0) return 0;
      nread=fread(&m, sizeof(m), 1, fptr);
      nread+=fread(&n, sizeof(n), 1, fptr);
      nread+=fread(&nel, sizeof(nel), 1, fptr);
      delete[] matrix;
      if (nel>0) array_size=nel; else array_size=1;
      matrix=new sparse_el<index_t, data_t>[array_size];
      nread+=fread(matrix, sizeof(sparse_el<index_t, data_t>), nel, fptr);
      update_flag=1;
    
      return nread;
    }
    
    template <class index_t, class data_t>
    size_t sparse<index_t, data_t>::write(FILE *fptr) {
      size_t nwritten;
    
      update();
      nwritten=0;
      nwritten+=fwrite(&m, sizeof(m), 1, fptr);
      nwritten+=fwrite(&n, sizeof(n), 1, fptr);
      nwritten+=fwrite(&nel, sizeof(nel), 1, fptr);
      if (nel > 0) {
        nwritten+=fwrite(matrix, sizeof(sparse_el<index_t, data_t>), nel, fptr);
      }
    
      return nwritten;
    }
    
    #define MAXLL 1000
    
    template <class index_t, class data_t>
    int sparse<index_t, data_t>::scan(FILE *fptr) {
      int32_t in, jn;
      data_t val;
      int err=0;
      char line[MAXLL];
      long fpos;
      char format[10];

      if (sizeof(data_t)==8) {
        strcpy(format, "%d %d %lg");
      } else {
        strcpy(format, "%d %d %g");
      }
    
      if (fgets(line, MAXLL, fptr)==NULL) return FILE_READ_ERROR;
      sscanf(line, "%d %d", &in, &jn);
      clear(in, jn);
    
      for (long i=0; i<m*n; i++) {
        fpos=ftell(fptr);
        if (fgets(line, MAXLL, fptr)==NULL) break;
        if (sscanf(line, format, &in, &jn, &val)<3) {
          fprintf(stderr, "sparse::scan -- error while scanning line %ld\n", (int64_t) i);
          fseek(fptr, fpos, SEEK_SET);
          err=FILE_READ_WARNING;
          break;
        }
        add_el(val, in, jn);
      }
      return 0;
    }
    
    template <class index_t, class data_t>
    void sparse<index_t, data_t>::print(FILE *fptr) {
      update();
    
      fprintf(fptr, "%d %d\n", m, n);
      for (long i=0; i<nel; i++) {
        fprintf(fptr, "%d %d %g\n", matrix[i].i, matrix[i].j, matrix[i].value);
      }
    }
    
    template <class index_t, class data_t>
    void sparse<index_t, data_t>::print_raw(FILE *fptr) {
    
      fprintf(fptr, "m=%d n=%d\n", m, n);
      for (long i=0; i<nel; i++) {
        fprintf(fptr, "%d %d %g\n", matrix[i].i, matrix[i].j, matrix[i].value);
      }
    }
    
    template <class index_t, class data_t>
    float sparse<index_t, data_t>::storage_ratio() {
      float ratio;
    
      ratio=1.0*array_size*sizeof(sparse_el<index_t, data_t>)/(m*n*sizeof(data_t));
    
      return ratio;
    }
    
    template <class index_t, class data_t>
    float sparse<index_t, data_t>::performance_ratio() {
      float cfull=3.25e-8;		//coefficient for full matrix
      float c2=5.47e-8;		//factor for calculating sparsity dependent coefficient
      float k=-0.64;		//exponent for sparsity dependent coeff.
      float sparsity;
      float csparse;
      float sfull, ssparse;		//speed of full and sparse matrices resp.
    
      //calculate perf. of equiv. full matrix:
      sfull=cfull*pow((m*n), (3./2));
    
    //  printf("Nominal matrix dimensions: %d %d\n", m, n);
    
      //calculate sparsity:
      sparsity=1.0*nel/(m*n);
    
      //calculate sparsity dependent coefficient:
      csparse=c2*pow(sparsity, k);
    
      //calculate the performance of the sparse matrix:
      ssparse=csparse*pow(nel, (3./2));
    
    //  printf("~s full: %f, sparsity: %f, c sparse: %g, ~s sparse: %f\n", sfull, sparsity, csparse, ssparse);
    
      return ssparse/sfull;
    }
    
    //some instantiations:
    //template class sparse<int16_t, float>;
    template class sparse<int32_t, float>;
    template class sparse<int32_t, double>;
    //template class sparse<long, complex<float>>;
    //template class sparse<long, complex<double>>;
   
  } //end namespace libsparse
} //end namespace libpetey
 
