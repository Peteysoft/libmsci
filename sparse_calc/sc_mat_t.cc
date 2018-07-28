#include "vector_s.h"
#include "sparse_array.h"

namespace libpetey {
  namespace libsparse {

/************************    sc_sparse_t  ************************/

    //constructors:
    template <class index_t, class data_t, class vector_t>
    void sc_sparse_t<index_t, data_t, vector_t>::sc_sparse_t(char *filename) {
      this->sparse<index_t, data_t>::sparse(filename);
    }
    
    template <class index_t, class data_t, class vector_t>
    void sc_sparse_t<index_t, data_t, vector_t>::sc_sparse_t(sc_mat_t<index_t, data_t, vector_t> *other) {
      this->sparse<index_t, data_t>::sparse(other);
    }

    //perform a vector multiplication:
    template <class index_t, class data_t, class vector_t>
    void sc_sparse_t<index_t, data_t, vector_t>::vect_mult2(vector_t &cand, vector_t &result) {
      for (index_t i=0; i<m; i++) result[i]=0;
      for (long i=0; i<nel; i++) {
        result[matrix[i].i]+=matrix[i].value*cand[matrix[i].j];
      }
    }

    //perform left vector multiplication:
    template <class index_t, class data_t, class vector_t>
    void sc_sparse_t<index_t, data_t, vector_t>::left_mult2(vector_t &cor, vector_t &result) {
      for (index_t i=0; i<m; i++) result[i]=0;
      for (long i=0; i<nel; i++) {
        result[matrix[i].j]+=matrix[i].value*cor[matrix[i].i];
      }
    }

    template <class index_t, class data_t, class vector_t>
    vector_t *sc_sparse_t<index_t, data_t, vector_t>::vect_mult2(vector_t *cand) {
      vector_t *result;

      if (cand->size() != this->n) {
        fprintf(stderr, "sc_sparse_t::vect_mult2: vector size (%d) and matrix inner dimension(%d) must match\n", cand->size(), this->n);
        return NULL;
      }
 
      result=new vector_t(m);
      vect_mult2(*cand, *result);
      return result;
    }

    template <class index_t, class data_t, class vector_t>
    vector_t * sc_sparse_t<index_t, data_t, vector_t>::left_mult2(vector_t *cor) {
      vector_t *result;
      if (cor->size() != this->m) {
        fprintf(stderr, "sc_sparse_t::left_mult2: vector size (%d) and matrix inner dimension(%d) must match\n", cor->size(), this->m);
        return NULL;
      }
 
      result=new vector_t(n);
      left_mult2(*cor, *result);
      return result;
    }

    template <class index_t, class data_t, class vector_t>
    void sc_sparse_t<index_t, data_t, vector_t>::get_row2 (index_t i, vector_t &row) {
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

    template <class index_t, class data_t, class vector_t>
    vector_t *sc_sparse_t<index_t, data_t, vector_t>::get_row2(index_t i) {
      vector_t *row;
      row=new vector_t(n);
      get_row(i, *row);
      return row;
    }

    //destructor:
    template <class index_t, class data_t, class vector_t>
    void sc_sparse_t<index_t, data_t, vector_t>::~sc_sparse_t() {
    }

/************************    sc_sparse_array_t  ************************/

    //constructors:
    template <class index_t, class data_t, class vector_t>
    void sc_sparse_array_t<index_t, data_t, vector_t>::sc_sparse_array_t(char *filename) {
      this->sparse_array<index_t, data_t>::sparse_array(filename);
    }
    
    template <class index_t, class data_t, class vector_t>
    void sc_sparse_array_t<index_t, data_t, vector_t>::sc_sparse_array_t(sc_mat_t<index_t, data_t, vector_t> *other) {
      this->sparse_array<index_t, data_t, sparse_t>::sparse_array(other);
    }

    template <class index_t, class real, class vector_t>
    void sc_sparse_array_t<index_t, real, vector_t>::vect_mult2(vector_t &cand, vector_t &result) {
      vector_t *med, *s, *r;
      med=new vector_t(m);

      sparse_a[0]->vect_mult2(cand, result);
      r=&result;
      for (long i=1; i<nsparse; i++) {
        sparse_a[i]->vect_mult2(*r, *med);
        s=med;
        med=r;
        r=s;
      }
      //this check should probably be applied everywhere I use this trick
      //and result is not allocated inside the function...
      if (nsparse % 2 == 0) {
        for (index_t i=0; i<m; i++) result[i]=(*med)[i];
      }
      delete med;
    }

    template <class index_t, class real, class vector_t>
    vector_t * sc_sparse_array_t<index_t, real, vector_t>::vect_mult2(vector_t *cand) {
      vector_t *result;
      if (cand->size() != this->n) {
        fprintf(stderr, "sc_sparse_array_t::vect_mult2: vector size (%d) and matrix inner dimension(%d) must match\n", cand->size(), this->n);
        return NULL;
      }
 
      result=new vector_t(m);
      vect_mult2(*cand, *result);
      return result;
    }

    template <class index_t, class real, class vector_t>
    void sc_sparse_array_t<index_t, real, vector_t>::left_mult2(vector_t &cor, vector_t &result) {
      vector_t *med, *s, *r;
      med=new vector_t(m);

      sparse_a[nsparse-1]->left_mult2(cor, result);
      r=&result;
      for (long i=nsparse-2; i>=0; i--) {
        sparse_a[i]->left_mult2(*r, *med);
        s=med;
        med=r;
        r=s;
      }
      if (nsparse % 2 == 0) {
        for (index_t i=0; i<m; i++) result[i]=(*med)[i];
      }
      delete med;
    }

    template <class index_t, class real, class vector_t>
    vector_t * sc_sparse_array_t<index_t, real, vector_t>::left_mult2(vector_t *cor) {
      vector_t *result;
      if (cor->size() != this->m) {
        fprintf(stderr, "sc_sparse_t::left_mult2: vector size (%d) and matrix inner dimension(%d) must match\n", cor->size(), this->m);
        return NULL;
      }
 
      result=new vector_t(m);
      left_mult2(*cor, *result);
      return result;
    }

    template <class index_t, class real, class vector_t>
    void sc_sparse_array_t<index_t, real, vector_t>::get_row2(index_t i, vector_t &row) {
      vector_t s1(m);
      for (index_t k=0; k<m; k++) s1[k]=0;
      s1[i]=1;
      left_mult2(s1, row);
    }

    template <class index_t, class real, class vector_t>
    vector_t *sc_sparse_array_t<index_t, real, vector_t>::get_row2(index_t i) {
      vector_t *row;
      row=new vector_t(m);
      get_row(i, *row);
      return row;
    }

    template <class index_t, class real, class vector_t>
    full_matrix<index_t, real, vector_t> * sparse_array<index_t, real, vector_t>::cmult2(vector_t &cand) {
      full_matrix<index_t, real, vector_t> * result;
      real cand2[this->m];

      if (cand->size() != this->n) {
        fprintf(stderr, "sc_sparse_t::vect_mult2: vector size (%d) and matrix inner dimension(%d) must match\n", cand->size(), this->m);
        return NULL;
      }
 
      result=new full_matrix<index_t, real, vector_t>(nsparse, m);

      for (index_t i=0; i<m; i++) cand2[i]=cand[i];
      sparse_a[0]->vect_mult(cand2, result->data[0]);
      for (long i=1; i<nsparse; i++) {
        sparse_a[i]->vect_mult(result->data[i-1], result->data[i]);
      }
      return result;
    }

    template <class index_t, class real, class vector_t>
    full_matrix<index_t, real, vector_t> * sparse_array<index_t, real, vector_t>::left_cmult2(vector_t &cor) {
      real cor2=real[this->m];
      full_matrix<index_t, real, vector_t> * result;

      if (car->size() != this->n) {
        fprintf(stderr, "sc_full_t::left_cmult2: vector size (%d) and matrix inner dimension(%d) must match\n", cor->size(), this->m);
        return NULL;
      }
 
      for (index_t i=0; i<m; i++) cor2[i]=cor[i];
      result=new full_matrix<index_t, real, vector_t>(nsparse, m);

      sparse_a[nsparse-1]->left_mult(cor2, result->data[0]);
      for (long i=nsparse-2; i>=0; i++) {
        sparse_a[i]->left_mult(result->data[nsparse-i-2], result->data[nsparse-i-1]);
      }

      return result;
    }

    //destructor:
    template <class index_t, class data_t, class vector_t>
    void sc_sparse_array_t<index_t, data_t, vector_t>::~sc_sparse_array_t() {
    }

/************************    sc_full_t  ************************/

    //constructors:
    template <class index_t, class data_t, class vector_t>
    void sc_full_t<index_t, data_t, vector_t>::sc_full_t(char *filename) {
      this->full_matrix<index_t, data_t>::full_matrix(filename);
    }
    
    template <class index_t, class data_t, class vector_t>
    void sc_full_t<index_t, data_t, vector_t>::sc_full_t(sc_mat_t<index_t, data_t, vector_t> *other) {
      this->full_matrix<index_t, data_t>::full_matrix(other);
    }

    //multiply with a vector:
    template <class index_t, class scalar, class vector_t>
    void sc_full_t<index_t, scalar, vector_t>::vect_mult2(vector_t &cand, vector_t &result) {
      for (index_t i=0; i<m; i++) {
        result[i]=0;
        for (index_t j=0; j<n; j++) result[i]+=data[i][j]*cand[j];
      }
    }

    template <class index_t, class scalar, class vector_t>
    vector_t * sc_full_t<index_t, scalar, vector_t>::vect_mult2(vector_t *cand) {
      vector_t *result;

      if (cand->size() != this->n) {
        fprintf(stderr, "sc_full_t::vect_mult2: vector size (%d) and matrix inner dimension(%d) must match\n", cand->size(), this->n);
        return NULL;
      }
 
      result=new vector_t(m);
      vect_mult2(*cand, *result);
      return result;
    }

    template <class index_t, class scalar, class vector_t>
    void sc_full_t<index_t, scalar, vector_t>::left_mult2(vector_t &cor, vector_t &result) {
      for (index_t i=0; i<n; i++) result[i]=0;
      for (index_t i=0; i<m; i++) {
        for (index_t j=0; j<n; j++) {
          result[j]+=cor[i]*data[i][j];
        }
      }
    }

    template <class index_t, class scalar, class vector_t>
    vector_t * sc_full_t<index_t, scalar, vector_t>::left_mult2(vector_t *cor) {
      vector_t *result;
      if (cor->size() != this->m) {
        fprintf(stderr, "sc_full_t::left_mult2: vector size (%d) and matrix inner dimension(%d) must match\n", cor->size(), this->m);
        return NULL;
      }
 
      result=new vector_t(n);
      left_mult2(*cor, *result);
      return result;
    }

    template <class index_t, class scalar, class vector_t>
    void sc_full_t<index_t, scalar, vector_t>::get_row2(index_t i, vector_t &row) {
      //pretty simple:
      for (index_t j=0; j<n; j++) row[j]=data[i][j];
    }

    template <class index_t, class scalar, class vector_t>
    vector_t *sc_full_t<index_t, scalar, vector_t>::get_row2(index_t i) {
      vector_t *row;
      row=new vector_t(n);
      get_row(i, *row);
      return row;
    }

    //destructor:
    template <class index_t, class data_t, class vector_t>
    void sc_full_t<index_t, data_t, vector_t>::~sc_full_t() {
    }


  }
}

