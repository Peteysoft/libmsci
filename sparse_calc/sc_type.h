#ifndef SC_TYPE__H 
#define SC_TYPE__H

#include "sparse_array.h"

using namespace libpetey;
using namespace libsparse;

namespace sparse_calc {

typedef double scalar_t;
typedef int32_t integer;

typedef sparse<integer, scalar_t> sparse_t;
typedef full_matrix<integer, scalar_t> full_t;
typedef sparse_array<integer, scalar_t> sparse_array_t; 

//variable type as integer code:
int sc_type_of(sc_t *v);

class sc_t {
  public:
    virtual ~sc_t()=0;
    virtual sc_t * copy();
    int tflag;

    //universal operators:
    virtual sc_t * mult(sc_t *cand);		//multiplication
    virtual sc_t * add(sc_t *b);		//addition
    virtual void neg()=0;			//negation
    virtual sc_t * norm();			//norm

    //operators specific to certain types:
    virtual sc_t *sub(sc_t *sub);		//1-d subscript	
    						//matrix | vector + vector | scalar)
    virtual sc_t *sub(sc_t *sub, sc_t *sub);	//2-d subscript 
    						//matrix + (vector + vector) | (scalar + scalar) (doesn't allow mixed subscripts)
    virtual sc_t *cprod(sc_t *cand);		//running product (sparse-array + (matrix | vector)

    virtual int load(FILE *fs);
    virtual int save(FILE *fs);
    virtual int read(FILE *fs);
    virtual int print(FILE *fs);

};

class sc_str_t: public sc_t {
  public:
    sc_str_t(char *str);
    virtual ~sc_str_t();
    virtual sc_t *sub(sc_t *sub);
    char *s;

};

template <class scalar_t>
class sc_scal_t: public sc_t {
  public:
    sc_scal_t(scalar_t val);
    virtual ~sc_scal_t();
    scalar_t value;

    virtual sc_t * mult(sc_t *cand);
    virtual sc_t * add(sc_t *b);
    virtual void neg();
    virtual sc_t * norm();
    virtual sc_t *sub(sc_t *sub);

    virtual int load(FILE *fs);
    virtual int save(FILE *fs);
    virtual int read(FILE *fs);
    virtual int print(FILE *fs);
};

template <class scalar_t>
class sc_vec_t: public sc_t {
  public:
    scalar_t *data;
    sc_int_t n;

    sc_vec_t();
    sc_vec_t(scalar_t start, scalar_t finish);
    sc_vec_t(sc_int_t n);
    sc_vec_t(char *filename);
    virtual ~sc_vec_t();

    virtual sc_t * mult(sc_t *cand);
    virtual sc_t * add(sc_t *b);
    virtual void neg();
    virtual sc_t * norm();
    virtual sc_t *sub(sc_t *sub);

    virtual int load(FILE *fs);
    virtual int save(FILE *fs);
    virtual int read(FILE *fs);
    virtual int print(FILE *fs);
};

template <class scalar_t>
class sc_mat_t: public sc_t {
  public:
    matrix_base<sc_int_t, scalar_t> *mat;

    sc_mat_t();
    virtual ~sc_mat_t();

    virtual sc_t * mult(sc_t *cand);
    virtual sc_t * add(sc_t *b);
    virtual void neg();
    virtual sc_t * norm();
    //virtual sc_t *sub(sc_t *sub);

    virtual int load(FILE *fs);
    virtual int save(FILE *fs);
    virtual int read(FILE *fs);
    virtual int print(FILE *fs);
};


template <class scalar_t>
class sc_list_t: public sc_t {
  public:
    vector_e<sc_t *> list;

    sc_list_t();
    virtual ~sc_list_t();

    virtual sc_t * mult(sc_t *cand);
    virtual sc_t * add(sc_t *b);
    virtual void neg();
    virtual sc_t * norm();
    virtual sc_t *sub(sc_t *sub);

    //multiply through:
    sc_t * resolve();

    //left multiply:
    sc_t * left_mult(sc_t *cor);

    virtual int load(FILE *fs);
    virtual int save(FILE *fs);
    virtual int read(FILE *fs);
    virtual int print(FILE *fs);
};

}

#endif

