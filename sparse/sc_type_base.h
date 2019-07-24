#include "sparse_def.h"

class sc_type_base {
  protected:
    int type_code;
  public:
    sc_type_base();
    ~sc_type_base();

    //Part I: functions that don't return another sc type:

    //I/O:
    virtual int load(FILE *fs);
    virtual int write(FILE *fs);

    virtual int scan(FILE *fs);
    virtual int print(FILE *fs);

    //informational:
    virtual int nel();
    virtual int typeof();

    //type conversion:
    virtual char * operator();
    virtual double operator();


    //Part II: returns another sc type:

    //informational:
    virtual sc_type_base *size();

    //modifies this class:
    virtual sc_type_base * negate();		//works for all
    virtual sc_type_base * transpose();		//only works for matrices
    
    //doesn't modify this class:
    //works for all:
    virtual sc_type_base *norm();
    //may fail if dimensions don't match and for certain type combinations:
    virtual sc_type_base *add(sc_type_base *other);
    virtual sc_type_base *multiply(sc_type_base *other);
    //only works for lists and sparse arrays:
    virtual sc_type_base *multiply_all(sc_type_base *other);

    //(need to list all the operators...)
    virtual sc_type_base *divide(sc_type_base *other);
    virtual sc_type_base *mod(sc_type_base *other);
    virtual sc_type_base *pow(sc_type_base *other);

    //works for scalars and literals:
    virtual sc_type_base *gt(sc_type_base *other);
    virtual sc_type_base *lt(sc_type_base *other);
    virtual sc_type_base *ge(sc_type_base *other);
    virtual sc_type_base *le(sc_type_base *other);
    virtual sc_type_base *eq(sc_type_base *other);
    virtual sc_type_base *ne(sc_type_base *other);

    //only works for literals:
    virtual sc_type_base *cat(sc_type_base *other);

    //vectors, lists and sparse_arrays:
    virtual sc_type_base *map(sc_type_literal *op);
    virtual sc_type_base *map(sc_type_literal *op, sc_type_base *other);

    //subscripting:
    //for matrices and lists:
    virtual sc_type_base *subscript(sc_type_base *sub);
    virtual sc_type_base *subscript(sc_type_base *sub1, sc_type_base *sub2);
    //doesn't work for sparse arrays:
    virtual sc_type_base *subscript_assign(sc_type_base *sub);
    virtual sc_type_base *subscript_assign(sc_type_base *sub1, sc_type_base *sub2);

    //for lists and sparse arrays:
    virtual sc_type_base *get_el(sc_type_base *el);

};




