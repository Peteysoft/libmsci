#include "sparse_def.h"

class sc_type_base {
  protected:
    int type_code;
  public:
    sc_type_base();
    ~sc_type_base();

    //informational:
    virtual sc_type_base *size();
    virtual sc_type_base *type();

    //modifies this class:
    virtual sc_type_base * negate();		//works for all
    virtual sc_type_base * transpose();		//only works for matrices
    
    //doesn't modify this class:
    //works for all:
    virtual sc_type_base *norm();
    //may fail if dimensions don't match:
    virtual sc_type_base *add(sc_type_base *other);
    virtual sc_type_base *multiply(sc_type_base *other);

    //only works for lists and sparse arrays:
    virtual sc_type_base *multiply_all(sc_type_base *other);

    //only works for scalars:
    virtual sc_type_base *scalar_op(sc_type_literal *op, sc_type_base *other);

    //only works for literals:
    virtual sc_type_base *cat(sc_type_base *other);

    //only works for vectors (?):
    virtual sc_type_base *distribute(sc_type_literal *op);
    virtual sc_type_base *distribute(sc_type_literal *op, sc_type_base *other);

    //I/O:
    virtual int read(FILE *fs);
    virtual int write(me);

    virtual int scan(sc_type_literal *fname);
    virtual int print(sc_type_literal *fname);

};




