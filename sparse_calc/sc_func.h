#ifndef SC_FUNC__H
#define SC_FUNC__H

#include "sc_type.h"

//forward declarations:
template <typename type>
vector_e;

class sc_sub_t {
  public:
    virtual ~sc_sub_t()=0;
    virtual sc_t * call(int narg, sc_t **args)=0;

};

class sc_func_t:public sc_sub_t {
  public:
    sc_func_t(sc_t * (*fp) (int, sc_t **));
    ~sc_func_t();
    virtual sc_t * call(vector_e<sc_t> *stack);
};

class sc_user_t:public sc_sub_t {
  protected:
    //environment:
    symbol_table *name;		//variable names
    vector_s<int> *vtype;	//local variables
    vector_s<scalar_t> *scalar; //list of scalar variables 
    vector_s<long> *code;	//byte code

    //function table:
    symbol_table *com;
    vector_e<*sc_sub_t> *sub;
  public:
    sc_user_t(char *filename, symbol_table *c1, vector_e<*sc_sub_t> *s1);
    ~sc_user_t();
    virtual sc_t * call(vector_e<sc_t> *stack);

    //push variable and function, respectively onto byte code:
    int push_var(char *symbol);
    int push_func(char *symbol);

    //add a variable without adding to code base:
    int add_var(char *symbol, sc_t *v);
    //clear the code while leaving variables intact:
    void clear_code();

};

#endif

