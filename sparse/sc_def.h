#ifndef "_SC_DEF__H_"
#define "_SC_DEF__H_"

#include <vector>

#include "sc_type_base.h"
#include "sc_type_literal.h"
#include "symbol_table.h"

#define SCALAR 1
#define LITERAL 2
#define VECTOR 3
#define FULL 4
#define SPARSE 5
#define SPARSE_ARRAY 6
#define LIST 7

char *type_name[8]={"NULL", "scalar", "literal", "vector", "full", "sparse", "sparse_array", "list"};

typedef double sc_scalar_t;
typedef int sc_int_t;

typedef sc_type_base * (*) (sc_type_list *) sc_fun_t;

//input stream (file or stdin):
extern FILE *sc_istream;

//session id
extern unsigned long sc_session_id;

//function table:
extern symbol_table<sc_literal *> sc_funtab;
extern vector<sc_fun_t> sc_flist;

//variable table:
extern symbol_table<sc_literal *> sc_vartab;
extern vector<int> sc_vartype;
extern vector<sc_type_base *> sc_vlist;
extern vector<char> sc_delflag;

//look up a variable and return it:
sc_type_base * sc_var_lookup(sc_type_literal *name);
//create a new variable:
int sc_var_new(sc_type_literal *name, int typecode);
//assign an existing expression to a variable:
int sc_var_assign(sc_type_literal *name, sc_type_base *var);

//look up a function and return it:
sc_fun_t sc_fun_lookup(sc_type_literal *name);

//add an argument:
int sc_add_arg(char *symbol, int typecode);
//get an argument:
sc_type_base * sc_getarg(int index);

#endif

