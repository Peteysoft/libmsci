#ifndef "_SC_DEF__H_"
#define "_SC_DEF__H_"

#include <vector>

#include "sc_type_base.h"
#include "sc_type_literal.h"
#include "symbol_table.h"

#define SCALAR_T 1
#define LITERAL_T 2
#define VECTOR_T 3
#define FULL_T 4
#define SPARSE_T 5
#define SPARSE_ARRAY_T 6
#define LIST_T 7

char *type_name[8]={"NULL", "scalar", "literal", "vector", "full", "sparse", "sparse_array", "list"};

typedef double sc_scalar_t;
typedef int sc_int_t;

typedef sc_type_base * (*) (sc_state_struct *, sc_type_list *) sc_fun_t;

struct sc_state_struct {

  //input stream (file or stdin):
  FILE *istream;

  //session id
  unsigned long session_id;

  //data path:
  sc_literal *data_path;

  //code path:
  sc_literal *code_path;

  //function table:
  symbol_table<sc_literal *> funtab;
  vector<sc_fun_t> flist;

  //variable table:
  symbol_table<sc_literal *> vartab;
  vector<int> vartype;
  vector<sc_type_base *> vlist;
  vector<char> delflag;
} sc_state;

sc_state_init(sc_state_struct *state);
void sc_add_arg(sc_state_struct *state, sc_type_base *arg);
int sc_state_destroy(sc_state_struct *state);

//look up a variable and return it:
sc_type_base * sc_var_lookup(sc_state_struct *state, sc_type_literal *name);
//assign an existing expression to a variable:
int sc_var_assign(sc_state_struct *state, sc_type_literal *name, sc_type_base *var);

//look up a function and return it:
sc_fun_t sc_fun_lookup(sc_state_struct *state, sc_type_literal *name);

//variable type code:
int sc_type_of(sc_type_base *var);

//read in a variable from a file:
sc_type_base *sc_read_var(sc_literal *name, int type);


#endif

