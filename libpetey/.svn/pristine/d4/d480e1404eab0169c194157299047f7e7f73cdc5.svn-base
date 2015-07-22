
//reserved op codes:
#define SC_CLOBBER_CODE -1	//clobber assignment
#define SC_ELBYEL_CODE -2	//element-by-element assignment 
#define SC_DEC_VEC_CODE -3	//declare a vector
#define SC_DEC_MAT_CODE -4	//declare a matrix
#define SC_END_DEC_CODE -5	//end declaration section

struct sc_func_table {
  symbol_table funtab;
  vector_e<sc_sub_t *> funs;
};

struct sc_code_env {
  int instr_ptr;
  int errcode;
  symbol_table *varname;
  vector_e<sc_t *> *var;
  vector_s<int> *vtype;
  vector_e<sc_t *> *stack;
};

void code_push_constant(scalar c);
void code_push_literal(char *c);
void code_push_narg(int n);

sc_error sc_code_interp(sc_code_env *env, sc_func_table *funtab, vector_s<long> *code);

