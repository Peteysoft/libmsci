#ifndef SC_OP__H
#define SC_OP__H

#include "sc_type.h"

//op codes:
extern long sc_neg_code;
extern long sc_sub_code;
extern long sc_norm_code;
extern long sc_assign_code;

//universal operators:
sc_t * sc_neg(int narg, sc_t **v);
sc_t * sc_prod(int narg, sc_t **v);
sc_t * sc_plus(int narg, sc_t **v);
sc_t * sc_minus(int narg, sc_t **v);
sc_t * sc_norm(int narg, sc_t **v);

//operators on non-scalars:
sc_t * sc_sub(int narg, sc_t **v);
sc_t * sc_range(int narg, sc_t **v);

//array and list operators:
sc_t * sc_cprod(int narg, sc_t **v);

//scalar operators:
sc_t * sc_div(int narg, sc_t **v);
sc_t * sc_pow(int narg, sc_t **v);
sc_t * sc_mod(int narg, sc_t **v);

//comparators (scalar):
sc_t * sc_gt(int narg, sc_t **v);
sc_t * sc_ge(int narg, sc_t **v);
sc_t * sc_lt(int narg, sc_t **v);
sc_t * sc_le(int narg, sc_t **v);
sc_t * sc_eq(int narg, sc_t **v);
sc_t * sc_ne(int narg, sc_t **v);

sc_t * sc_assign(int narg, sc_t **v);

//adds all the operators and built-ins to the function table:
long sc_op_init(sc_func_table *funtab);

#endif

