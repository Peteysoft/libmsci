#include "sc_fun.h"

sc_type_base *sc_builtin_vector(sc_state_struct *state, sc_list *args) {
  sc_type_base *result;
  //check number and type of arguments:
  int nel=args->nel();

  if (nel==1) {
    sc_type_base *el1=args->get_el(0);
    if (sc_type_of(el1) == LITERAL_T) {
      state->vartab.add(el1);
      state->vartype[id]=VECTOR_T;
      result=new sc_type_vector(el1);
    } else {
      result=el1->vector();
    }
  } else {
    fprintf("Wrong number of arguments to function, vector\n");
    result=NULL;
  }
  return result;
}

sc_type_base *sc_builtin_full(sc_state_struct *state, sc_list *args) {
  sc_type_base *result;
  //check number and type of arguments:
  int nel=args->nel();

  if (nel==1) {
    sc_type_base *el1=args->get_el(0);
    if (sc_type_of(el1) == LITERAL_T) {
      state->vartab.add(el1);
      state->vartype[id]=FULL_T;
      result=new sc_type_full(el1);
    } else {
      result=el1->full();
    }
  } else {
    fprintf("Wrong number of arguments to function, full\n");
    result=NULL;
  }
  return result;
}

sc_type_base *sc_builtin_sparse(sc_state_struct *state, sc_list *args) {
  sc_type_base *result;
  //check number and type of arguments:
  int nel=args->nel();

  if (nel==1) {
    sc_type_base *el1=args->get_el(0);
    if (sc_type_of(el1) == LITERAL_T) {
      state->vartab.add(el1);
      state->vartype[id]=SPARSE_T;
      result=new sc_type_sparse(el1);
    } else {
      result=el1->sparse();
    }
  } else {
    fprintf("Wrong number of arguments to function, sparse\n");
    result=NULL;
  }
  return result;
}

sc_type_base *sc_builtin_sparse_array(sc_state_struct *state, sc_list *args) {
  sc_type_base *result;
  //check number and type of arguments:
  int nel=args->nel();

  if (nel==1) {
    sc_type_base *el1=args->get_el(0);
    if (sc_type_of(el1) == LITERAL_T) {
      state->vartab.add(el1);
      state->vartype[id]=SPARSE_ARRAY_T;
      result=new sc_type_sparse_array(el1);
    } else {
      fprintf("Can't convert to sparse_array\n");
      result=NULL;
    }
  } else {
    fprintf("Wrong number of arguments to function, sparse\n");
    result=NULL;
  }
  return result;
}

