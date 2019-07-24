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
      result=new sc_type_matrix(el1, FULL_T);
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
      result=new sc_type_matrix(el1, SPARSE_T);
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
      result=new sc_type_matrix(el1, SPARSE_ARRAY_T);
    } else {
      fprintf("Can't convert to sparse_array\n");
      result=NULL;
    }
  } else {
    fprintf("Wrong number of arguments to function, sparse_array\n");
    result=NULL;
  }
  return result;
}

sc_type_base *sc_builtin_list(sc_state_struct *state, sc_list *args) {
  sc_type_base *result;
  //check number and type of arguments:
  int nel=args->nel();

  if (nel==1) {
    sc_type_base *el1=args->get_el(0);
    if (sc_type_of(el1) == LITERAL_T) {
      state->vartab.add(el1);
      state->vartype[id]=LIST_T;
      result=new sc_type_list(el1);
    } else {
      fprintf("Can't convert to list: use cat\n");
      result=NULL;
    }
  } else {
    fprintf("Wrong number of arguments to function, list\n");
    result=NULL;
  }
  return result;
}

sc_type_base *sc_builtin_cat(sc_state_struct *state, sc_list *args) {
  //pretty simple:
  return args;
}

sc_type_base *sc_builtin_size(sc_state_struct *state, sc_list *args) {
  sc_type_base *result;
  //check number and type of arguments:
  int nel=args->nel();
  if (nel==1) {
    sc_type_base *el1=args->get_el(0);
    result=el1->size();
  } else {
    fprintf("Wrong number of arguments to function, size\n");
    result=NULL;
  }
  return result;
}

sc_type_base *sc_builtin_delete(sc_state_struct *state, sc_list *args) {
  sc_type_base *result;
  //check number and type of arguments:
  int nel=args->nel();
  if (nel==1) {
    sc_type_base *el1=args->get_el(0);
    if (sc_type_of(el1) == LITERAL_T) {
      long id=state->vartab.lookup(el1);
      if (id==-1) {
        fprintf(stderr, "Variable not found\n");
	result=NULL;
      } else {
        state->delflag[id]=1;
	result=new sc_type_scalar(0.);
      }
    } else {
      fprintf("delete: argument must be must variable\n");
      result=NULL;
    }
  } else {
    fprintf("Wrong number of arguments to function, delete\n");
    result=NULL;
  }
  return result;
}
  
sc_type_base *sc_builtin_print(sc_state_struct *state, sc_list *args) {
  sc_type_base *result;
  int err;
  //check number and type of arguments:
  int nel=args->nel();
  if (nel==1) {
    err=args->get_el(0)->print(stdout);
    result=new sc_type_scalar(err);
  } else if (nel==2) {
    sc_type_base *fname=args->get_el(1);
    if (fname->typeof() != LITERAL) {
      fprintf(stdout, "Error in print: second argument must be file name\n");
      result=NULL;
    } else {
      FILE *fs;
      char *name=(char *) fname;
      if (fname[0]=='>') {
        fs=fopen(name, "a");
      } else {
        fs=fopen(name, "w");
      }
      err=args->get_el(0)->print(fs);
      fclose(fs);
      result=new sc_type_scalar(err);
    }
  } else {
    fprintf("Wrong number of arguments to function, print\n");
    result=NULL;
  }

  return result;
}

sc_type_base *sc_builtin_save(sc_state_struct *state, sc_list *args) {
  sc_type_base *result;
  //check number and type of arguments:
  int nel=args->nel();
  if (nel==1) {
    sc_type_base *el1=args->get_el(0);
    if (sc_type_of(el1) == LITERAL_T) {
      long id=state->vartab.lookup(el1);
      if (id==-1) {
        fprintf(stderr, "Variable not found\n");
	result=NULL;
      } else {
        state->delflag[id]=0;
	result=new sc_type_scalar(0.);
      }
    } else {
      fprintf("save: argument must be variable name\n");
      result=NULL;
    }
  } else {
    fprintf("Wrong number of arguments to function, delete\n");
    result=NULL;
  }
  return result;
}
  
sc_type_base *sc_builtin_read(sc_state_struct *state, sc_list *args) {
  sc_type_base *result;
  int err;
  //check number and type of arguments:
  int nel=args->nel();
  if (nel==2) {
    sc_type_base *type0=args->get_el(0);
    sc_type_base *fname=args->get_el(1);
    if (type->typeof() != LITERAL_T || fname->typeof() != LITERAL_T) {
      fprintf("Wrong type(s) to builtin, read\n");
      result=NULL;
    } else {
      FILE *fs;
      char *name=(char *) fname;
      fs=fopen(name, "w");
      if (fs == NULL) {
        fprintf("Unable to open file, %s\n", name);
	result=NULL;
      }
      if (type0 == "scalar") {
        result=new sc_type_scalar();
      } else if (type0 == "literal") {
        result=new sc_type_literal();
      } else if (type0 == "vector") {
        result=new sc_type_vector();
      } else if (type0 == "full") {
        result=new sc_type_matrix(FULL_T);
      } else if (type0 == "sparse") {
        result=new sc_type_matrix(SPARSE_T);
      } else if (type0 == "sparse_array") {
        result=new sc_type_matrix(SPARSE_ARRAY_T);
      } else if (type0 == "list") {
        result=new sc_type_list();
      } else {
        fprintf("Unrecognized type name\n");
	result=NULL;
      }
      err=result->load(fs);
      fclose(fs);
      if (err!=0) {
        delete result;
	result=NULL;
      }
    }
  } else {
    fprintf("Wrong number of arguments to function, print\n");
    result=NULL;
  }

  return result;
}

sc_type_base *sc_builtin_empty(sc_state_struct *state, sc_list *arg) {
  sc_type_base *result;
  int nel=arg->nel();
  if (nel!=2) {
    fprintf("Wrong number of arguments to builtin, empty\n");
    result=NULL;
  } else {
    if (arg->get_el(0)->typeof() != SCALAR_T || arg->get_el(1)->typeof() != SCALAR_T) {
      fprintf("Wrong type of arguments to builtin, empty\n");
      result=NULL;
    } else {
      int m=(double) arg->get_el(0);
      int n=(double) arg->get_el(1);
      result=new sc_type_matrix(m, n);
    }
  }
  return result;
}
  
sc_type_base *sc_builtin_identity(sc_state_struct *state, sc_list *arg) {
  sc_type_base *result;
  int nel=arg->nel();
  if (nel!=2) {
    fprintf("Wrong number of arguments to builtin, identity\n");
    result=NULL;
  } else {
    if (arg->get_el(0)->typeof() != SCALAR_T || arg->get_el(1)->typeof() != SCALAR_T) {
      fprintf("Wrong type of arguments to builtin, identity\n");
      result=NULL;
    } else {
      int m=(double) arg->get_el(0);
      int n=(double) arg->get_el(1);
      result=new sc_type_matrix(m, n, 1);
    }
  }
  return result;
}
  
sc_type_base *sc_builtin_transpose(sc_state_struct *state, sc_list *arg) {
  sc_type_base *result;
  int nel=arg->nel();
  if (nel!=1) {
    fprintf("Wrong number of arguments to builtin, transpose\n");
    result=NULL;
  } else {
    result=arg->get_el(0)->transpose();
  }
  return result;
}
