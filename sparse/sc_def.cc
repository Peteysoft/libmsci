
sc_state_init(sc_state_struct *state) {
  state->data_path=new sc_literal(".");
  state->code_path=new sc_literal(".");

  //I/O:
  int id=state->funtab.add("print");
  flist[id]=& sc_builtin_print;
  id=state->funtab.add("save");
  flist[id]=& sc_builtin_save;
  id=state->funtab->add("read");
  flist[id]=& sc_builtin_read;

  //type definition/conversion:
  id=state->funtab->add("vector");
  flist[id]=& sc_builtin_vector;
  id=state->funtab->add("full");
  flist[id]=& sc_builtin_full;
  id=state->funtab->add("sparse");
  flist[id]=& sc_builtin_sparse;
  id=state->funtab->add("sparse_array");
  flist[id]=& sc_builtin_sparse_array;
  id=state->funtab.add("list");
  flist[id]=& sc_builtin_list;
  id=state->funtab.add("cat");
  flist[id]=& sc_builtin_cat;
  id=state->funtab.add("size");
  flist[id]=& sc_builtin_size;
  id=state->funtab.add("delete");
  flist[id]=& sc_builtin_delete;

  //math:
  id=state->funtab.add("empty");
  flist[id]=& sc_builtin_empty;
  id=state->funtab.add("identity");
  flist[id]=& sc_builtin_identity;
  id=state->funtab.add("transpose");
  flist[id]=& sc_builtin_transpose;
  id=state->funtab.add("solve");
  flist[id]=& sc_builtin_solve;
  id=state->funtab.add("evd");
  flist[id]=& sc_builtin_evd;
  id=state->funtab.add("svd");
  flist[id]=& sc_builtin_svd;

  //system:
  id=state->funtab.add("sys");
  flist[id]=& sc_builtin_sys;
  id=state->funtab.add("com");
  flist[id]=& sc_builtin_com;
  id=state->funtab.add("path");
  flist[id]=& sc_builtin_path;
  id=state->funtab.add("run");
  flist[id]=& sc_builtin_run;
  id=state->funtab.add("help");
  flist[id]=& sc_builtin_help;

}

void sc_add_arg(sc_state_struct *state, sc_type_base *arg) {
  long id=state->vartab.add("$_");
  vartype[id]=LIST_T;
  vlist[id]=NULL;
  delflag[id]=0;
}

sc_type_base * sc_read_var(sc_literal *name, int type) {
  sc_type_base *var;
  switch (type) {
    VECTOR_T:
      var=new sc_type_vector(name);
      break;
    FULL_T:
      var=new sc_type_full(name);
      break;
    SPARSE_T:
      var=new sc_type_sparse(name);
      break;
    SPARSE_ARRAY_T:
      var=new sc_type_sparse_array(name);
      break;
    LIST_T:
      var=new sc_type_list(name);
      break;
    default:
      fprintf("Error: variable type not stored on drive\n");
      var=NULL;
  }
  return var;
}

sc_type_base * sc_var_lookup(sc_state_struct *state, sc_literal *name) {
  sc_type_base *var;

  long id=state->vartab.lookup(name);
  if (id == -1) {
    var=NULL;
  } else if (state->vartype[id] == SCALAR_T || state->vartype[id] == LITERAL_T) {
    var=vlist[id];
  } else {
    var=sc_read_var(name, state->vartype[id]);
  }

  return var;
}

int sc_var_assign(sc_state_struc *state, sc_type_literal *name, sc_type_base *var) {
  int type=sc_type_of(var);
  long id=state->vartab.lookup(name);
  if (id == -1) {
    id=state->vartab.lookup(name);
  }
  state->vartype[id]=type;
  if (type == SCALAR || type == LITERAL) {
    state->vlist[id]=var;
  } else {
    FILE *fs;
    char *nm;

    state->vlist[id]=NULL;
    nm=(char *) name;
    fs=fopen(nm, "r");
    if (fs == NULL) {
      fprintf(stderr, "Error opening file, %s\n", nm);
      return ERROR_OPENING_FILE;
    }
    //(what's the exit-code?):
    var->write(fs);
    state->delflag[id]=1;
  }
}


sc_fun_t sc_fun_lookup(sc_state_struct *state, sc_literal *name) {
  sc_type_base * (*fun) (sc_state_struct *, sc_list *args);
  long id=state->funtab.lookup(name);
  if (id == -1) {
    fun=NULL;
  } else {
    fun=state->flist[id];
  }
  return fun;
}

int sc_type_of(sc_type_base *var) {
  int type;
  if (typeid(*var) == typeid(sc_scalar)) {
    type=SCALAR_T;
  } else if (typeid(*var) == typeid(sc_literal)) {
    type=LITERAL_T;
  } else if (typeid(*var) == typeid(sc_type_full)) {
    type=FULL_T;
  } else if (typeid(*var) == typeid(sc_sparse)) {
    type=SPARSE_T;
  } else if (typeid(*var) == typeid(sc_sparse_array)) {
    type=SPARSE_ARRAY_T;
  } else if (typeid(*var) == typeid(sc_list)) {
    type=LIST_T;
  } else {
    type=-1;
  }
  return type;
}

