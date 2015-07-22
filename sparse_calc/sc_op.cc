#include "sc_glob.h"
#include "sc_op.h"

sc_t *sc_neg(int narg, sc_t **v) {
  int t;
  sc_t *result=NULL;

  if (narg != 1) return NULL;

  int t=sc_type_of(v[0]);

  if (t==SC_SCALAR_T || t==SC_VECTOR_T || t==SC_FULL_T || t==SC_SPARSE_T 
		  || t==SC_SPARSE_ARRAY || t==SC_STR_T) {
    result=v[0]->copy();
    result->neg();
  } else {
    fprintf(stderr, "Type mismatch in negation operator (-).\n");
  } 
  return result;
}

sc_t *sc_prod(int narg, sc_t **v) {
  sc_t *result=NULL;

  if (narg != 2) return NULL;

  result=v[0]->mult(v[1]);
  if (result==NULL) {
    fprintf(stderr, "Type mismatch in product operator (*).\n");
  }

  return result;
}

sc_t *sc_add(int narg, sc_t **v) {
  sc_t *result=NULL;

  if (narg != 2) return NULL;

  result=v[0]->add(v[1]);
  if (result==NULL) {
    fprintf(stderr, "Type mismatch in addition operator (+).\n");
  }

  return result;
}

sc_t *sc_minus(int narg, sc_t **v) {
  sc_t *result=NULL;

  if (narg != 2) return NULL;

  v[1].neg();
  result=v[0]->add(v[1]);
  if (result==NULL) {
    fprintf(stderr, "Type mismatch in subtract operator (-).\n");
  }

  return result;
}

sc_t *sc_div(int narg, sc_t **v) {
  sc_t *result=NULL;
  int t1, t2;

  if (narg!=2) return NULL;

  t1=sc_type_of(v[0]);
  t2=sc_type_of(v[1]);

  if (t1==SC_SCALAR_T && t2==SC_SCALAR_T) {
    sc_scal_t *r1;
    sc_scal_t *v1=(sc_scal_t *) v[0];
    sc_scal_t *v2=(sc_scal_t *) v[1];
    r1=new sc_scal_t(v1->val/v2->val);
    result=r1;
  } else {
    fprintf(stderr, "Type mismatch in division operator (/).\n");
  }

  return result;
}
  
sc_t *sc_sub(int narg, sc_t **v) {
  int t1, t2, t3;
  sc_t *result=NULL;

  if (narg==2) {
    result=v[0]->sub(v[1]);
  } else if (narg==3) {
    result=v[0]->sub(v[1], v[2]);
  }
  if (result == NULL) {
    fprintf(stderr, "Type mismatch in subscript operator ([]).\n");
  }

  return result;
}

sc_t *sc_range(int narg, sc_t **v) {
  int t1, t2;

  if (narg==2) {
    sc_scal_t *v1;
    sc_scal_t *v2;
    sc_vec_t *r1;
    t1=sc_type_of(v[0]);
    t2=sc_type_of(v[1]);
    if (t1 != SC_SCALAR_T && t2 != SC_SCALAR_T) {
      fprintf(stderr, "Type mismatch in range operator (..).\n");
      return NULL;
    }
    v1=(sc_scal_t *) v[0];
    v2=(sc_scal_t *) v[1];
    r1=new sc_vec_t(v1->value, v2->value);
  } else {
    return NULL;
  }

  return r1;
}

sc_t *sc_cprod(int narg, sc_t **v) {
  int t1, t2;
  sc_t *result=NULL;

  if (narg==2) {
    result=v[0]->cprod(v[1]);
  } else {
    assert(0);
  }
  if (result == NULL) {
    fprintf(stderr, "Type mismatch in cumulative product operator (#).\n");
  }

  return result;
}

sc_t *sc_norm(int narg, sc_t **v) {
  int t1, t2;
  sc_t *result=NULL;

  if (narg==1) {
    result=v[0]->norm();
  } else {
    return NULL;
  }
  if (result == NULL) {
    fprintf(stderr, "Type mismatch in norm operator (||).\n");
  }

  return result;
}

sc_t *sc_pow(int narg, sc_t **v) {
  sc_t *result=NULL;
  int t1, t2;

  if (narg!=2) return NULL;

  t1=sc_type_of(v[0]);
  t2=sc_type_of(v[1]);

  if (t1==SC_SCALAR_T && t2==SC_SCALAR_T) {
    sc_vec_t *r1;
    sc_scal_t *v1=(sc_scal_t *) v[0];
    sc_scal_t *v2=(sc_scal_t *) v[1];
    r1=new sc_scal_t(pow(v1->val, v2->val));
    result=r1;
  } else {
    fprintf(stderr, "Type mismatch in exponentiation operator (^).\n");
  }

  return result;
}
  
sc_t *sc_gt(int narg, sc_t **v) {
  sc_t *result=NULL;
  int t1, t2;

  if (narg!=2) return NULL;

  t1=sc_type_of(v[0]);
  t2=sc_type_of(v[1]);

  if (t1==SC_SCALAR_T && t2==SC_SCALAR_T) {
    sc_vec_t *r1;
    sc_scal_t *v1=(sc_scal_t *) v[0];
    sc_scal_t *v2=(sc_scal_t *) v[1];
    r1=new sc_scal_t(v1->val > v2->val);
    result=r1;
  } else {
    fprintf(stderr, "Type mismatch in greater-than operator (>).\n");
  }

  return result;
}

long sc_op_init(sc_func_table *funtab) {
  long id;
  id=funtab->funtab->add("\-");
  funtab->funs[id]=new sc_func_t(&sc_neg);

  id=funtab->funtab.add("*");
  funtab->funs[id]=new sc_func_t(&sc_prod);

  id=funtab->funtab.add("+");
  funtab->funs[id]=new sc_func_t(&sc_plus);

  id=funtab->funtab.add("-");
  funtab->funs[id]=new sc_func_t(&sc_minus);

  id=funtab->funtab.add("[]");
  funtab->funs[id]=new sc_func_t(&sc_sub);

  id=funtab->funtab.add("..");
  funtab->funs[id]=new sc_func_t(&sc_range);

  id=funtab->funtab.add("#");
  funtab->funs[id]=new sc_func_t(&sc_cprod);

  id=funtab->funtab.add("/");
  funtab->funs[id]=new sc_func_t(&sc_div);

  id=funtab->funtab.add("^");
  funtab->funs[id]=new sc_func_t(&sc_pow);

  return id;
}

