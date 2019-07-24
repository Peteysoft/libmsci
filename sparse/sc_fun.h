#ifndef _SC_FUN_H_
#define _SC_FUN_H_

  #include "sc_type_base.h"
  #include "sc_def.h"

  //type definition/conversion:
  sc_type_base *sc_builtin_vector(sc_state_struct *state, sc_list * args);
  sc_type_base *sc_builtin_full(sc_state_struct *state, sc_list * args);
  sc_type_base *sc_builtin_sparse(sc_state_struct *state, sc_list * args);
  sc_type_base *sc_builtin_sparse_array(sc_state_struct *state, sc_list * args);
  sc_type_base *sc_builtin_list(sc_state_struct *state, sc_list *arg);
  sc_type_base *sc_builtin_cat(sc_state_struct *state, sc_list *arg);
  sc_type_base *sc_builtin_size(sc_state_struct *state, sc_list *arg);
  sc_type_base *sc_builtin_delete(sc_state_struct *state, sc_list *arg);

  //I/O:
  sc_type_base *sc_builtin_print(sc_state_struct *state, sc_list * args);
  sc_type_base *sc_builtin_save(sc_state_struct *state, sc_list * args);
  sc_type_base *sc_builtin_read(sc_state_struct *state, sc_list * args);

  //math:
  sc_type_base *sc_builtin_empty(sc_state_struct *state, sc_list *arg);
  sc_type_base *sc_builtin_identity(sc_state_struct *state, sc_list *arg);
  sc_type_base *sc_builtin_transpose(sc_state_struct *state, sc_list *arg);
  sc_type_base *sc_builtin_solve(sc_state_struct *state, sc_list *arg);
  sc_type_base *sc_builtin_evd(sc_state_struct *state, sc_list *arg);
  sc_type_base *sc_builtin_svd(sc_state_struct *state, sc_list *arg);
  sc_type_base *sc_builtin_map(sc_state_struct *state, sc_list *arg);

  //system:
  sc_type_base *sc_builtin_com(sc_state_struct *state, sc_list *arg);
  sc_type_base *sc_builtin_sys(sc_state_struct *state, sc_list *arg);
  sc_type_base *sc_builtin_run(sc_state_struct *state, sc_list *arg);
  sc_type_base *sc_builtin_path(sc_state_struct *state, sc_list *arg);
  sc_type_base *sc_builtin_help(sc_state_struct *state, sc_list *arg);

#endif

