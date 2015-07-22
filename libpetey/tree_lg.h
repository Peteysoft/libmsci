#ifndef TREE_LG_H_INCLUDED
#define TREE_LG_H_INCLUDED 1

#include <stdio.h>

namespace libpetey {

  template <class type>
  struct tree_lg_el {

    type value;
    tree_lg_el *left;
    tree_lg_el *right;
  };

  template<class type>
  class tree_lg {
    protected:
    
      tree_lg_el<type> * trunk;
      long n;
      char *fcode;		//format code

      void delete_el(tree_lg_el<type> *tel);
      void decompose(tree_lg_el<type> *t, type *sarr, long nd, long &iter);

      void set_fcode();
      void print(FILE *fs, tree_lg_el<type> *t, long depth);

    public:
      tree_lg();
      tree_lg(type data, long nd);
      ~tree_lg();

      long add(type data);
      long add_member(type data);
      long nel();
      void decompose(type *sarr, long nd);

      long delete_least();
      long delete_greatest();

      void print(FILE *fs);
  };

}

#endif

