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

      //adds a new element, returns number of elements:
      long add(type data);
      //adds a unique element (member), returns number of elements:
      long add_member(type data);
      //return number of elements:
      long nel();
      //decompose tree into a sorted array:
      void decompose(type *sarr, long nd);

      //delete smallest element:
      long delete_least();
      //delete largest element:
      long delete_greatest();

      //print out the tree:
      void print(FILE *fs);
  };

}

#endif

