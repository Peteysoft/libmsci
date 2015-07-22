#ifndef TREE_LGI_H_INCLUDED
#define TREE_LGI_H_INCLUDED 1

namespace libpetey {

  template <class type>
  struct tree_lgi_el {

    type value;
    tree_lgi_el *left;
    tree_lgi_el *right;
    long ind;
  };

  template<class type>
  class tree_lgi {
    protected:

      tree_lgi_el<type> * trunk;
      long n;

      void delete_el(tree_lgi_el<type> *tel);
      void decompose(tree_lgi_el<type> *t, type *sarr, long *ind, long nd, long &iter);

    public:
      tree_lgi();
      tree_lgi(type data, long nd);
      ~tree_lgi();

      long add(type data, long ind);
      long add_member(type data, long ind);
      long nel();
      void decompose(type *sarr, long *ind, long nd);

      long delete_least();
      long delete_greatest();
  };

}

#endif

