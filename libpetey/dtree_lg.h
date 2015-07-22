#ifndef DTREE_LG_H_INCLUDED
#define DTREE_LG_H_INCLUDED 1

namespace libpetey {

  template <class type>
  struct dtree_lg_el {

    type value;
    dtree_lg_el *left;
    dtree_lg_el *right;
  };

  template<class type>
  class dtree_lg {
    protected:
    
      dtree_lg_el<type> * trunk;
      long n;

      void delete_el(dtree_lg_el<type> *tel);
      void decompose(dtree_lg_el<type> *t, type *sarr, long nd, long &iter);

    public:
      dtree_lg();
      dtree_lg(type data, long nd);
      ~dtree_lg();

      long add(type data);
      long add_member(type data);
      long nel();
      void decompose(type *sarr, long nd);

      long delete_least();
      long delete_greatest();
  };
}

#endif

