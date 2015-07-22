#ifndef TREE_TMP_H_INCLUDED
#define TREE_TMP_H_INCLUDED 1

#pragma interface

namespace libpetey {

  //recursive data type which holds a binary tree:
  //implemented as a template
  template <class primitive>
  class tree {
    protected:
    //    data_type type;
      primitive *value;
      tree *left;
      tree *right;
      long index;		//exists to facilitate tracking indices
    public:
      tree();
      long add_member(primitive new_val, long new_ind);
      long add(primitive new_val, long new_ind);
      void decompose(primitive *data_arr, long *index_arr, long n, long &i);
      ~tree();
  };

  template <class primitive>
  primitive *treesort(primitive *data, long n, long **index);

}

#endif

