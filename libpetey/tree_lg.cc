#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#include <vector>

#include "tree_lg.h"

using namespace std;

namespace libpetey {

template<class type>
tree_lg<type>::tree_lg() {
  trunk=NULL;
  n=0;
  set_fcode();
}

template<class type>
tree_lg<type>::~tree_lg() {
  delete [] fcode;
  delete_el(trunk);
}


template<class type>
void tree_lg<type>::delete_el(tree_lg_el<type> *tel) {
  if (tel == NULL) return;
  delete_el(tel->right);
  delete_el(tel->left);
  delete tel;
}

template<class type>
void tree_lg<type>::decompose(tree_lg_el<type> *t, type *sarr, long nd, long &iter) {
  if (t == NULL) return;
  if (iter >= nd) return;
  decompose(t->left, sarr, nd, iter);
  sarr[iter]=t->value;
  iter++;
  decompose(t->right, sarr, nd, iter);
  
}


template<class type>
long tree_lg<type>::add(type data) {
  tree_lg_el<type> *current;
  
  if (trunk==NULL) {
    trunk=new tree_lg_el<type>;
    trunk->value=data;
    trunk->left=NULL;
    trunk->right=NULL;
    n=1;
    return n;
  }

  current=trunk;
  for (;;) {
    if (data < current->value) {
      if (current->left == NULL) {
        current->left=new tree_lg_el<type>;
	current=current->left;
	break;
      } else {
        current=current->left;
      }
    } else {
      if (current->right == NULL) {
        current->right=new tree_lg_el<type>;
	current=current->right;
	break;
      } else {
        current=current->right;
      }
    }
  }

  current->value=data;
  current->right=NULL;
  current->left=NULL;

  n++;

  return n;
}


//searches for an element
//if found, returns an error code
//if not, it's added to the tree...
template<class type>
long tree_lg<type>::add_member(type data) {
  tree_lg_el<type> *current;
  
  if (trunk==NULL) {
    trunk=new tree_lg_el<type>;
    trunk->value=data;
    trunk->left=NULL;
    trunk->right=NULL;
    n=1;
    return n;
  }

  current=trunk;
  for (;;) {
    if (data == current->value) {
      return -n;
    } else if (data < current->value) {
      if (current->left == NULL) {
        current->left=new tree_lg_el<type>;
	current=current->left;
	break;
      } else {
        current=current->left;
      }
    } else {
      if (current->right == NULL) {
        current->right=new tree_lg_el<type>;
	current=current->right;
	break;
      } else {
        current=current->right;
      }
    }
  }

  current->value=data;
  current->right=NULL;
  current->left=NULL;

  n++;

  return n;
}

template<class type>
long tree_lg<type>::nel() {
  return n;
}

template<class type>
void tree_lg<type>::decompose(type *sarr, long nd) {
  long i=0;
  decompose(trunk, sarr, nd, i);
}

template<class type>
long tree_lg<type>::delete_least() {
  tree_lg_el<type> *current, *previous;

  if (trunk == NULL) return n;

  if (trunk->left == NULL) {
    if (trunk->right == NULL) {
      delete trunk;
      trunk=NULL;
    } else {
      current=trunk->right;
      delete trunk;
      trunk=current;
    }
  } else {
    current=trunk;
    while (current->left != NULL) {
      previous=current;
      current=current->left;
    }
    previous->left=current->right;
    delete current;

  }
  n--;
  return n;
}

template<class type>
long tree_lg<type>::delete_greatest() {
  tree_lg_el<type> *current, *previous;

  if (trunk == NULL) return n;

  if (trunk->right == NULL) {
    if (trunk->left == NULL) {
      delete trunk;
      trunk=NULL;
    } else {
      current=trunk->left;
      delete trunk;
      trunk=current;
    }
  } else {
    current=trunk;
    while (current->right != NULL) {
      previous=current;
      current=current->right;
    }
    previous->right=current->left;
    delete current;

  }
  n--;
  return n;
}

template <class type>
void tree_lg<type>::set_fcode() {
  fcode=new char[1];
  strcpy(fcode, "");
}

template <>
void tree_lg<float>::set_fcode() {
  fcode=new char[3];
  strcpy(fcode, "%g");
}

template <>
void tree_lg<double>::set_fcode() {
  fcode=new char[4];
  strcpy(fcode, "%lg");
}

template <>
void tree_lg<int64_t>::set_fcode() {
  fcode=new char[4];
  strcpy(fcode, "%ld");
}

template <class type>
void tree_lg<type>::print(FILE *fs, tree_lg_el<type> *t, long depth) {
  if (t==NULL) return;
  print(fs, t->left, depth+1);
  for (long i=0; i<depth; i++) fprintf(fs, "  ");
  fprintf(fs, fcode, t->value);
  fprintf(fs, "\n");
  print(fs, t->right, depth+1);
}

template <>
void tree_lg<vector<int> >::print(FILE *fs, tree_lg_el<vector<int> > *t, long depth) {
  if (t==NULL) return;
  print(fs, t->left, depth+1);
  for (int i=0; i<t->value.size(); i++) {
    for (long j=0; j<depth; j++) fprintf(fs, "  ");
    fprintf(fs, "%d ", t->value[i]);
  }
  //fprintf(fs, "\n");
  fprintf(fs, "\n");
  print(fs, t->right, depth+1);
}

template <class type>
void tree_lg<type>::print(FILE *fs) {
  print(fs, trunk, 0);
}

template class tree_lg<float>;
template class tree_lg<double>;
template class tree_lg<int64_t>;
template class tree_lg<vector<int> >;

} //end namespace libpetey
