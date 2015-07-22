#include <stdlib.h>
#include "tree_lgi.h"

namespace libpetey {

template<class type>
tree_lgi<type>::tree_lgi() {
  trunk=NULL;
  n=0;
}

template<class type>
tree_lgi<type>::~tree_lgi() {
  delete_el(trunk);
}


template<class type>
void tree_lgi<type>::delete_el(tree_lgi_el<type> *tel) {
  if (tel == NULL) return;
  delete_el(tel->right);
  delete_el(tel->left);
  delete tel;
}

template<class type>
void tree_lgi<type>::decompose(tree_lgi_el<type> *t, type *sarr, long *ind, long nd, long &iter) {
  if (t == NULL) return;
  if (iter >= nd) return;
  decompose(t->left, sarr, ind, nd, iter);
  sarr[iter]=t->value;
  ind[iter]=t->ind;
  iter++;
  decompose(t->right, sarr, ind, nd, iter);

}


template<class type>
long tree_lgi<type>::add(type data, long ind) {
  tree_lgi_el<type> *current;

  if (trunk==NULL) {
    trunk=new tree_lgi_el<type>;
    trunk->value=data;
    trunk->left=NULL;
    trunk->right=NULL;
    trunk->ind=ind;
    n=1;
    return n;
  }

  current=trunk;
  for (;;) {
    if (data < current->value) {
      if (current->left == NULL) {
        current->left=new tree_lgi_el<type>;
	current=current->left;
	break;
      } else {
        current=current->left;
      }
    } else {
      if (current->right == NULL) {
        current->right=new tree_lgi_el<type>;
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
  current->ind=ind;

  n++;

  return n;
}

template<class type>
long tree_lgi<type>::nel() {
  return n;
}

template<class type>
void tree_lgi<type>::decompose(type *sarr, long *ind, long nd) {
  long i=0;
  decompose(trunk, sarr, ind, nd, i);
}


template<class type>
long tree_lgi<type>::delete_least() {
  tree_lgi_el<type> *current, *previous;

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
long tree_lgi<type>::delete_greatest() {
  tree_lgi_el<type> *current, *previous;

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


template class tree_lgi<float>;
template class tree_lgi<double>;

} //end namespace libpetey

