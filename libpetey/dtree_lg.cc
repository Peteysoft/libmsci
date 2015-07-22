#include <stdlib.h>
#include "string_operators.h"
#include "dtree_lg.h"

namespace libpetey {

template<class type>
dtree_lg<type>::dtree_lg() {
  trunk=NULL;
  greatest=NULL;
  least=NULL;
  n=0;
}

template<class type>
dtree_lg<type>::~dtree_lg() {
  delete_el(trunk);
}


template<class type>
void dtree_lg<type>::delete_el(dtree_lg_el<type> *tel) {
  if (tel == NULL) return;
  delete_el(tel->right);
  delete_el(tel->left);
  delete tel;
}

template<class type>
void dtree_lg<type>::decompose(dtree_lg_el<type> *t, type *sarr, long nd, long &iter) {
  if (t == NULL) return;
  if (iter >= nd) return;
  decompose(t->left, sarr, nd, iter);
  sarr[iter]=t->value;
  iter++;
  decompose(t->right, sarr, nd, iter);
  
}


template<class type>
long dtree_lg<type>::add(type data) {
  dtree_lg_el<type> *current;

  if (trunk==NULL) {
    trunk=new dtree_lg_el<type>;
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
          current->left=new dtree_lg_el<type>;
          current=current->left;
          break;
        } else {
          current=current->left;
        }
      } else {
        if (current->right == NULL) {
          current->right=new dtree_lg_el<type>;
          current=current->right;
          break;
        } else {
          current=current->right;
        }
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
long dtree_lg<type>::nel() {
  return n;
}

template<class type>
void dtree_lg<type>::decompose(type *sarr, long nd) {
  long i=0;
  decompose(trunk, sarr, nd, i);
}


template<class type>
long dtree_lg<type>::delete_least() {
  dtree_lg_el<type> *current;

  if (trunk == NULL) return n;

  if (trunk->left == NULL) {
    if (trunk->right == NULL) {
      delete trunk;
      trunk=NULL;
      n=0;
    } else {
      current=trunk->right;
      delete trunk;
      trunk=current;
      trunk->parent=NULL;

      //then step down until we find the new least:
      while (current->left != NULL) {
        current=current->left;
      }
      least=current;
    }
  } else {
    //step up one:
    current=least->parent;
    current->left=least->right;
    if (current->left != NULL) current->left->parent=current;
    delete least;
    //then step down until we find the new least:
    while (current->left != NULL) {
      current=current->left;
    }
    least=current;

  }
  n--;
  return n;
}

template<class type>
long dtree_lg<type>::delete_greatest() {
  dtree_lg_el<type> *current;

  if (trunk == NULL) return n;

  if (trunk->right == NULL) {
    if (trunk->left == NULL) {
      delete trunk;
      trunk=NULL;
      n=0;
    } else {
      current=trunk->left;
      delete trunk;
      trunk=current;
      trunk->parent=NULL;
      while (current->right != NULL) {
        current=current->right;
      }
      greatest=current;
    }
  } else {
    current=greatest->parent;
    current->right=greatest->left;
    if (current->right != NULL) current->right->parent=current;
    delete greatest;
    while (current->right != NULL) {
      current=current->right;
    }
    greatest=current;

  }
  n--;
  return n;
}

template class dtree_lg<float>;
template class dtree_lg<char *>;

} //end namespace libpetey

