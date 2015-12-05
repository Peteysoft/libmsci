#include <math.h>
#include <stdio.h>

#ifndef CLUSTER_LIB_H
#define CLUSTER_LIB_H 1

namespace libagf {

  //clustering analysis using knn and a threshold density:
  template <class real, class cls_t>
  void cluster_knn(real **x, 		//sample vectors
		dim_ta d,		//dimensionality
		nel_ta n, 		//number of samples
		nel_ta k, 		//number of nearest neighbours
		real pth, 		//threshold density
		cls_t *cnum);		//returned cluster number

  template <class leaf_t, class metric_t>
  class cluster_node;

  template <class leaf_t, class metric_t>
  union cluster_branch {
    cluster_node<leaf_t, metric_t> *branch[2];
    leaf_t leaf;
  };

  template <class leaf_t, class metric_t>
  class cluster_node {
    protected:
      cluster_branch<leaf_t, metric_t> child;
      cluster_node<leaf_t, metric_t> *parent;
      float d;
    public:
      //initialize with a pair of points:
      cluster_node(metric_t d1,
		    cluster_node<leaf_t, metric_t> *child1, 
		    cluster_node<leaf_t, metric_t> *child2);
      //"leaf"= location of a sample:
      cluster_node(leaf_t leaf);

      //check integrity of dendrogram:
      cluster_node<leaf_t, metric_t> *check();

      void print(FILE *fs, int depth, char *opt=NULL);

      //recursive delete:
      void del();

      //pointers to parent, children and top-most node:
      cluster_node<leaf_t, metric_t> *up();
      cluster_node<leaf_t, metric_t> *left();
      cluster_node<leaf_t, metric_t> *right();
      cluster_node<leaf_t, metric_t> *top();

      metric_t val(leaf_t &leaf);
      int count();		//number of "leaves"
      //locations of leaves below this node:
      int gather(leaf_t *vec, int maxn);
      //split into classes (?) -- not finished
      int split(cluster_node<leaf_t, metric_t> **root, int ncls);
  };

  template <class real, class cls_t>
  class cluster_tree {
    protected:
      int n;		//number of vectors
      cluster_node<int, real> **leaf;		//top "leaf" nodes
      cluster_node<int, real> *root;
    public:
      cluster_tree();
      ~cluster_tree();

      //build the tree using all the distances:
      int build_all(real *d, int nvec, int nvar);
      int build_all(real **vec, int nvec, int nvar);

      //gets class labels once they've been marked:
      int get_classes(cls_t *cls, cls_t ncls);

      //interactive browsing utility:
      int browse(cls_t *cls, FILE *fs=NULL);

      //print out the tree:
      void print(FILE *fs, char *opt=NULL);
  };
  
}

#endif

