#include <math.h>
  
#include "gsl/gsl_sf.h"

#include "error_codes.h"
#include "kextreme.h"
#include "quicksort.h"
#include "agf_metric2.h"

#include "peteys_tmpl_lib.h"

#include "cluster_lib.h"

using namespace std;
using namespace libpetey;

namespace libagf {

template <class leaf_t, class metric_t>
cluster_node<leaf_t, metric_t>::cluster_node(metric_t d1, 
			cluster_node *child1, 
			cluster_node *child2) {
  cluster_node<leaf_t, metric_t> *swp;
  d=d1;
  child.branch[0]=child1;
  //search to the root of the tree:
  while (child.branch[0]->parent!=NULL) {
    //printf("left: up\n");
    child.branch[0]=child.branch[0]->parent;
  }
  child.branch[1]=child2;
  while (child.branch[1]->parent!=NULL) {
    //printf("right: up\n");
    child.branch[1]=child.branch[1]->parent;
  }

  //left branch always has more leaves:
  if (child.branch[0]!=child.branch[1]) {
    child.branch[0]->parent=this;
    child.branch[1]->parent=this;
    if (child.branch[0]->count() < child.branch[1]->count()) {
      swp=child.branch[0];
      child.branch[0]=child.branch[1];
      child.branch[1]=swp;
    }
  }

  parent=NULL;

}

//check if it's a redundant node:
//(not that most of this code isn't pretty fucking redundant...)
template <class leaf_t, class metric_t>
cluster_node<leaf_t, metric_t> *cluster_node<leaf_t, metric_t>::check() {
  if (child.branch[0]!=child.branch[1]) {
    return this;
  } else {
    return child.branch[0];
  }
}

template <class leaf_t, class metric_t>
cluster_node<leaf_t, metric_t>::cluster_node(leaf_t leaf) {
  child.leaf=leaf;
  d=-1;
  parent=NULL;
}

template <class leaf_t, class metric_t>
void cluster_node<leaf_t, metric_t>::print(FILE *fs, int depth, char *opt) {
  for (int i=0; i<depth; i++) printf("  ");
  if (d<0) {
    fprintf(fs, "%d\n", child.leaf);
  } else {
    if (opt==NULL) {
      fprintf(fs, "%g\n", d);
    } else {
      //generate multi-borders control file:
      fprintf(fs, "%s {\n", opt);
    }
    child.branch[0]->print(fs, depth+1, opt);
    child.branch[1]->print(fs, depth+1, opt);
    for (int i=0; i<depth; i++) printf("  ");
    fprintf(fs, "}\n", opt);
  }
}

template <class leaf_t, class metric_t>
cluster_node<leaf_t, metric_t> *cluster_node<leaf_t, metric_t>::up() {
  return parent;
}

template <class leaf_t, class metric_t>
cluster_node<leaf_t, metric_t> *cluster_node<leaf_t, metric_t>::left() {
  if (d < 0) return NULL;
  return child.branch[0];
}

template <class leaf_t, class metric_t>
cluster_node<leaf_t, metric_t> *cluster_node<leaf_t, metric_t>::right() {
  if (d < 0) return NULL;
  return child.branch[1];
}

template <class leaf_t, class metric_t>
cluster_node<leaf_t, metric_t> *cluster_node<leaf_t, metric_t>::top() {
  if (parent==NULL) return this; else return parent->top();
}

template <class leaf_t, class metric_t>
metric_t cluster_node<leaf_t, metric_t>::val(leaf_t &leaf) {
  if (d<0) leaf=child.leaf;
  return d;
}

template <class leaf_t, class metric_t>
void cluster_node<leaf_t, metric_t>::del() {
  if (d>=0) {
    child.branch[0]->del();
    child.branch[1]->del();
  }
}

template <class leaf_t, class metric_t>
int cluster_node<leaf_t, metric_t>::count() {
  //printf("d=%g\n", d);
  if (d<0) {
    //printf("item=%d\n", child.leaf);
    return 1;
  }
  return child.branch[0]->count()+child.branch[1]->count();

}

template <class leaf_t, class metric_t>
int cluster_node<leaf_t, metric_t>::gather(leaf_t *vec, int maxn) {
  int n;

  if (d<0) {
    vec[0]=child.leaf;
    n=1;
  } else {
    n=child.branch[0]->gather(vec, maxn);
    if (n>=maxn) return n;
    n+=child.branch[1]->gather(vec+n, maxn-n);
  }

  return n;
}

template <class leaf_t, class metric_t>
int cluster_node<leaf_t, metric_t>::split(cluster_node<leaf_t, metric_t> **root, int ncls) {
  int nlb2;		//n log base 2...
  int ncls2;		//number to search for one level up
  int nfound;		//number actually found

  //either we're at the max. depth or we've found a single-member class:
  if (d<0 || ncls==1) {
    root[0]=this;
    return 1;
  } 

  nlb2=log(ncls)/log(2);
  ncls2=pow(2, nlb2-1);

  nfound=child.branch[0]->split(root, ncls2);
  nfound+=child.branch[1]->split(root+nfound, ncls2);

//**** not finished ****
  if (ncls-2*ncls2 > 0) {
    metric_t d[2*ncls2];
    long *sind;
    for (int i=0; i<2*ncls2; i++) d[i]=root[i].d;
    sind=heapsort(d, 2*ncls2);
    for (int i=2*ncls2-1; i>=0; i++) {}
  }
//**** not finished ****

  return nfound;

}

void trimat_coord(int n, int &j, int &k) {
  j=(-1.+sqrt(1.+8.*n))/2.;
  k=n - (j+1)*j/2;
  j++;
}

template <class real, class cls_t>
cluster_tree<real, cls_t>::cluster_tree() {
  n=0;
  leaf=NULL;
  root=NULL;
}

template <class real, class cls_t>
cluster_tree<real, cls_t>::~cluster_tree() {
  root->del();
  delete [] leaf;
}

template <class real, class cls_t>
int cluster_tree<real, cls_t>::build_all(real *d, int nvec, int nvar) {
  cluster_node<int, real> *check;
  int nd;
  long *sind;
  int j, k;
 
  n=nvec;

  nd=n*(n-1)/2; 

  sind=heapsort(d, nd);

  leaf=new cluster_node<int, real> *[n];
  //tree stores only indices at the leaf-level, not entire vectors:
  for (int i=0; i<n; i++) leaf[i]=new cluster_node<int, real>(i);

  printf("Building tree:         %6.1f%%", 0.);
  for (int i=0; i<nd; i++) {
    trimat_coord(sind[i], j, k);
    //printf("%d %d; %g\n", j, k, d[sind[i]]);
    check=new cluster_node<int, real>(d[sind[i]], leaf[j], leaf[k]);
    root=check->check();
    if (root!=check) delete check;
    printf("\b\b\b\b\b\b\b%6.1f%%", 100.*(i+1)/nd);
    fflush(stdout);
  }
  printf("\n");

  return 0;
}

template <class real, class cls_t>
int cluster_tree<real, cls_t>::build_all(real **vec, int nvec, int nvar) {
  cluster_node<int, real> *check;
  real *d;
  int nd;
  long *sind;
  int j, k;
 
  n=nvec;

  nd=n*(n-1)/2; 
  d=new real[nd];

  printf("Calculating distances: %6.1f%%", 0.);

  for (int i=0; i<nd; i++) {
    trimat_coord(i, j, k);
    d[i]=sqrt(metric2(vec[j], vec[k], nvar));
    //printf("%d %d; %g\n", j, k, d[i]);
    printf("\b\b\b\b\b\b\b%6.1f%%", 100.*(i+1)/nd);
    fflush(stdout);
  }
  printf("\n");

  build_all(d, n, nvar);

  return 0;
}

template <class real, class cls_t>
int cluster_tree<real, cls_t>::get_classes(cls_t *cls, cls_t ncls) {
  cluster_node<int, real> *subroot[ncls];		//roots for each class
  int nal2;
  cluster_node<int, real> **sri;

  nal2=log(ncls)/log(2);
  sri=new cluster_node<int, real> *[nal2];

}

//#include <ncurses/curses.h>

template <class real, class cls_t>
int cluster_tree<real, cls_t>::browse(cls_t *cls, FILE *fs) {
  cluster_node<int, real> *cur;		//current node
  int level=0;				//depth
  char c=0;				//character read from keyboard
  int ind0;				//index of current leaf, if applicable
  real d;				//distance between clusters
  int ind[n];				//indices of items in current cluster
  int nfound;				//number of items in current cluster
  cls_t nwcl;				//label for items in current cluster

  for (int i=0; i<n; i++) cls[i]=0;

  cur=root;

  //initscr();
  //cbreak();
  //timeout(-1);

  do {
    if (c!='\n') {
      printf("\n");
      printf("At level %d; %d items total\n", level, n);
      d=cur->val(ind0);
      if (d < 0) {
        printf("Bottom-most node: item number %d has been labelled %d\n", ind0, cls[ind0]);
      } else {
        printf("%d items in left sub-tree; %d items in right sub-tree\n",
		    (cur->left())->count(), (cur->right())->count());
        printf("Clusters are %g apart\n", d);
      }
      printf("Choose one of: 0-9, c, h, l, r, u; h is for help\n");
    }
    //refresh();
    c=getc(stdin);
    //c=getch();
    switch (c) {
      case('l'):
        if (cur->left()==NULL) {
          printf("Already at bottom, cannot go down\n");
        } else {
          cur=cur->left();
	  level++;
	}
	break;
      case('r'):
        if (cur->right()==NULL) {
          printf("Already at bottom, cannot go down\n");
        } else {
          cur=cur->right();
	  level++;
	}
	break;
      case('u'):
        if (cur->up()==NULL) {
          printf("Already at top, cannot go up\n");
        } else {
          cur=cur->up();
	  level--;
	}
	break;
      case('0'):
	nfound=cur->gather(ind, n);
	printf("Labelling %d items as %d\n", nfound, 0);
	for (int i=0; i<nfound; i++) cls[ind[i]]=0;
	break;
      case('1'):
	nfound=cur->gather(ind, n);
	printf("Labelling %d items as %d\n", nfound, 1);
	for (int i=0; i<nfound; i++) cls[ind[i]]=1;
	break;
      case('n'):
        break;
      case('2'):
	nfound=cur->gather(ind, n);
	printf("Labelling %d items as %d\n", nfound, 2);
	for (int i=0; i<nfound; i++) cls[ind[i]]=2;
	break;
      case('3'):
	nfound=cur->gather(ind, n);
	printf("Labelling %d items as %d\n", nfound, 3);
	for (int i=0; i<nfound; i++) cls[ind[i]]=3;
	break;
      case('4'):
	nfound=cur->gather(ind, n);
	printf("Labelling %d items as %d\n", nfound, 4);
	for (int i=0; i<nfound; i++) cls[ind[i]]=4;
	break;
      case('5'):
	nfound=cur->gather(ind, n);
	printf("Labelling %d items as %d\n", nfound, 5);
	for (int i=0; i<nfound; i++) cls[ind[i]]=5;
	break;
      case('6'):
	nfound=cur->gather(ind, n);
	printf("Labelling %d items as %d\n", nfound, 6);
	for (int i=0; i<nfound; i++) cls[ind[i]]=6;
	break;
      case('7'):
	nfound=cur->gather(ind, n);
	printf("Labelling %d items as %d\n", nfound, 7);
	for (int i=0; i<nfound; i++) cls[ind[i]]=7;
	break;
      case('8'):
	nfound=cur->gather(ind, n);
	printf("Labelling %d items as %d\n", nfound, 8);
	for (int i=0; i<nfound; i++) cls[ind[i]]=8;
	break;
      case('9'):
	nfound=cur->gather(ind, n);
	printf("Labelling %d items as %d\n", nfound, 9);
	for (int i=0; i<nfound; i++) cls[ind[i]]=9;
	break;
      case ('c'):
	//nocbreak();
	printf("Enter class number:");
	scanf("%d", &nwcl);
	printf("\n");
	nfound=cur->gather(ind, n);
	printf("Labelling %d items as %d\n", nfound, nwcl);
	for (int i=0; i<nfound; i++) cls[ind[i]]=nwcl;
	//cbreak();
	break;
      case ('q'):
	break;
      case ('h'):
	printf("0-9, c = choose class label for current cluster\n");
	printf("h      = print this help screen\n");
	printf("l      = move down to left sub-tree\n");
	printf("r      = move down to right sub-tree\n");
	printf("u      = move up to parent node\n");
        printf("n      = mark node (in conjunction with log file)\n");
	continue;
      case ('\n'):
	continue;
      default:
	printf("Unrecognized option: type h for help.\n");
        continue;
    }
    if (fs!=NULL) fprintf(fs, "%c\n", c);
  } while (c!='q');

  //endwin();

}

template <class real, class cls_t>
void cluster_tree<real, cls_t>::print(FILE *fs, char *opt) {
  root->print(fs, 0, opt);
}


template class cluster_tree<float, int32_t>;
template class cluster_tree<double, int32_t>;

//recursively fill the clusters:
template <class real, class cls_t>
void fill_cluster(real **x, 
		dim_ta d,
		nel_ta n,
		nel_ta k,
		real pth,
		cls_t &cnum,
		nel_ta sub,
		cls_t *cnum_all,
		flag_a *calc, 
		flag_a inflag)
{
  real *d2;
  real *kl;
  long ind[k+1];
  real pdf;
  real r;
  real sqrtpir;
  real V;

  if (calc[sub]==1) return;

  //calculate distances:
  d2=new real[n];
  for (nel_ta j=0; j<n; j++) d2[j]=metric2(x[sub], x[j], d);
  //find nearest-neighbours:
  kl=new real[k+1];
  KLEAST_FUNC(d2, n, k+1, kl, ind);
  delete [] d2;			//run out of heap pretty quick otherwise

  //calculate pdf:
  r=(sqrt(kl[k-1])+sqrt(kl[k]))/2;
  sqrtpir=sqrt(M_PI)*r;
  V=1;
  for (dim_ta i=0; i<d; i++) V*=sqrtpir;
  V=V/gsl_sf_gamma(d/2.+1);
  pdf=(real) k/V/(real) n;

  calc[sub]=1;

  delete [] kl;

  if (pdf >= pth) {
    if (inflag) cnum++;
    cnum_all[sub]=cnum;
    for (nel_ta j=0; j<=k; j++) 
	  fill_cluster(x, d, n, k, pth, cnum, ind[j], cnum_all, calc, 0);
  } else {
    cnum_all[sub]=0;
  }
}

//clustering analysis using knn and a threshold density:
template <class real, class cls_t>
void cluster_knn(real **x, 		//sample vectors
		dim_ta d,		//dimensionality
		nel_ta n, 		//number of samples
		nel_ta k, 		//number of nearest neighbours
		real pth, 		//threshold density
		cls_t *cnum)		//cluster number
{
  flag_a calc[n];		//sample has been evaluated?
  cls_t cluster=0;		//cluster number

  for (nel_ta j=0; j<n; j++) calc[j]=0;

  for (nel_ta i=0; i<n; i++) {
    //if (calc[i] == 1) continue;
    //cluster++;
    fill_cluster(x, d, n, k, pth, cluster, i, cnum, calc, 1);
  }

}

template void cluster_knn<float, cls_ta>(float **x, dim_ta d, nel_ta n,
		nel_ta k, float pth, cls_ta *cnum);

template void cluster_knn<double, cls_ta>(double **x, dim_ta d, nel_ta n,
		nel_ta k, double pth, cls_ta *cnum);

}

