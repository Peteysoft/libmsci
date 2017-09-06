#ifndef MULTICLASS_METHODS_H
#define MULTICLASS_METHODS_H
//methods for solving for multi-class

namespace libagf {
  template <class real>
  void prep_nonstrict(real **a, 	//coding matrix
		  int m, 		//number of binary classifiers
		  int n, 		//number of classes
		  real *r,		//binary classification results
		  gsl_matrix *q,	//linear problem
		  gsl_vector *b);

  //normalize class probabilities:

  //standard normalization: zero negative probalities; divide by total
  template <class real>
  void p_renorm(real *p, 
		  int n);

  //motion is always equal in all probabilities; negative probabilities
  //are zeroed and removed from the calculation
  template <class real>
  void p_constrain_renorm1(real *p, int n);

  template <class real>
  void p_constrain_renorm1a(real *p, int n);

  //basic least-squares solution:
  template <class real>
  void solve_class_scratch(real **a, 	//coding matrix
		  int m, 		//number of binary classifiers
		  int n, 		//number of classes
		  real *r, 		//binary class results
		  real *p);		//returned probabilities

  //vote from labels:
  template <class real>
  void solve_class_vote(real **a, int m, int n, real *r, real *p);

  //vote from probabilities:
  template <class real>
  void solve_class_vote_pdf(real **a, int m, int n, real *r, real *p);

  //least square applying only normalization constraint:
  template <class real>
  void solve_class_norm1(real **a, int m, int n, real *r, real *p);

  //appropriate for 1 vs. 1:
  template <class real>
  void solve_class_norm2(real **a, int m, int n, real *r, real *p);

  //two attempts at fully constrained, general solutions:
  template <class real>
  void solve_class_constrained1(real **a, int m, int n, real *r, real *p);

  template <class real>
  void solve_class_constrained2(real **a, int m, int n, real *r, real *p);

  //method I developed for orthogonal coding matrices (applied as general
  //method):
  template <class real>
  void solve_class_renorm(real **a, int m, int n, real *r, real *p);

  //method I developed for orthogonal coding matrices (mainly appropriate for
  //same):
  template <class real>
  void solve_class_vote_pdf2(real **a, int m, int n, real *r, real *p);

  //one-versus-rest: basic solution; if you ignore "the rest" then the 
  //binary probabilities are the same as the multi-class probs.
  template <class real>
  void solve_class_1vR(real **a, int m, int n, real *r, real *p);

}

#endif

