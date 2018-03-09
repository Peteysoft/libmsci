#ifndef MULTICLASS_METHODS_H
#define MULTICLASS_METHODS_H
//methods for solving for multi-class

namespace libagf {
  template <typename code_t, typename real>
  void prep_nonstrict(code_t **a, 	//coding matrix
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
  template <typename code_t, typename real>
  void solve_class_scratch(code_t **a, 	//coding matrix
		  int m, 		//number of binary classifiers
		  int n, 		//number of classes
		  real *r, 		//binary class results
		  real *p);		//returned probabilities

  //vote from labels:
  template <typename code_t, typename real>
  void solve_class_vote(code_t **a, int m, int n, real *r, real *p);

  //vote from probabilities:
  template <typename code_t, typename real>
  void solve_class_vote_pdf(code_t **a, int m, int n, real *r, real *p);

  //least square applying only normalization constraint:
  template <typename code_t, typename real>
  void solve_class_norm1(code_t **a, int m, int n, real *r, real *p);

  //appropriate for 1 vs. 1:
  template <typename code_t, typename real>
  void solve_class_norm2(code_t **a, int m, int n, real *r, real *p);

  //two attempts at fully constrained, general solutions:
  template <typename code_t, typename real>
  void solve_class_constrained1(code_t **a, int m, int n, real *r, real *p);

  template <typename code_t, typename real>
  void solve_class_constrained2(code_t **a, int m, int n, real *r, real *p);

  //method I developed for orthogonal coding matrices (applied as general
  //method):
  //template <typename code_t, typename real>
  //void solve_class_renorm(code_t **a, int m, int n, real *r, real *p);
  template <typename code_t>
  void solve_class_renorm(code_t **a, int m, int n, float *r, float *p);
  template <typename code_t>
  void solve_class_renorm(code_t **a, int m, int n, double *r, double *p);

  //method I developed for orthogonal coding matrices (mainly appropriate for
  //same):
  template <typename code_t, typename real>
  void solve_class_vote_pdf2(code_t **a, int m, int n, real *r, real *p);

  //one-versus-rest: basic solution; if you ignore "the rest" then the 
  //binary probabilities are the same as the multi-class probs.
  template <typename code_t, typename real>
  void solve_class_1vR(code_t **a, int m, int n, real *r, real *p);

  //(Zadrozny 2002)
  template <typename code_t, typename real>
  void solve_class_Zadrozny(code_t **a, int m, int n, real *r, real *p);

  //convert to logarithms, solve by Newton's method:
  template <typename code_t, typename real>
  void solve_class_interior(code_t **a, int m, int n, real *r, real *p);

}

#endif

