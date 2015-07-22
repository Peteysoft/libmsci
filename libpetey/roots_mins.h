#ifndef ROOTS_MINS_H_INCLUDED
#define ROOTS_MINS_H_INCLUDED 1

namespace libpetey {

  //minimization routine that uses bisection:
  template <class real>
  real min_golden(void (*func) (real, void *, real *),
		void *params,	//function parameters
                real x1,	//first bracket
		real x2,	//second bracket
		real x3,
		real tol,	//desired tolerance
		long maxiter,	//maximum number of iterations
		long &err,	//returns number of iterations
		real &y0);	//to avoid redundant calculation

  //root-finding function which brackets the root
  //and interpolates between the two brackets:
  template <class real>
  real root_false_position(void (*func) (real, void *, real *),
		void *params,	//function parameters
                real x1,	//first bracket
		real x2,	//second bracket
		real xtol,	//desired (absolute) tolerance
		real ytol,
		long maxiter,	//maximum number of iterations
		long &err,	//returns number of iterations
		real &y1,	//to avoid redundant calculation
		real y2);

  template <class real>
  inline real root_false_position(void (*func) (real, void *, real *),
		void *params,	//function parameters
                real x1,	//first bracket
		real x2,	//second bracket
		real xtol,	//desired (absolute) tolerance
		real ytol,
		long maxiter,	//maximum number of iterations
		long &err)	//returns number of iterations
  {

    real y1, y2;

    (*func) (x1, params, &y1);
    (*func) (x2, params, &y2);

    return root_false_position(func, params, x1, x2, xtol, ytol, maxiter, err,
		  y1, y2);
  }

}

#endif

