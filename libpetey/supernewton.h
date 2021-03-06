#ifndef SUPERNEWTON_H_INCLUDED
#define SUPERNEWTON_H_INCLUDED 1

namespace libpetey {

  //complicates things somewhat, but we need to know both
  //the number of iterations and the number of bisection steps
  //and to separate these from the error status:
  struct supernewton_stat {
    int niter;			//total number of iterations
    int nbis;			//number of bisection steps
    int code;			//error status
  };

  //root-finding function which brackets the root
  //and then approximates root by fitting a third-
  //order polynomial
  template <class real>
  real supernewton(void (*funcd) (real, void *, real *, real *),
		void *params,		//function parameters
                real x1,		//first bracket
		real x2,		//second bracket
		real xtol,		//desired (absolute) tolerance
		real ytol,
		long maxiter,		//maximum number of iterations
		supernewton_stat *err,	//returns error status
					//returns 0 or less for failure
		real &y1,		//to avoid redundant calculation
		real &dydx1,		//returns final values in these two
		real y2,
		real dydx2);

  template <class real>
  inline real supernewton(void (*funcd) (real, void *, real *, real *),
		void *params,		//function parameters
                real x1,		//first bracket
		real x2,		//second bracket
		real xtol,		//desired (absolute) tolerance
		real ytol,
		long maxiter,		//maximum number of iterations
		supernewton_stat *err)	//number of iterations/error code
  {

    real y1, y2, dydx1, dydx2;

    (*funcd) (x1, params, &y1, &dydx1);
    (*funcd) (x2, params, &y2, &dydx2);

    return supernewton(funcd, params, x1, x2, xtol, ytol, maxiter, err,
		  y1, dydx1, y2, dydx2);
  }

}

#endif

