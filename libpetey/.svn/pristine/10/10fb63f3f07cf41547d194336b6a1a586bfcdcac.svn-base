#include <stdio.h>

#include <math.h>

#include "error_codes.h"
#include "roots_mins.h"

namespace libpetey {

template <class real>
real min_golden(void (*func) (real, void *, real *),
		void *params,	//fucntion parameters
		real ax,	//first bracket
		real bx,	//second bracket
		real cx,
		real tol,	//desired tolerance 
		long maxiter,	//maximum number of iterations
		long &k,		//error code and number of iterations
		real &ymin)	//to avoid re-calculating

{
  real C = (3-sqrt(5))/2;
  real R = 1-C;

  real x0, x1, x2, x3;
  real f1, f2;
  real xmin;
 
  x0 = ax;
  x3 = cx;
  if (fabs(cx-bx) < fabs(bx-ax)) {
    x1 = bx;
    x2 = bx + C*(cx-bx);
  } else {
    x2 = bx;
    x1 = bx - C*(bx-ax);
  }
  (*func)(x1, params, &f1);
  (*func)(x2, params, &f2);
 
  k = 1;
  while (fabs(x3-x0) > tol*(fabs(x1)+fabs(x2))) {
    printf("%g\n", x0);
    printf("%g %g\n", x1, f1);
    printf("%g %g\n", x2, f2);
    printf("%g\n", x3);
    printf("\n");
    if (f2 < f1) {
      x0 = x1;
      x1 = x2;
      x2 = R*x1 + C*x3;   // x2 = x1+c*(x3-x1)
      f1 = f2;
     (*func)(x2, params, &f2);
    } else {
      x3 = x2;
      x2 = x1;
      x1 = R*x2 + C*x0;   // x1 = x2+c*(x0-x2)
      f2 = f1;
      (*func)(x1, params, &f1);
    }
    k = k+1;
    if (k>maxiter) {
      fprintf(stderr, "min_golden: maximum number of iterations (%d) exceeded\n", maxiter);
      break;
    }
  }
 
  if (f1 < f2) {
    xmin = x1;
    ymin = f1;
  } else {
    xmin = x2;
    ymin = f2;
  }

  return xmin;
}


//minimization function which brackets the root
//and interpolates between the brackets:
template <class real>
real root_false_position(void (*func) (real, void *, real *),
		void *params,	//fucntion parameters
		real x1,	//first bracket
		real x2,	//second bracket
		real xtol,	//desired tolerance in x direction
		real ytol,	//desired tolerance in y direction
		long maxiter,	//maximum number of iterations
		long &err,	//error code and number of iterations
		real &y1,	//to avoid re-calculating
		real y2)
{

  real x0;		//solution
  real y0;		//solution value for y

  //x and y error:
  real xerr, yerr;

  if (y1*y2>0) {
    fprintf(stderr, "supernewton: ordinates must have opposite sign\n");
    err=-1;
    return (x1+x2)/2;
  }

  err=0;

  do {
/*
    printf("Brackets: x: [%g, %g] dx=%g\n", x1, x2, xerr);
    printf("          [%g, %g]\n", y1, y2);
    printf("          [%g, %g]\n", dydx1, dydx2);
*/

    x0=x1+y1*(x2-x1)/(y2-y1);
    (*func) (x0, params, &y0);

    //test for convergence:
    yerr=fabs(y0);
    if (yerr < ytol) break;

    //printf("Brackets: [%f, %f]\n", x1, x2);
    //rebracket the "true" root:
    if (y0*y1 > 0) {
      x1=x0;
      y1=y0;
    } else {
      x2=x0;
      y2=y0;
    }

    //test for convergence:
    xerr=fabs(2*(x1-x2)/(x1+x2));
    if (xerr< xtol) break;

    err++;
    if (err > maxiter) {
      fprintf(stderr, "root_false_position: Maximum number of iterations exceeded (%d)\n", maxiter);
      break;
      printf("Brackets: [%f, %f]\n", x1, x2);
    }
  } while(1);

  //place the function value in the second last parameter:
  y1=y0;

  //printf("root_false_position: %d iterations required to reach convergence\n", i);

  return x0;

}

template float min_golden<float>(void (*func) (float, void *, float *),
		void *params,
		float x1,
		float x2,
		float x3,
		float xtol,
		long maxiter,
		long &err,
		float &y2);

template double min_golden<double>(void (*func) (double, void *, double *),
		void *params,
		double x1,
		double x2,
		double x3,
		double xtol,
		long maxiter,
		long &err,
		double &y1);

template float root_false_position<float>(void (*func) (float, void *, float *),
		void *params,
		float x1,
		float x2,
		float xtol,
		float ytol,
		long maxiter,
		long &err,
		float &y1,
		float y2);

template double root_false_position<double>(void (*func) (double, void *, double *),
		void *params,
		double x1,
		double x2,
		double xtol,
		double ytol,
		long maxiter,
		long &err,
		double &y1,
		double y2);

} //end namespace libpetey

