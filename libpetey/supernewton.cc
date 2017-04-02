#include <stdio.h>

#include <math.h>

#include "gsl/gsl_linalg.h"
#include "gsl/gsl_poly.h"

#include "error_codes.h"
#include "supernewton.h"

namespace libpetey {

//note: for smooth ("well-behaved") functions with local non-zero third
//moments, the algorithm has a tendency to land repeatedly on one side
//of the root, slowing down convergence.
//This issue has been addressed by throwing in a bisection step if the 
//approximated roots are getting repeatedly stuck on one side:

//maximum number of re-brackets on the same side of the root:
int glob_supnewt_maxnrebr=3;

template <class real>
real nearest_of_three(real value,
			real x1,
			real x2,
			real x3) {

  real diff1, diff2, diff3;
  real result;

  diff1=fabs(value-x1);
  diff2=fabs(value-x2);
  diff3=fabs(value-x3);

  if (diff1 > diff2) {
    if (diff2 > diff3) result = x3; else result = x2;
  } else {
    if (diff1 > diff3) result = x3; else result = x1;
  }

  return result;

}

//minimization function which brackets the root
//and then approximates root by fitting a third-
//order polynomial

template <class real>
real supernewton(void (*funcd) (real, void *, real *, real *),
		void *params,	//fucntion parameters
		real x1,	//first bracket
		real x2,	//second bracket
		real xtol,	//desired tolerance in x direction
		real ytol,	//desired tolerance in y direction
		long maxiter,	//maximum number of iterations
		supernewton_stat *err,	//error code and number of iterations
		real &y1,	//to avoid re-calculating
		real &dydx1,	//these items
		real y2,
		real dydx2)
{

  //polynomial coefficients:
  double a1, b1;	//might as well carry over from double to double

  real x0;		//solution

  real y0;		//solution value for y
  real dydx0;		//solution value for dydx

  //intermediates in calculation:
  real dx;
  long nroot=1;

  //roots of polynomial:
  double dx0_1, dx0_2, dx0_3;	//these need to be double because they're 
  				//returned from a GSL routine

  real x0_4;

  //x and y error:
  real xerr, yerr;
  real xerr1, xerr2;

  //approximate slope of the function:
  real m;

  int gsl_status;	//error state from GSL calls

  int nss=0;		//number of re-brackets on the same side

  //initialize status indicator:
  err->code=0;
  err->niter=0;
  err->nbis=0;

  do {
    if (y1*y2>0) {
      fprintf(stderr, "supernewton: ordinates must have opposite sign\n");
      err->code=PARAMETER_OUT_OF_RANGE;
      return (x1+x2)/2;
    }

/*
    printf("Brackets: x: [%g, %g] dx=%g\n", x1, x2, xerr);
    printf("          [%g, %g]\n", y1, y2);
    printf("          [%g, %g]\n", dydx1, dydx2);
*/

    if (abs(nss) >= glob_supnewt_maxnrebr) {
      //if the root falls one too many times on the same side, we just take a bisection step:
      x0=(x2+x1)/2;
      err->nbis++;		//count number of bisection steps
      nss=0;
    } else {
      //a simple coordinate shift vastly simplifies this
      //calculation (rather than using a matrix solver):
      dx=x2-x1;
      b1=(3*(y2-y1)/dx-2*dydx1-dydx2)/dx;
      a1=(dydx2-dydx1-2*b1*dx)/3/dx/dx;

      {
        //printf("Polynomial coeffs.:%f, %f, %f, %f\n", a1, b1, dydx1, y1);

        //solve the cubic:
        nroot=gsl_poly_solve_cubic(b1/a1, dydx1/a1, y1/a1, &dx0_1, &dx0_2, &dx0_3);

        if (nroot==1) {
          //if there is only one root AND it's finite AND it advances the solution
          //we use that one:
          x0=(x1+x2)/2;
          //printf("Root: %f; central value: %f\n", x1+dx0_1, x0);
          if (isfinite(dx0_1) && dx0_1>0 && x1+dx0_1 < x2) {
            x0=x1+dx0_1;
          } else {
            err->nbis++;
          }
        } else if (nroot==3) {
          real xdiff, xdiffmin;

          x0_4=(x1+x2)/2;
          //printf("Roots: %f, %f, %f; central value: %f\n", x1+dx0_1, x1+dx0_2, x1+dx0_3, x0_4);
          //x0=nearest_of_three(x0_4, (real) x0_1, (real) x0_2, (real) x0_3);

          //if there are three,

          //use the root that's closest to the centre and still inside the brackets:
          //(what about evaluating each of them and seeing which is closest to 0?)
          x0=x0_4;
          xdiffmin=fabs((x1-x2)/2);
          if (isfinite(dx0_1)) {
            xdiff=fabs(x1+dx0_1-x0_4);
            if (xdiff < xdiffmin) {
              xdiffmin=xdiff;
              x0=x1+dx0_1;
            }
          }
          if (isfinite(dx0_2)) {
            xdiff=fabs(x1+dx0_2-x0_4);
            if (xdiff < xdiffmin) {
              xdiffmin=xdiff;
              x0=x1+dx0_2;
            }
          }
          if (isfinite(dx0_3)) {
            xdiff=fabs(x1+dx0_3-x0_4);
            if (xdiff < xdiffmin) {
              xdiffmin=xdiff;
              x0=x1+dx0_3;
            }
          }
          if (x0==x0_4) err->nbis++;
        } else {
          x0=(x1+x2)/2;
          err->nbis++;
        }
      }

    }

    //evaluate the function at the approximated root
    (*funcd) (x0, params, &y0, &dydx0);
    //test for y-convergence:
    yerr=fabs(y0);
    if (yerr < ytol) break;

    //printf("Brackets: [%f, %f]\n", x1, x2);
    //rebracket the "true" root:
    if (y0*y1 > 0) {
      if (nss > 0) nss++; else nss=1;	//count how many times root falls on same side
      x1=x0;
      y1=y0;
      dydx1=dydx0;
    } else {
      if (nss < 0) nss--; else nss=-1;
      x2=x0;
      y2=y0;
      dydx2=dydx0;
    }

    //test for x-convergence:
    xerr=fabs(2*(x1-x2)/(x1+x2));
    if (xerr < xtol) break;

    err->niter++;		//count number of iterations
    //test for maximum iterations:
    if (err->niter > maxiter) {
      fprintf(stderr, "Maximum number of iterations exceeded in supernewton\n");
      err->code = MAX_ITER_EXCEEDED;
      break;
      printf("Brackets: [%f, %f]\n", x1, x2);
      printf("Polynomial coeffs.:%f, %f, %f, %f\n", a1, b1, dydx1, y1);
      printf("Roots: %f, %f, %f\n", dx0_1, dx0_2, dx0_3);
    }
  } while(1);

  //place the function values in the fourth and third last parameters:
  y1=y0;
  dydx1=dydx0;

  //printf("Supernewton: %d iterations required to reach convergence\n", i);
  //if (err->nbis>0) fprintf(stderr, "Supernewton: %d bisection steps\n", nbis);

  return x0;

}

template float supernewton<float>(void (*funcd) (float, void *, float *, float *),
		void *params,
		float x1,
		float x2,
		float xtol,
		float ytol,
		long maxiter,
		supernewton_stat *err,
		float &y1,
		float &dydx1,
		float y2,
		float dydx2);

template double supernewton<double>(void (*funcd) (double, void *, double *, double *),
		void *params,
		double x1,
		double x2,
		double xtol,
		double ytol,
		long maxiter,
		supernewton_stat *err,
		double &y1,
		double &dydx1,
		double y2,
		double dydx2);

} //end namespace libpetey

