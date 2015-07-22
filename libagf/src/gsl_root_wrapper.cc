#include <stdio.h>
#include <math.h>

#include "gsl/gsl_roots.h"

double gsl_root_wrapper(double x1, 
			double x2, 
			double tol, 
			long MAXIT, 
			int err, 
			double (* funcd) (double, void *), 
			void *params,
			gsl_root_fsolver_type *T=gsl_root_fsolver_brent);
  gsl_root_fsolver *sol;
  gsl_function F;
  int status;
  long iter;
  double r;

  //initialize the GSL solver:
  sol=gsl_root_fsolver_alloc(T);

  //try using one of the GSL bracketing routines:
  F.function=funcd;
  F.params=params;				//set function parameters
  gsl_root_fsolver_set(sol, &F, x1, x2);	//initialize solver
  iter=0;
  err=0;
  do {
    iter++;
    status = gsl_root_fsolver_iterate(sol);
    x1=gsl_root_fsolver_x_lower(sol);
    x2=gsl_root_fsolver_x_upper(sol);
    r = gsl_root_fsolver_root(sol);
    //status = gsl_root_test_residual(fabs(bfind2(t0, &t0)), tol);
    if (iter > MAXIT) {
      fprintf(stderr, "GSL solver failed to converge\n");
      err=1;
      return r;
    }
  } while (fabs(x2-x1) > tol && fabs(r) > tol);
  printf("Converged in %ld iterations\n", iter);

  return r;

}


