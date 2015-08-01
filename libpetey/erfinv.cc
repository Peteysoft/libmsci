#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_errno.h>
//#include <gsl/gsl_sf_erf.h>

#include "peteys_tmpl_lib.h"
#include "error_codes.h"
#include "supernewton.h"
#include "erfinv.h"

namespace libpetey {

  const double ERFINV_XTOL=1e-14;
  const double ERFINV_YTOL=0.;
  //works best if we just skip the complementary error function and do a straight inverse:
  const double ERFINV_XTRAN=gsl_sf_erf(ERFINV_YTRAN);
  const double ERFINV_EPS=1e-9;
  const double ERFINV_MAXITER=200;
  const double ERFINV_MAXY=6;           //value (approximately) at which erf returns exactly 1

  const double ERFCINV_XTRAN=gsl_sf_erfc(ERFCINV_YTRAN);
  const double ERFCINV_MAXY=27.5;         //value at which erfc returns exactly 0

  int erfinv_niter;		//evil global variable: number of iterations
  				//taken by root-finder
  int erfinv_nbis;		//number of bisection steps taken

  void erf_plus_deriv(double x, void *param, double *y, double *dydx) {
    *y=gsl_sf_erf(x)-*(double *) param; 
    *dydx=M_2_SQRTPI*exp(-x*x);
  }

  void erfc_plus_deriv(double x, void *param, double *y, double *dydx) {
    *y=gsl_sf_erfc(x)-*(double *) param; 
    *dydx=-M_2_SQRTPI*exp(-x*x);
  }

  double erfinv(double x) {
    double absx;			//absolute value of x (function is anti-symmetric)
    double sgnx;			//sign of x
    double y;			//returned inverse
    double y1, y2;		//root brackets
    supernewton_stat err;		//supernewton error/iterations

    gsl_error_handler_t *old_handler;
    old_handler=gsl_set_error_handler_off();

    absx=fabs(x);
    sgnx=x/absx;
    if (absx>=1) return sgnx*INFINITY;

    //not as fast as the look-up table approach, but simpler and still effective:
    if (absx>ERFINV_XTRAN) {
      x=1-absx;
      void *param=&x;
      y1=ERFINV_YTRAN;
      y2=ERFINV_MAXY;
      y=supernewton(&erfc_plus_deriv, param, y1, y2, ERFINV_XTOL, ERFINV_YTOL, ERFINV_MAXITER, &err);
    } else {
      y1=-ERFINV_EPS;
      y2=ERFINV_YTRAN+ERFINV_EPS;
      void *param=&absx;

      y=supernewton(&erf_plus_deriv, param, y1, y2, ERFINV_XTOL, ERFINV_YTOL, ERFINV_MAXITER, &err);
    }
    gsl_set_error_handler(old_handler);

    if (err.code!=0) {
      fprintf(stderr, "supernewton root-finder returned %d error code\n", err.code);
      exit(err.code);
    }

    erfinv_niter=err.niter;
    erfinv_nbis=err.nbis;

    return sgnx*y;

  }

  double erfcinv(double x) {
    double y;			//returned inverse
    double y1, y2;		//root brackets
    supernewton_stat err;		//supernewton error/iterations

    gsl_error_handler_t *old_handler;
    old_handler=gsl_set_error_handler_off();

    if (x>=2) return -INFINITY;
    if (x<=0) return INFINITY;

    //not as fast as the look-up table approach, but simpler and still effective:
    if (x>1) {
      return erfinv(1-x);
    } else if (x<ERFCINV_XTRAN) {
      void *param=&x;
      y1=ERFCINV_YTRAN;
      y2=ERFCINV_MAXY;
      y=supernewton(&erfc_plus_deriv, param, y1, y2, ERFINV_XTOL, ERFINV_YTOL, ERFINV_MAXITER, &err);
    } else {
      x=1-x;
      void *param=&x;
      y1=-ERFINV_EPS;
      y2=ERFCINV_YTRAN+ERFINV_EPS;

      y=supernewton(&erf_plus_deriv, param, y1, y2, ERFINV_XTOL, ERFINV_YTOL, ERFINV_MAXITER, &err);
    }
    gsl_set_error_handler(old_handler);

    if (err.code!=0 && err.code!=MAX_ITER_EXCEEDED) {
      fprintf(stderr, "supernewton root-finder returned %d error code\n", err.code);
      exit(err.code);
    }

    erfinv_niter=err.niter;
    erfinv_nbis=err.nbis;

    return y;

  }


  //unit tests:
  int test_erfinv(double miny, double maxy, int nx, int complementary) {
    double x[nx];
    double y[nx];
    double xcalc[nx];
    double ycalc[nx];
    double dx[nx];
    double dy[nx];
    int tniter=0;
    int tnbis=0;
    double (*func1) (double);
    double (*func2) (double);

    if (complementary==0) {
      printf("Transitions to complementary error function at: %g=erf(%g)\n", ERFINV_XTRAN, ERFINV_YTRAN);
      func1=&erfinv;
      func2=&gsl_sf_erf;
    } else {
      printf("Transitions to complementary error function at: %g=erfc(%g)\n", ERFCINV_XTRAN, ERFCINV_YTRAN);
      func1=&erfcinv;
      func2=&gsl_sf_erfc;
    }


    printf("y           y2=erfinv(x) x=erf(y)   x2=erf(y2)  dy          dx          niter nbis\n");  

    for (int i=0; i<nx; i++) {
      y[i]=(maxy-miny)*i/(nx-1)+miny;
      //x[i]=gsl_sf_erf(y[i]);
      x[i]=(*func2)(y[i]);
      //ycalc[i]=erfinv(x[i]);
      ycalc[i]=(*func1)(x[i]);
      //xcalc[i]=gsl_sf_erf(ycalc[i]);
      xcalc[i]=(*func2)(ycalc[i]);
      dx[i]=x[i]-xcalc[i];
      dy[i]=y[i]-ycalc[i];
      printf("%12.6g%12.6g%12.6g%12.6g%12.6g%12.6g%6d%4d\n", y[i], ycalc[i], x[i], xcalc[i], dy[i], dx[i], erfinv_niter, erfinv_nbis);
      tniter+=erfinv_niter;
      tnbis+=erfinv_nbis;
    }
    printf("\n");
    printf("                                                                        %6d%5d\n", tniter, tnbis);  

  }


} //end namespace libpetey

