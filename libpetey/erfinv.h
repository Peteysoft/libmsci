#include <gsl/gsl_sf_erf.h>

namespace libpetey {

  const double ERFINV_YTRAN=6.;		//larger than this we invert the complementary error func.
  const double ERFCINV_YTRAN=3.;	//larger than this we invert the complementary error func.
 

  int erfinv_niter;			//sets number of iterations after each call
  int erfinv_nbis;			//number of bisection steps

  double erfinv(double x);

  double erfcinv(double x);

  int test_erfinv(double miny, double maxy, int nx, int complementary);

} //end namespace libpetey

