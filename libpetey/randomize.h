#ifndef RANDOMIZE_H
#define RANDOMIZE_H 1

//provides a simple and uniform interface for random number generation...
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

namespace libpetey {

  extern gsl_rng *libpetey_ran;

  unsigned long seed_from_clock();

  void ran_init(const gsl_rng_type *rt=gsl_rng_mt19937);

  //uniform deviates:
  inline double ranu() {
    return gsl_rng_uniform(libpetey_ran);
  }

  //gaussian deviates:
  inline double rang() {
    return gsl_ran_ugaussian(libpetey_ran);
  }

  //for permuting arrays:
  long *randomize(long n);

  void ran_end();

} //end namespace libpetey

#endif

