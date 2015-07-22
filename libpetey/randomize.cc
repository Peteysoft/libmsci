
//provides a simple and uniform interface for random number generation...
//#include <sys/timeb.h>
#include <sys/time.h>
#include <unistd.h>

#include "peteys_tmpl_lib.h"
#include "randomize.h"

namespace libpetey {

  long libpetey_ran_seed;
  gsl_rng *libpetey_ran;

  unsigned long seed_from_clock() {
    unsigned long seed;
    timeval now;
    pid_t pid;

    //timezone tz;
    //tz.tz_minuteswest=0;
    gettimeofday(&now, NULL);
    pid=getpid();
    //seed= (unsigned long) now.time + ((unsigned long) now.millitm) << 7;
    seed=now.tv_usec*10000+pid;
    //printf("micro-s=%d; pid=%d\n", now.tv_usec, pid);
    //printf("seed = %d\n", seed);
    return seed;
  }

  void ran_init(const gsl_rng_type *rt) {
    libpetey_ran=gsl_rng_alloc(rt);
    gsl_rng_set(libpetey_ran, seed_from_clock());
  }

  long *randomize(long n) {
    double *flt;
    long *ind;

    flt=new double[n];
    for (long i=0; i<n; i++) {
      flt[i]=ranu();
      //printf("%g\n", flt[i]);
    }
    ind=heapsort(flt, n);

    delete [] flt;
    return ind;
  }

  //faster version of randomize O(n):
  long *randomize_f(long n) {
    long *ind;
    long *rind;
    long sub;

    ind=new long[n];
    for (long i=0; i<n; i++) ind[i]=i;

    rind=new long[n];
    for (long i=n; i>0; i--) {
      long swp;
      sub=i*ranu();
      //inclusive or non-inclusive? let's assume inclusive:
      if (sub==i) sub=i-1;
      rind[n-i]=ind[sub];
      swp=ind[i-1];
      ind[i-1]=ind[sub];
      ind[sub]=swp;
    }
    delete [] ind;

    return rind;
  }

  //not really necessary, since one should only call "ran_init"
  //once in a given program...
  void ran_end() {
    gsl_rng_free(libpetey_ran);
  }

} //end namespace libpetey

