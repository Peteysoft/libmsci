#include "gsl/gsl_rng.h"
#include "gsl/gsl_interp.h"

#define NSAMPLE_DEFAULT 1000

namespace libagf {

template <class real>
class sample_class2_obj {
  protected:
    double *xs;          //location of spine
    double *ys;
    int ns;
    real width;        //nominal width

    //cubic spline fit:
    double *s;
    gsl_interp * xinterp;
    gsl_interp * yinterp;
    gsl_interp_accel * saccel;
    //real *ddx;
    //real *ddy;

    //spline constraints (boundary conditions):
    real ddx0, ddxf;
    real ddy0, ddyf;

    //long idum;         //random number seed
    gsl_rng *rann;      //gsl random number generator

    void set_spine();   //sets the locations of the spine

    //to avoid duplicate code:
    real pdf(real x, real y, real *dpdx, real *dpdy, real *Lap, int nsamples);

  public:
    sample_class2_obj(real w=0.1);
    ~sample_class2_obj();

    void sample(real &x, real &y);
    //why isn't nsamples stored in the object class??
    real pdf(real x, real y, int nsamples=NSAMPLE_DEFAULT);
    real pdf(real x, real y, real &dpdx, real &dpdy, int nsamples=NSAMPLE_DEFAULT);
    real pdf(real x, real y, real &Lap, int nsamples=NSAMPLE_DEFAULT);

};

}

//typedef sample_class2_obj<real_a> sample_class1_obj;

