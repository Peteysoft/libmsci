
#include "gsl/gsl_rng.h"

namespace libagf {

template <class real>
class sample_class1_obj {
  protected:
    real x0;          //location of centre
    real y0;

    real w1;		//width in x-dir (before rotation)
    real w2;		//width in y-dir (       "       )

    //rotation angle as unit vector:
    real cosang;
    real sinang;

//    long idum;         //random number seed
    gsl_rng *rann;	//gsl random number generator

  public:
    sample_class1_obj();
    ~sample_class1_obj();
    sample_class1_obj(real x0u, real y0u, real w1u, real w2u, real angle);

    void sample(real &x, real &y);
    real pdf(real x, real y);

    //returns the derivatives of the conditional probabilities:
    real pdf(real x, real y, real &dpdx, real &dpdy);

};

}

//typedef sample_class1_obj<real_a> sample_class1_obj;

