#ifndef CTRAJ_VFIELD_2D_H
#define CTRAJ_VFIELD_2D_H 1

#include "time_class.h"
#include "simple_temp.h"
#include "dependent_swap.h"
#include "composite_dataset.h"

#include "ctraj_vfield_base.h"

namespace ctraj {
  using namespace libpetey;
  using namespace datasets;

  template <class real>
  class ctraj_vfield_2d:public ctraj_vfield_single_domain<real> {

    protected:
      composite_dataset all;
      simple<time_class> *t;
      simple<real> *x;
      simple<real> *y;
      simple<real> *z;
      dependent_swap<real> *U;
      dependent_swap<real> *V;
      FILE *fstream;
    public:
      ctraj_vfield_2d();
      virtual ~ctraj_vfield_2d();

      int init(char *name, int64_t pagesize);

      virtual void help(FILE *fs);
      virtual int setup(int argc, char **argv);

      virtual int v(int32_t domain, double tind, real *x, real *v);

      virtual int32_t ndim();
      virtual double maxt();

      virtual double get_tind(char *date);
      virtual int get_t(double tind, char *date);

      double get_tind(time_class date);
      time_class get_t(double tind);

  };

}

#endif

