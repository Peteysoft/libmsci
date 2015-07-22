#ifndef CTRAJ_VFIELD_ANAL_H
#define CTRAJ_VFIELD_ANAL_H 1

#include "ctraj_vfield_base.h"

namespace ctraj {

  template <class real>
  class ctraj_vfield_anal:public ctraj_vfield_single_domain<real> {
    protected:
      int32_t ndims;
      void *handle;
      int (* vfield) (double t, real *x, real *v, void *param);
      void *param;
      int32_t nparam;

    public:
      ctraj_vfield_anal(int32_t nd);
      virtual ~ctraj_vfield_anal();

      int init(char *file, char *sym, void *p);
      virtual void help(FILE *fs);
      virtual int setup(int argc, char **argv);

      virtual int v(int32_t domain, double tind, real *x, real *v);

      virtual int32_t ndim();
      virtual double maxt();

  };
} //end namespace ctraj

#endif

