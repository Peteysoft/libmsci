#ifndef ctraj_METRIC_BASE__H
#define ctraj_METRIC_BASE__H 1

namespace ctraj {

  template <class real>
  class metric_base {
    public:
      virtual real ds2(real *x1, real *x2)=0;
      virtual void mcoef2(real *x, real *c)=0;
      virtual real gcurv(real *x)=0;
      virtual int ndim()=0;

  };

  template <class real>
  class Cart_t:public metric_base<real> {
    protected:
      int n;
    public:
      Cart_t(int n);

      virtual real ds2(real *x1, real *x2);
      virtual void mcoef2(real *x, real *c);
      virtual real gcurv(real *x);
      virtual int ndim();

  };

} //end namespace ctraj

#endif

