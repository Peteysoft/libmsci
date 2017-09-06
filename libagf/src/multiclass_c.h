#include "multiclass_hier.h"

namespace libagf {
  //multi-class classification adapted for continuum retrievals:
  template <typename real, typename cls_t>
  class multiclass_c:public multiclass_hier<real, cls_t> {
    protected:
      real *y0;		//abscissa
    public:
      multiclass_c(const char *file, int tp);
      virtual ~multiclass_c();

      //continuum retrieval:
      real ret(real *x, real &err);
		

  };

}

