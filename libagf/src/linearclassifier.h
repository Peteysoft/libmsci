#ifndef _LIBAGF_LINEARCLASSIFIER__H
#define _LIBAGF_LINEARCLASSIFIER__H

#include "binaryclassifier.h"

namespace libagf {
  template <typename real, typename cls_t>
  class linearclassifier:public binaryclassifier<real, cls_t> {
    protected:
      char **header;		//contains rest of header for liblinear
      real *coef;
      real offset;
    public:
      linearclassifier();
      linearclassifier(char *fname);
      virtual ~linearclassifier();

      real decision(real *x);

      virtual int load(FILE *fs);
      virtual int save(FILE *fs);

  };

}

#endif
