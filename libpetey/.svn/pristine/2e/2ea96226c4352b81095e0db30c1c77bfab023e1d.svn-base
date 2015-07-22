#ifndef DEPENDENT_INTC__H
#define DEPENDENT_INTC__H 1

#include <stdio.h>

#include "simple_temp.h"
#include "dependent_dataset.h"

namespace libpetey {
namespace datasets {

typedef simple_dataset stype;

//dummy dataset class that stores no data, but is simply used to
//calculate interpolation coefficients:

class dependent_intc:public dependent_dataset {
  public:
    dependent_intc(stype **s, rank_type nd);

    dependent_intc(stype *s1);
    dependent_intc(stype *s1, stype *s2);
    dependent_intc(stype *s1, stype *s2, stype *s3);
    dependent_intc(stype *s1, stype *s2, stype *s3, stype *s4);

    errtype interpol_coeff(interpol_index *indices, sub_1d_type *subscripts, double *coeffs);

    //these should never be used, but we must have them in order for compilation
    //to take place:
    virtual long read(FILE *fs);
    virtual long write(FILE *fs);

    virtual void print_meta();

    virtual errtype insert(rank_type r, ind_type index, ind_type ni);
    virtual errtype del(rank_type r, ind_type index, ind_type ni);
};

inline long dependent_intc::read(FILE *fs) {
  return 0;
}

inline long dependent_intc::write(FILE *fs) {
  return 0;
}

inline errtype dependent_intc::insert(rank_type r, ind_type index, ind_type ni) {
  return NO_DATA;
}

inline errtype dependent_intc::del(rank_type r, ind_type index, ind_type ni) {
  return NO_DATA;
}

} //end namespace datasets
} //end namespace libpetey

#endif

