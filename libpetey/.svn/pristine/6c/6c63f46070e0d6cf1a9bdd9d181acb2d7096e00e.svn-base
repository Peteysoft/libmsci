#include "dependent_dataset.h"

namespace libpetey {
namespace datasets {

errtype INDEX_OUT_OF_RANGE[20] = {
        INDEX1_OUT_OF_RANGE,
        INDEX2_OUT_OF_RANGE,
        INDEX3_OUT_OF_RANGE,
        INDEX4_OUT_OF_RANGE,
        INDEX5_OUT_OF_RANGE,
        INDEX6_OUT_OF_RANGE,
        INDEX7_OUT_OF_RANGE,
        INDEX8_OUT_OF_RANGE,
        INDEX9_OUT_OF_RANGE,
        INDEX10_OUT_OF_RANGE,
        INDEX11_OUT_OF_RANGE,
        INDEX12_OUT_OF_RANGE,
        INDEX13_OUT_OF_RANGE,
        INDEX14_OUT_OF_RANGE,
        INDEX15_OUT_OF_RANGE,
        INDEX16_OUT_OF_RANGE,
        INDEX17_OUT_OF_RANGE,
        INDEX18_OUT_OF_RANGE,
        INDEX19_OUT_OF_RANGE,
        INDEX20_OUT_OF_RANGE};

dependent_dataset::~dependent_dataset() {
  rank_type i;
  for (i=0;i<rank;i++) {
    if (dependencies[i]!=NULL) dependencies[i]->del_dependent(this);
  }
  if (dependencies != NULL) delete [] dependencies;
  delete [] dim;
}

void dependent_dataset::indices(sub_1d_type oned, ind_type *indices) {
  sub_1d_type div;
  sub_1d_type remainder;

  remainder=oned;
  for (rank_type i=0; i<rank-1; i++) {
    div=remainder/dim[i];
    indices[i]=remainder-div*dim[i];
    remainder=div;
    //indices[i]=remainder % dim[i];
    //remainder=remainder/dim[i];
  }
  indices[rank-1]=remainder;
}

dependent_dataset::dependent_dataset() {
  type=DEPENDENT;
  rank=0;
  dim=NULL;
  dependencies=NULL;
}

} //end namespace datasets
} //end namespace libpetey

