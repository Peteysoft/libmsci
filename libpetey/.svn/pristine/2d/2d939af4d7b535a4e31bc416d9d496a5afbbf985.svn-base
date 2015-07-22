#ifndef DEPENDENT_DATASET_DEFINED
#define DEPENDENT_DATASET_DEFINED

#include "dataset.h"

namespace libpetey {
  namespace datasets {
    class dependent_dataset;
  }
}

#include "simple_dataset.h"
//class simple_dataset;
//#include "error.h"

namespace libpetey {
namespace datasets {

class composite_dataset;

#define stype simple_dataset

//simpler version of the "dependent" dataset with no messing
//around with linked lists etc.  Any new values are simply
//directly inserted into the array.

class dependent_dataset: public dataset {
  friend class composite_dataset;
  protected:
//    dtype *data;
    rank_type rank;
    ind_type *dim;
    stype **dependencies;
  public:
    dependent_dataset();
    virtual ~dependent_dataset()=0;

    virtual long read(FILE *fileptr)=0;
    virtual long write(FILE *fileptr)=0;
    virtual void print_meta()=0;

    virtual errtype insert(rank_type r, ind_type index, ind_type ni)=0;	//inserts slab at specified rank
    virtual errtype del(rank_type r, ind_type index, ind_type ni)=0;	//deletes slab at specified rank

    //calculates the one-dim. index into arry:
    sub_1d_type calc_sub(ind_type *indices);
    sub_1d_type calc_sub(ind_type index1);
    sub_1d_type calc_sub(ind_type index1, ind_type index2);
    sub_1d_type calc_sub(ind_type index1, ind_type index2, ind_type index3);
    sub_1d_type calc_sub(ind_type index1, ind_type index2, 
		    ind_type index3, ind_type index4);
    
    //given a one dimensional index, calculates the n-dim. index:
    void indices(sub_1d_type oned, ind_type *indices);

    //extract array elements:

    //informational routines:
    rank_type get_rank();
    void get_dim(ind_type *indices);
    sub_1d_type size_of();

/*    //conversion routines:
    template <class dtype>
    virtual operator dependent<dtype> ()=0;
*/

};

//inline everything to keep it all in one file:

inline rank_type dependent_dataset::get_rank() {
	return rank;
}

inline void dependent_dataset::get_dim(ind_type *indices) {
	for (long i=0; i<rank; i++) indices[i]=dim[i];
}

inline sub_1d_type dependent_dataset::calc_sub(ind_type *indices) {
	sub_1d_type sub;
	ind_type i, j;

	sub=0;
	j=1;
	for (i=0;i<rank;i++) {
		if (indices[i]<0 || indices[i]>=dim[i]) return INDEX_OUT_OF_RANGE[i];
		sub=sub+j*indices[i];
		j=j*dim[i];
	}
	return sub;
}

inline sub_1d_type dependent_dataset::calc_sub(ind_type index1) {
  return index1;
}

inline sub_1d_type dependent_dataset::calc_sub(ind_type index1, ind_type index2) {
  return index1+dim[0]*index2;
}

inline sub_1d_type dependent_dataset::calc_sub(ind_type index1, ind_type index2, 
		ind_type index3) {
  return index1+dim[0]*(index2+dim[1]*index3);
}


inline sub_1d_type dependent_dataset::calc_sub(ind_type index1, ind_type index2,
		ind_type index3, ind_type index4) {
  return index1+dim[0]*(index2+dim[1]*(index3+dim[2]*index4));
}

} //end namespace datasets
} //end namespace libpetey

#endif
