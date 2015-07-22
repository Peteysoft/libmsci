#ifndef DEPENDENT_SWAP_INCLUDED
#define DEPENDENT_SWAP_INCLUDED

#include <stdio.h>
#include "dependent_temp.h"

namespace libpetey {
namespace datasets {

template <class dtype>
class dependent_swap: public dependent<dtype> {
  private:
    void set_type();
    errtype init(stype **s, rank_type ndep);

  protected:
    FILE *swap;			//stream pointing to swap file
    sub_1d_type swap_start;		//byte in file where data starts
    sub_1d_type page_start;	//index where the current loaded page starts
    sub_1d_type page_size;	//number of data elements comprising the page
    char fchange;		//file has been changed, current page must be re-loaded
    char dchange;		//data in ram has changed, page must be written
    sub_1d_type check_page(sub_1d_type ind, char debug_flag=0);	//corrects 1d index, checks for page faults

  public:
    virtual ~dependent_swap();
    dependent_swap();

    void set_page_size(sub_1d_type size, int mb=0);

    dependent_swap(stype **s, rank_type ndep);
    dependent_swap(stype *s1);
    dependent_swap(stype *s1, stype *s2);
    dependent_swap(stype *s1, stype *s2, stype *s3);
    dependent_swap(stype *s1, stype *s2, stype *s3, stype *s4);

    //read/write from a binary file:
    virtual long read(FILE *fileptr);
    virtual long write(FILE *fileptr);

    virtual void print_meta();

    //change the specified elements:
    errtype cel_1d(dtype new_data, sub_1d_type sub);
    errtype cel(dtype new_data, ind_type index1);
    errtype cel(dtype new_data, ind_type index1, ind_type index2);
    errtype cel(dtype new_data, ind_type index1, ind_type index2, ind_type index3);
    errtype cel(dtype new_data, ind_type index1, ind_type index2, 
		    ind_type index3, ind_type index4);
    errtype cel(dtype new_data, ind_type *indices);

    virtual errtype insert(rank_type r, ind_type index,
		   ind_type ni);		//inserts slab at specified rank
    virtual errtype del(rank_type r, ind_type index, 
		    ind_type ni);		//deletes slab at specified rank

    //extract array elements:
    errtype get_1d(dtype &value, sub_1d_type sub);
    errtype get(dtype &value, ind_type *indices);
    errtype get(dtype &value, ind_type index1);
    errtype get(dtype &value, ind_type index1, ind_type index2);
    errtype get(dtype &value, ind_type index1, ind_type index2, ind_type index3);
    errtype get(dtype &value, ind_type index1, ind_type index2, 
		    ind_type index3, ind_type index4);

    errtype interpol(dtype &value, interpol_index *indices);
    errtype interpol(dtype &value, interpol_index index1);
    errtype interpol(dtype &value, interpol_index index1, interpol_index index2);
    errtype interpol(dtype &value, interpol_index index1, interpol_index index2, 
		    interpol_index index3);
    errtype interpol(dtype &value, interpol_index index1, interpol_index index2,
    		interpol_index index3, interpol_index index4);

    //copy constructor:
    dependent_swap(const dependent_swap<dtype> &old);

    //preload data values:
    virtual errtype preload(dtype *new_data, sub_1d_type n, sub_1d_type offset);
    //virtual long preload(dependent<dtype> *new_data, sub_1d_type offset);
    virtual errtype read_chunk(FILE *fileptr, sub_1d_type offset, sub_1d_type n);
};

template <class dtype>
inline errtype dependent_swap<dtype>::cel_1d(dtype new_data, sub_1d_type sub) {
  sub_1d_type new_sub;

  new_sub=check_page(sub);

  if (new_sub != -1) {
    dchange=1;
    this->data[new_sub]=new_data;
  }

  return new_sub;
}

template <class dtype>
inline errtype dependent_swap<dtype>::cel(dtype new_data, ind_type index1) {
  if (this->rank != 1) return RANK_ERROR;
  return cel_1d(new_data, index1);
}

template <class dtype>
inline errtype dependent_swap<dtype>::cel(dtype new_data, ind_type index1, 
		ind_type index2) {
  if (this->rank != 2) return RANK_ERROR;
  sub_1d_type sub=index1+index2*this->dim[0];
  return cel_1d(new_data, sub);
}

template <class dtype>
inline errtype dependent_swap<dtype>::cel(dtype new_data, ind_type index1, 
		ind_type index2, ind_type index3) {
  if (this->rank != 3) return RANK_ERROR;
  sub_1d_type sub=index1+this->dim[0]*(index2+this->dim[1]*index3);
  return cel_1d(new_data, sub);
}

template <class dtype>
inline errtype dependent_swap<dtype>::cel(dtype new_data, ind_type index1, 
		ind_type index2, ind_type index3, ind_type index4) {
  if (this->rank != 4) return RANK_ERROR;
  sub_1d_type sub=index1+this->dim[0]*(index2+this->dim[1]*(index3+this->dim[2]*index4));
  //printf("%f\n", new_data);
  return cel_1d(new_data, sub);
}

template <class dtype>
inline errtype dependent_swap<dtype>::cel(dtype new_data, ind_type *indices) {
  sub_1d_type sub;

  sub=this->calc_sub(indices);
  return cel_1d(new_data, sub);
}

template <class dtype>
inline errtype dependent_swap<dtype>::get_1d(dtype &value, sub_1d_type sub) {
  sub_1d_type new_sub;
  new_sub=check_page(sub);
  if (new_sub >= 0) value=this->data[new_sub];
  return new_sub;
}

template <class dtype>
inline errtype dependent_swap<dtype>::get(dtype &value, ind_type *indices) {
  sub_1d_type sub;

  sub=this->calc_sub(indices);
  return get_1d(value, sub);
}

template <class dtype>
inline errtype dependent_swap<dtype>::get(dtype &value, ind_type index1) {
  if (this->rank != 1) return RANK_ERROR;
  return get_1d(value, index1);
}

template <class dtype>
inline errtype dependent_swap<dtype>::get(dtype &value, ind_type index1, 
		ind_type index2) {
  if (this->rank != 2) return RANK_ERROR;
  sub_1d_type sub=index1+index2*this->dim[0];
  return get_1d(value, sub);
}

template <class dtype>
inline errtype dependent_swap<dtype>::get(dtype &value, ind_type index1, 
		ind_type index2, ind_type index3) {
  if (this->rank != 3) return RANK_ERROR;
  sub_1d_type sub=index1+this->dim[0]*(index2+this->dim[1]*index3);
  return get_1d(value, sub);
}

template <class dtype>
inline errtype dependent_swap<dtype>::get(dtype &value, ind_type index1, 
		ind_type index2, ind_type index3, ind_type index4) {
  if (this->rank != 4) return RANK_ERROR;
  sub_1d_type sub=index1+this->dim[0]*(index2+this->dim[1]*(index3+this->dim[2]*index4));
  get_1d(value, sub);
}

template <class dtype>
inline errtype dependent_swap<dtype>::interpol(dtype &value, interpol_index index1) {
  if (this->rank != 1) return RANK_ERROR;
  return interpol(value, &index1);
}

template <class dtype>
inline errtype dependent_swap<dtype>::interpol(dtype &value, interpol_index index1,
	       interpol_index index2) {
  interpol_index indices[2];

  if (this->rank != 2) return RANK_ERROR;
  indices[0]=index1;
  indices[1]=index2;
  return interpol(value, indices);
}

template <class dtype>
inline errtype dependent_swap<dtype>::interpol(dtype &value, interpol_index index1, 
		interpol_index index2, interpol_index index3) {
  interpol_index indices[3];

  if (this->rank != 3) return RANK_ERROR;
  indices[0]=index1;
  indices[1]=index2;
  indices[2]=index3;
  return interpol(value, indices);
}

template <class dtype>
inline errtype dependent_swap<dtype>::interpol(dtype &value, interpol_index index1, 
		interpol_index index2, interpol_index index3, interpol_index index4) {
  interpol_index indices[4];

  if (this->rank != 4) return RANK_ERROR;
  indices[0]=index1;
  indices[1]=index2;
  indices[2]=index3;
  indices[3]=index4;

  return interpol(value, indices);
}

} //end namespace datasets
} //end namespace libpetey

#endif

