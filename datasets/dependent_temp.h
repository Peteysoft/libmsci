#ifndef DEPENDENT_TEMP_INCLUDED
#define DEPENDENT_TEMP_INCLUDED

#include "dependent_dataset.h"
//#include "composite_dataset.h"

//template <class dtype>
//class dependent_swap<dtype>;

namespace libpetey {
namespace datasets {

typedef simple_dataset stype;

//simpler version of the "dependent" dataset with no messing
//around with linked lists etc.  Any new values are simply
//directly inserted into the array.

template <class dtype>
class dependent: public dependent_dataset {
  //friend dependent_swap<dtype>;
  private:
    void set_type();
    errtype init(stype **s, rank_type ndep);
  protected:
    dtype *data;
    dtype missing;

  public:
    virtual ~dependent();
    dependent();

    dependent(stype **s, rank_type ndep);
    dependent(stype *s1);
    dependent(stype *s1, stype *s2);
    dependent(stype *s1, stype *s2, stype *s3);
    dependent(stype *s1, stype *s2, stype *s3, stype *s4);

    //read/write from a binary file:
    virtual long read(FILE *fileptr);
    virtual long write(FILE *fileptr);

    virtual void print_meta();

    //change the specified elements:
    errtype cel_1d(dtype new_data, sub_1d_type sub);
    errtype cel(dtype new_data, ind_type index1);
    errtype cel(dtype new_data, ind_type index1, ind_type index2);
    errtype cel(dtype new_data, ind_type index1, ind_type index2, 
		    ind_type index3);
    errtype cel(dtype new_data, ind_type index1, ind_type index2, 
		    ind_type index3, ind_type index4);
    errtype cel(dtype new_data, ind_type *indices);
    virtual errtype insert(rank_type r, ind_type index, ind_type ni);	//inserts slab at specified rank
    virtual errtype del(rank_type r, ind_type index, ind_type ni);	//deletes slab at specified rank

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
    errtype interpol(dtype &value, interpol_index index1, interpol_index index2, interpol_index index3);
    errtype interpol(dtype &value, interpol_index index1, interpol_index index2,
    		interpol_index index3, interpol_index index4);

    //returns a dataset whose dependents are the union of the two operands:
    dependent<dtype> * joint(dependent_dataset *other);

    //for iterating over each dependent:
    errtype iter_get_next(ind_type *indices, rank_type &rank_cur);

    //get interpolation coefficients:
    dependent<interpol_index> *intcoeff(stype *s);

    //interpolate based on new datasets:
    errtype interpol(dependent<dtype> *result);

    dependent<dtype> * interpolate(dependent<interpol_index> **c);
    errtype interpolate(dependent<interpol_index> **c, dependent<dtype> *result);

    errtype interpol_old(dtype &value, interpol_index *indices);
    errtype interpol_debug(dtype &value, interpol_index *indices);
    
    //copy constructor:
    dependent(const dependent<dtype> &old);

    //preload data values:
    virtual errtype preload(dtype *new_data, sub_1d_type n, sub_1d_type offset);
    virtual errtype read_chunk(FILE *fileptr, sub_1d_type offset, sub_1d_type n);

    //assignment operator:
    dependent<dtype> operator = (const dependent<dtype> &old);

/*
    //conversion routines:
    template <class dtype2>
    virtual operator dependent<dtype2> ();
*/
/*
    //binary operators:
    virtual dataset *add(const dataset &other);
    template <class dtype2>
    dependent<dtype> *add (const dependent<dtype2> &other);

    virtual dataset *subtract(const dataset &other);
    template <class dtype2>
    dependent<dtype> *subtract (const dependent<dtype2> &other);

    virtual dataset *multiply(const dataset &other);
    template <class dtype2>
    dependent<dtype> *multiply (const dependent<dtype2> &other);

    virtual dataset *divide(const dataset &other);
    template <class dtype2>
    dependent<dtype> *divide (const dependent<dtype2> &other);
*/

};

template <class dtype>
inline errtype dependent<dtype>::cel_1d(dtype new_data, sub_1d_type sub) {
	if (sub < 0 || sub >= n_data) return SUBSCRIPT_OUT_OF_RANGE;

	data[sub]=new_data;

	return 0;

}

template <class dtype>
inline errtype dependent<dtype>::cel(dtype data_el, ind_type *indices) {
	sub_1d_type sub;

	sub=calc_sub(indices);
	if (sub < 0) return sub;
	data[sub]=data_el;
	return 0;
}

template <class dtype>
inline errtype dependent<dtype>::cel(dtype new_data, ind_type index1) {
	if (rank != 1) return -1;
	if (index1 < 0 || index1 >= dim[0]) return INDEX1_OUT_OF_RANGE;

	data[index1]=new_data;

	return 0;

}

template <class dtype>
inline errtype dependent<dtype>::cel(dtype new_data, ind_type index1, ind_type index2) {
	sub_1d_type sub;

        if (rank != 2) return RANK_ERROR;
	if (index1<0 || index1 >= dim[0]) return INDEX1_OUT_OF_RANGE;
	if (index2<0 || index2 >= dim[1]) return INDEX2_OUT_OF_RANGE;

//	cout << "Changing element at: " << index1 << " " << index2 <<
//			" to: " << new_data << "\n";

	sub=index1+index2*dim[0];
	data[sub]=new_data;

	return 0;
}

template <class dtype>
inline errtype dependent<dtype>::cel(dtype new_data, ind_type index1, ind_type index2,
		ind_type index3) {

  if (rank != 3) return RANK_ERROR;
  if (index1<0 || index1 >= dim[0]) return INDEX1_OUT_OF_RANGE;
  if (index2<0 || index2 >= dim[1]) return INDEX2_OUT_OF_RANGE;
  if (index3<0 || index3 >= dim[2]) return INDEX3_OUT_OF_RANGE;

  data[index1+dim[0]*(index2+index3*dim[1])]=new_data;

  return 0;
}

template <class dtype>
inline errtype dependent<dtype>::cel(dtype new_data, ind_type index1, ind_type index2,
		ind_type index3, ind_type index4) {
  sub_1d_type sub;

  if (rank != 4) return RANK_ERROR;
  if (index1<0 || index1 >= dim[0]) return INDEX1_OUT_OF_RANGE;
  if (index2<0 || index2 >= dim[1]) return INDEX2_OUT_OF_RANGE;
  if (index3<0 || index3 >= dim[2]) return INDEX3_OUT_OF_RANGE;
  if (index4<0 || index4 >= dim[3]) return INDEX4_OUT_OF_RANGE;

  sub=index1+dim[0]*(index2+dim[1]*(index3+index4*dim[2]));
  data[sub]=new_data;

  return 0;
}

template <class dtype>
inline errtype dependent<dtype>::get_1d(dtype &value, sub_1d_type sub) {
  if (sub<0 || sub >= n_data) return SUBSCRIPT_OUT_OF_RANGE;
  value=data[sub];

  return sub;
}

template <class dtype>
inline errtype dependent<dtype>::get(dtype &value, ind_type *indices) {
	sub_1d_type sub;

	sub=calc_sub(indices);
	if (sub < 0) return sub;

	value=data[sub];

	return sub;
}

template <class dtype>
inline errtype dependent<dtype>::get(dtype &value, ind_type index1) {
  if (rank != 1L) return RANK_ERROR;
  if (index1<0 || index1 >= dim[0]) return INDEX1_OUT_OF_RANGE;
  value=data[index1];

  return index1;
}

template <class dtype>
inline errtype dependent<dtype>::get(dtype &value, ind_type index1, ind_type index2) {
  sub_1d_type sub;

  if (rank != 2L) return RANK_ERROR;
  if (index1<0 || index1 >= dim[0]) return INDEX1_OUT_OF_RANGE;
  if (index2<0 || index2 >= dim[1]) return INDEX2_OUT_OF_RANGE;
  sub=index1+dim[0]*index2;
  value=data[sub];

  return sub;
}

template <class dtype>
inline errtype dependent<dtype>::get(dtype &value, ind_type index1, ind_type index2, 
			ind_type index3) {
  sub_1d_type sub;

  if (rank != 3L) return RANK_ERROR;
  if (index1<0 || index1 >= dim[0]) return INDEX1_OUT_OF_RANGE;
  if (index2<0 || index2 >= dim[1]) return INDEX2_OUT_OF_RANGE;
  if (index3<0 || index3 >= dim[2]) return INDEX3_OUT_OF_RANGE;
  sub=index1+dim[0]*(index2+dim[1]*index3);
  value=data[sub];

  return sub;
}

template <class dtype>
inline errtype dependent<dtype>::get(dtype &value, ind_type index1, ind_type index2,
			ind_type index3, ind_type index4) {
  sub_1d_type sub;

  if (rank != 4L) return RANK_ERROR;
  if (index1<0 || index1 >= dim[0]) return INDEX1_OUT_OF_RANGE;
  if (index2<0 || index2 >= dim[1]) return INDEX2_OUT_OF_RANGE;
  if (index3<0 || index3 >= dim[2]) return INDEX3_OUT_OF_RANGE;
  if (index4<0 || index4 >= dim[3]) return INDEX4_OUT_OF_RANGE;
  sub=index1+dim[0]*(index2+dim[1]*(index3+index4*dim[2]));
  value=data[sub];

  return sub;
}

template <class dtype>
inline errtype dependent<dtype>::interpol(dtype &value, interpol_index index1) {
  if (rank != 1) return RANK_ERROR;
  return interpol(value, &index1);
}

template <class dtype>
inline errtype dependent<dtype>::interpol(dtype &value, interpol_index index1, 
			interpol_index index2) {
  interpol_index indices[2];

  if (rank != 2) return RANK_ERROR;
  indices[0]=index1;
  indices[1]=index2;
  return interpol(value, indices);
}

template <class dtype>
inline errtype dependent<dtype>::interpol(dtype &value, interpol_index index1, 
		interpol_index index2, interpol_index index3) {
  interpol_index indices[3];

  if (rank != 3) return RANK_ERROR;
  indices[0]=index1;
  indices[1]=index2;
  indices[2]=index3;
  return interpol(value, indices);
}

template <class dtype>
inline errtype dependent<dtype>::interpol(dtype &value, interpol_index index1, 
		interpol_index index2, interpol_index index3, interpol_index index4) {
  interpol_index indices[4];

  if (rank != 4) return RANK_ERROR;
  indices[0]=index1;
  indices[1]=index2;
  indices[2]=index3;
  indices[3]=index4;

  return interpol(value, indices);
}

} //end namespace datasets
} //end namespace libpetey

#endif

