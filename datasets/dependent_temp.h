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
        virtual errtype cel_1d(dtype new_data, sub_1d_type sub);
        virtual errtype cel(dtype new_data, ind_type index1);
        virtual errtype cel(dtype new_data, ind_type index1, ind_type index2);
        virtual errtype cel(dtype new_data, ind_type index1, ind_type index2, 
		    ind_type index3);
        virtual errtype cel(dtype new_data, ind_type index1, ind_type index2, 
		    ind_type index3, ind_type index4);
        virtual errtype cel(dtype new_data, ind_type *indices);

        virtual errtype insert(rank_type r, ind_type index, ind_type ni);	//inserts slab at specified rank
        virtual errtype del(rank_type r, ind_type index, ind_type ni);	//deletes slab at specified rank

        //extract array elements:
        virtual errtype get_1d(dtype &value, sub_1d_type sub);
        virtual errtype get(dtype &value, ind_type *indices);
        virtual errtype get(dtype &value, ind_type index1);
        virtual errtype get(dtype &value, ind_type index1, ind_type index2);
        virtual errtype get(dtype &value, ind_type index1, ind_type index2, ind_type index3);
        virtual errtype get(dtype &value, ind_type index1, ind_type index2, 
		    ind_type index3, ind_type index4);

        virtual errtype interpol(dtype &value, interpol_index *indices);
        virtual errtype interpol(dtype &value, interpol_index index1);
        virtual errtype interpol(dtype &value, interpol_index index1, interpol_index index2);
        virtual errtype interpol(dtype &value, interpol_index index1, interpol_index index2, interpol_index index3);
        virtual errtype interpol(dtype &value, interpol_index index1, interpol_index index2,
    		interpol_index index3, interpol_index index4);

        //returns a dataset whose dependents are the union of the two operands:
        dependent<dtype> * joint(dependent_dataset *other);

        //for iterating over each dependent:
        virtual errtype iter_get_next(ind_type *indices, rank_type &rank_cur);

        //get interpolation coefficients:
        //virtual dependent<interpol_index> *intcoeff(stype *s);

        //interpolate based on new datasets:
        //virtual errtype interpol(dependent<dtype> *result);
        //virtual dependent<dtype> * interpolate(dependent<interpol_index> **c);
        //virtual errtype interpolate(dependent<interpol_index> **c, dependent<dtype> *result);

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

    };

  } //end namespace datasets
} //end namespace libpetey

#endif

