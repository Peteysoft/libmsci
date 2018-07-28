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

	//set size of page held in RAM:
        void set_page_size(sub_1d_type size,	//page size in data elements
			int mb=0);		//make pages "snap to" grid

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
        virtual errtype cel_1d(dtype new_data, sub_1d_type sub);
        virtual errtype cel(dtype new_data, ind_type index1);
        virtual errtype cel(dtype new_data, ind_type index1, ind_type index2);
        virtual errtype cel(dtype new_data, ind_type index1, ind_type index2, ind_type index3);
        virtual errtype cel(dtype new_data, ind_type index1, ind_type index2, 
		    ind_type index3, ind_type index4);
        virtual errtype cel(dtype new_data, ind_type *indices);

        virtual errtype insert(rank_type r, ind_type index,
		   ind_type ni);		//inserts slab at specified rank
        virtual errtype del(rank_type r, ind_type index, 
		    ind_type ni);		//deletes slab at specified rank

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
        virtual errtype interpol(dtype &value, interpol_index index1, interpol_index index2, 
		    interpol_index index3);
        virtual errtype interpol(dtype &value, interpol_index index1, interpol_index index2,
    		interpol_index index3, interpol_index index4);

        //copy constructor:
        dependent_swap(const dependent_swap<dtype> &old);

        //preload data values:
        virtual errtype preload(dtype *new_data, sub_1d_type n, sub_1d_type offset);
        //virtual long preload(dependent<dtype> *new_data, sub_1d_type offset);
        virtual errtype read_chunk(FILE *fileptr, sub_1d_type offset, sub_1d_type n);
    };


  } //end namespace datasets
} //end namespace libpetey

#endif

