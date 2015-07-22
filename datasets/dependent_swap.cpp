#include "dependent_swap.h"

namespace libpetey {
namespace datasets {

template <class dtype>
dependent_swap<dtype>::dependent_swap() {
  this->rank=0;
  this->dim=NULL;
  this->dependencies=NULL;

  set_type();
  this->data=NULL;

  fchange=0;
  dchange=1;

  swap=NULL;
}

template <class dtype>
dependent_swap<dtype>::~dependent_swap() {
  sub_1d_type nrw;
  size_t nwritten;

  //if data has changed, must be swapped back to the file:
  if (dchange && !fchange) {
    if (page_start+page_size > this->n_data) nrw=this->n_data-page_start;
	    else nrw=page_size;
    fseek(swap, swap_start+page_start*sizeof(dtype), SEEK_SET);
    nwritten=fwrite(this->data, sizeof(dtype), nrw, swap);
    fprintf(datasets_log, "Writing page (%ld bytes) at location %ld (%ld offset)\n", 
		(int64_t) nrw*sizeof(dtype), (int64_t) page_start, (int64_t) (swap_start+page_start*sizeof(dtype)));
    fprintf(datasets_log, "Wrote %ld data elements\n", (int64_t) nwritten);
    //for (long i=0L; i<page_size; i++) if (this->data[i]!=0) printf("%g \n", this->data[i]);
    //for (long i=0L; i<page_size; i++) printf("%g \n", this->data[i]);
  }

}

template <class dtype>
void dependent_swap<dtype>::set_page_size(sub_1d_type size, int mb) {
  ind_type last_dim=this->dim[this->rank-1];
  sub_1d_type goal_size, approx, one_grid;

  if (mb  != 0) goal_size=size/sizeof(dtype);
		else goal_size=50000000/sizeof(dtype);
  //default page size of 50 mb

  if (size==0 || mb != 0) {
    //goal_size=50000000/sizeof(dtype);	//we aim for a page size of ~50 Mb
    one_grid=this->n_data/last_dim;
    approx=(goal_size/one_grid)*one_grid;
    if (last_dim <2) {
      if (this->n_data<goal_size) page_size=this->n_data; else page_size=goal_size;
    } else {
      //use approx 50 Mb or two grids, whichever is bigger:
      if (this->n_data < approx) page_size=this->n_data;
        else if (one_grid*2 > approx) page_size=one_grid*2;
        else page_size=approx;
    }
  } else {
    if (size > this->n_data) page_size=this->n_data; else page_size=size;
  }

  fprintf(datasets_log, "Using page size of %ld data elements (%ld bytes)\n", 
		(int64_t) page_size, (int64_t) page_size*sizeof(dtype));
}

template <class dtype>
sub_1d_type dependent_swap<dtype>::check_page(sub_1d_type ind, char debug_flag) {
  sub_1d_type nrw;
  sub_1d_type new_ind;
  sub_1d_type one_grid=this->n_data/this->dim[this->rank-1];
  int read_flag;
  size_t nwritten;

  if (ind < 0 || ind >= this->n_data) {
    new_ind=SUBSCRIPT_OUT_OF_RANGE;
  } else {
    //if (debug_flag) printf("page_start=%d, ind=%d\n", page_start, ind);
    //if (dchange == 1 && debug_flag) printf("dchange flag set\n");
    read_flag=0;
    if (this->data==NULL) {
      //load up the first page:
      this->data=new dtype[page_size];
      page_start=ind-page_size/10;
      page_start=(page_start/one_grid)*one_grid;
      if (page_start<0) page_start=0;
      read_flag=1;
    }
    new_ind=ind-page_start;

    if (fchange) {
      dchange=0;		//file changes always take precedence over buffer changes
      fchange=0;
      //if the file has changed, read or re-read regardless:
      read_flag=1;
    }

    if (new_ind < 0L || new_ind >= page_size) {
      fprintf(datasets_log, "Page fault: start %ld, size %ld, abs %ld, rel %ld\n", 
		(int64_t) page_start, (int64_t) page_size, (int64_t) ind, (int64_t) new_ind);
      //if buffer has changed, but file has not, then write to file:
      if (dchange && !fchange) {
        if (page_start+page_size > this->n_data) nrw=this->n_data-page_start;
		else nrw=page_size;
        fseek(swap, swap_start+page_start*sizeof(dtype), SEEK_SET);
        nwritten=fwrite(this->data, sizeof(dtype), nrw, swap);
        dchange=0;
        fprintf(datasets_log, "Writing page (%ld bytes) at location %ld (%ld offset)\n", 
			(int64_t) nrw*sizeof(dtype), (int64_t) page_start, 
			(int64_t) (swap_start+page_start*sizeof(dtype)));
	fprintf(datasets_log, "Wrote %ld data elements\n", (int64_t) nwritten);
      }
      //calculate page start:
      if (new_ind<0L) {
        page_start=ind-9*page_size/10;
        page_start=(page_start/one_grid)*one_grid;
        //printf("new_ind=%d, new page start calculated at: %d\n", new_ind, page_start);
        //printf("ind=%d\n", ind);
      } else {
        page_start=ind-page_size/10;
        page_start=(page_start/one_grid)*one_grid;
        //printf("ind=%d\n", ind);
        //printf("new_ind=%d, new page start calculated at: %d\n", new_ind, page_start);
      }
      if (page_start<0) page_start=0;
      read_flag=1;
    }

    //if the read flag is set, read in new data at page_start:
    if (read_flag) {
      if (page_start+page_size > this->n_data) {
        nrw=this->n_data-page_start;
      } else {
        nrw=page_size;
      }
      fseek(swap, swap_start+page_start*sizeof(dtype), SEEK_SET);
      fread(this->data, sizeof(dtype), nrw, swap);
      new_ind=ind-page_start;
      fprintf(datasets_log, "Swapping in new page at location %d\n", page_start);
    }
  }

  //if (debug_flag) printf("new page_start=%d, new_ind=%d\n", page_start, new_ind);

  return new_ind;

}

template <class dtype>
dependent_swap<dtype>::dependent_swap(stype **s, rank_type ndep) {
  init(s, ndep);
}

template <class dtype>
errtype dependent_swap<dtype>::init(stype **s, rank_type ndep) {
  rank_type i;

  set_type();
  swap=NULL;
  this->data=NULL;

  this->rank=ndep;
  this->dim=new ind_type[this->rank];
  this->dependencies=new stype *[this->rank];

  this->n_data=1;
  for (i=0; i<this->rank; i++) {
    this->dependencies[i]=s[i];
    s[i]->add_dependent(this, i);
    this->dim[i]=this->dependencies[i]->nel();
    this->n_data*=this->dim[i];
  }

  set_page_size(0);
  fchange=0;
  dchange=0;
}

template <class dtype>
dependent_swap<dtype>::dependent_swap(stype *s1) {
  init(&s1, 1);
}


template <class dtype>
dependent_swap<dtype>::dependent_swap(stype *s1, stype *s2) {
  stype *s[2]={s1, s2};
  init(s, 2);
}

template <class dtype>
dependent_swap<dtype>::dependent_swap(stype *s1, stype *s2, stype *s3) {
  stype *s[3]={s1, s2, s3};
  init(s, 3);
}

template <class dtype>
dependent_swap<dtype>::dependent_swap(stype *s1, stype *s2, stype *s3, stype *s4) {
  stype *s[4]={s1, s2, s3, s4};
  init(s, 4);
}

template <class dtype>
long dependent_swap<dtype>::read(FILE *fileptr) {
  //set the swap pointer to fileptr and get the position in the file:
  swap=fileptr;
  swap_start=ftell(swap);
  fchange=1;

  //advance the file pointer to the end of the dataset data:
  fseek(swap, this->n_data*sizeof(dtype), SEEK_CUR);

}

template <class dtype>
long dependent_swap<dtype>::write(FILE *fileptr) {
  sub_1d_type npages;
  sub_1d_type nwritten=0;
  sub_1d_type sub;

  //if the current swap file stream is null, then set it to fileptr
  //and advance the byte location forward the size of the dataset:
  if (swap == NULL) {
    swap=fileptr;
    swap_start=ftell(swap);
    fseek(swap, this->n_data*sizeof(dtype), SEEK_CUR);
    return this->n_data*sizeof(dtype);
  }

  //write out all the data:
  npages=this->n_data/page_size;
  for (sub_1d_type i=0; i<npages; i++) {
    sub=i*page_size;
    check_page(sub);
    nwritten+=fwrite(this->data, sizeof(dtype), page_size, fileptr);
  }
  sub=npages*page_size;
  if (check_page(sub) != -1) {
    nwritten+=fwrite(this->data, sizeof(dtype), this->n_data-sub, fileptr);
  }

  return nwritten;

}

template <class dtype>
void dependent_swap<dtype>::print_meta() {
  printf("Type: %d\n", this->type);
  for (rank_type i=0; i<this->rank-1; i++) printf("%d X ", this->dim[i]);
  printf("%d = %d\n", this->dim[this->rank-1], this->n_data);
  for (rank_type i=0; i<this->rank; i++) {
    this->dependencies[i]->print();
  }
}

template <class dtype>
errtype dependent_swap<dtype>::interpol(dtype &value, interpol_index *indices) {
  ind_type lind[this->rank], hind[this->rank];
  ind_type l, h;
  sub_1d_type sub=0, mult;
  long n=1 << this->rank;
  long order;
  interpol_index frac;
  long norm;
  double weight;
  dtype data_sub;
  errtype result=0;

  if (this->n_data == 0) return NO_DATA;

  //determine the upper and lower values for the indices:
  norm=1;
  for (rank_type i=0; i<this->rank; i++) {
    if (this->dim[i]==1) {
//      l=0;
//      h=0;
      norm*=2;		//also determine the normalization coeff.
    } else {
      l=(ind_type) indices[i];
      //check for out-of-range indices (do extrapolation):
      if (l>=this->dim[i]-1) {
        l=this->dim[i]-2; 
	result=INDEX_OUT_OF_RANGE[i];
      } else if (l<0) {
	l=0;
	result=INDEX_OUT_OF_RANGE[i];
      }
      lind[i]=l;
      hind[i]=l+1;
    }
  }

  //calculate the weights and apply them:
  value=0;
  for (long i=0; i<n; i++) {
    weight=1;
    sub=0;
    mult=1;
    for (rank_type j=0; j<this->rank; j++) {
      if (this->dim[j] == 1) continue;
      order= 1 << j;
      if ((i & order) == 0) {
        frac=(interpol_index) hind[j]-indices[j];
        sub=sub+mult*lind[j];
      } else {
        frac=indices[j]-(interpol_index) lind[j];
        sub=sub+mult*hind[j];
      }
      weight*=frac;
      mult=mult*this->dim[j];
   }
   get_1d(data_sub, sub);
   value+=(dtype) weight*data_sub;
  }
  value/=(dtype) norm;

//  printf("Normalization coef: %f\n", norm);

  return result;
}

template <class dtype>
errtype dependent_swap<dtype>::insert(rank_type r, ind_type index, 
		ind_type ni) {
	sub_1d_type j, k;
	sub_1d_type n;
	sub_1d_type stride, offset, gap;
	dtype *new_data;
	sub_1d_type d;

	//this will be difficult, therefore we put it off as long as possible...

	return -1;

}

template <class dtype>
errtype dependent_swap<dtype>::del(rank_type r, ind_type index, 
		ind_type ni) {
	long i, j;
	long nold, n;
	long stride, offset, gap;
	dtype *new_data;

	return -1;

}

template <class dtype>
errtype dependent_swap<dtype>::preload(dtype *new_data, sub_1d_type n, 
		sub_1d_type offset) {
  if (offset >= this->n_data) return 0;

  fseek(swap, swap_start+offset*sizeof(dtype), SEEK_SET);
  if (offset+n > this->n_data) n=this->n_data-offset;
  fchange=1;
  dchange=0;
  return fwrite(new_data, sizeof(dtype), n, swap);
}

/*
template <class dtype>
long dependent_swap<dtype>::preload(dependent<dtype> *new_data,  
		sub_1d_type offset) {
  sub_1d_type n;
  if (offset >= this->n_data) return 0;
  n=new_data->n_data;

  fseek(swap, swap_start+offset*sizeof(dtype), SEEK_SET);
  if (offset+n > this->n_data) n=this->n_data-offset;
  fchange=1;
  dchange=0;
  return fwrite(new_data->data, sizeof(dtype), n, swap);
}
*/

template <class dtype>
errtype dependent_swap<dtype>::read_chunk(FILE *fileptr, sub_1d_type offset, 
		sub_1d_type n) {
  long nread;
  dtype *new_data;

//  printf("Read chunk called for swapping dependent\n");

//  printf("this->n_data: %d, Offset: %d, n read: %d\n", this->n_data, offset, n);

  if (offset >= this->n_data || offset < 0) return 0;

  fseek(swap, swap_start+offset*sizeof(dtype), SEEK_SET);
  if (offset+n > this->n_data) n=this->n_data-offset;

  new_data=new dtype[n];
  nread=fread(new_data, sizeof(dtype), n, fileptr);

  fwrite(new_data, n*sizeof(dtype), 1, swap);
  fchange=1;
  dchange=0;

  delete [] new_data;

  return nread;
}

template<>
void dependent_swap<float>::set_type() {
  this->missing=FLOAT_MISSING;
  type=DEPENDENT_FLOAT_S;
}

template<>
void dependent_swap<double>::set_type() {
  this->missing=DOUBLE_MISSING;
  type=DEPENDENT_DOUBLE_S;
}

template<>
void dependent_swap<int32_t>::set_type() {
  this->missing=INT32_MISSING;
  type=DEPENDENT_LONG_S;
}

template class dependent_swap<float>;
template class dependent_swap<double>;
template class dependent_swap<int32_t>;

} //end namespace datasets
} //end namespace libpetey

