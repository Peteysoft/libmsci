#include <stdlib.h>
//#include <stdio.h>
#include <iostream>
#include <math.h>
#include <stdint.h>

#include "dependent_temp.h"

namespace libpetey {
namespace datasets {

template <class dtype>
dependent<dtype>::dependent() {
  set_type();
  data=NULL;
  rank=0;
  dim=NULL;
  dependencies=NULL;
}

template <class dtype>
dependent<dtype>::~dependent() {
  delete [] data;
}

template <class dtype>
dependent<dtype>::dependent(stype **s, rank_type ndep) {
	init(s, ndep);
}
	
template <class dtype>
errtype dependent<dtype>::init(stype **s, rank_type ndep) {
  rank_type i;

  set_type();

  rank=ndep;
  dim=new ind_type[rank];
  dependencies=new stype *[rank];

  n_data=1;
  for (i=0; i<rank; i++) {
    dependencies[i]=s[i];
    s[i]->add_dependent(this, i);
    dim[i]=s[i]->nel();
    n_data=n_data*dim[i];
  }

  if (n_data==0) return 0;
  data=new dtype[n_data];
  for (sub_1d_type i=0; i<n_data; i++) data[i]=missing;
  //printf("%d\n", n_data);

  return n_data;
}

template <class dtype>
dependent<dtype>::dependent(stype *s1) {

	init(&s1, 1);
}


template <class dtype>
dependent<dtype>::dependent(stype *s1, stype *s2) {
	stype *dep[2]={s1, s2};
	init(dep, 2);
}

template <class dtype>
dependent<dtype>::dependent(stype *s1, stype *s2, stype *s3) {
	stype *dep[3]={s1, s2, s3};
	init(dep, 3);
}

template <class dtype>
dependent<dtype>::dependent(stype *s1, stype *s2, stype *s3, stype *s4) {
	stype *dep[4]={s1, s2, s3, s4};
	init(dep, 4);
}

template <class dtype>
long dependent<dtype>::read(FILE *fileptr) {
  //read in the data:
  return fread(data, sizeof(dtype)*n_data, 1, fileptr);

}

template <class dtype>
long dependent<dtype>::write(FILE *fileptr) {
  //write out all the data:
  return fwrite(data, sizeof(dtype)*n_data, 1, fileptr);
}


template <class dtype>
errtype dependent<dtype>::interpol(dtype &value, interpol_index *indices) {
  ind_type lind[rank], hind[rank];
  ind_type l, h;
  sub_1d_type sub, mult;
  long n=1 << rank;
  long order;
  interpol_index frac;
  long norm;
  interpol_index weight;
  errtype result=0;

  if (n_data == 0) return NO_DATA;

  //determine the upper and lower values for the indices:
  norm=1;
  for (rank_type i=0; i<rank; i++) {
    if (dim[i]==1) {
//      l=0;
//      h=0;
      norm*=2;		//also determine the normalization coeff.
    } else {
      l=(ind_type) indices[i];
      if (l>=dim[i]-1) {
        l=dim[i]-2;
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
    for (rank_type j=0; j<rank; j++) {
      if (dim[j] == 1) continue;
      order= 1 << j;
      if ((i & order) == 0) {
        frac=(interpol_index) hind[j]-indices[j];
        sub=sub+mult*lind[j];
      } else {
        frac=indices[j]-(interpol_index) lind[j];
        sub=sub+mult*hind[j];
      }
      weight*=frac;
      mult=mult*dim[j];
    }
    value+=(dtype) weight*data[sub];
  }
  value/=(dtype) norm;

//  printf("Normalization coef: %f\n", norm);

  return result;
}

template <class dtype>
errtype dependent<dtype>::insert(rank_type r, ind_type index, ind_type ni) {
	sub_1d_type j, k;
	sub_1d_type n;
	sub_1d_type stride, offset, gap;
	dtype *new_data;
	ind_type d;

	//cout << "Attempt to insert hyperslab at dim " << r << ", index " <<
	//		index << "\n";
	fprintf(datasets_log, "Attempt to insert hyperslab at dim %d, index, %d\n", r,	index);

	if (r >= rank || r < 0) return RANK_ERROR;
	if (index < 0 || index > dim[r]) return INDEX_OUT_OF_RANGE[r];

	//cout << "Current dimensions: ";
	fprintf(datasets_log, "Current dimensions: ");
	for (rank_type i=0; i<rank;i++) fprintf(datasets_log, "%d ", dim[i]);
	//cout << ":" << n_data << "\n";
	fprintf(datasets_log, ": %d \n", n_data);

	d=dim[r];
	dim[r]+=ni;

        if (d != 0) {
	  n=n_data/d*(d+ni);
	} else {
	  n=1;
	  for (rank_type i=0;i<rank;i++) n=n*dim[i];
	}

	new_data=new dtype[n];

	gap=ni;
	for (rank_type i=0; i<r; i++) {
		gap=gap*dim[i];
	}
	stride=gap*d;
//	stride=gap*dim[r];
	offset=gap*index;

	//cout << "Offset: " << offset <<"; stride: " << stride <<
	//		"; gap: " << gap << "\n";
	fprintf(datasets_log, "Offset: %d; stride: %d; gap %d\n", offset, stride, gap);

	j=0;
	for (sub_1d_type i=0;i<n_data;i++) {
		//printf("%4d%4d\n", i, j);
		new_data[j]=data[i];
		if ((i-offset+1) % stride == 0) {
			k=j+gap;
			for (j=j+1; j<=k; j++) new_data[j]=missing;
		} else {
			j++;
		}
	}

	//for (i=0; i< n ; i++) cout << new_data[i] << " ";
	//cout << "\n";

	n_data=n;
	if (data != NULL) delete [] data;
	data=new_data;

	return NO_PROBLEM;

}

template <class dtype>
errtype dependent<dtype>::del(rank_type r, ind_type index, ind_type ni) {
	sub_1d_type j;
	sub_1d_type nold, n;
	sub_1d_type stride, offset, gap;
	dtype *new_data;

	if (r >= rank || r < 0) return RANK_ERROR;
	if (index < 0 || index >= dim[r]) return INDEX_OUT_OF_RANGE[r];
	if (ni > dim[r]) return MISC_INDEXING_ERROR;

        n=n_data/dim[r]*(dim[r]-ni);

	if (n == 0) {
	  delete [] data;
	  data=NULL;
	  dim[r]--;
	  n_data=n;
	  return 0;
	}

	new_data=new dtype[n];

	gap=ni;
	for (rank_type i=0; i<r; i++) {
		gap=gap*dim[i];
	}
	stride=gap*dim[r];
	offset=gap*index-1;

	j=0;
	for (sub_1d_type i=0;i<n;i++) {
		new_data[i]=data[j];
		if ((i-offset+1) % stride == 0) {
			j=j+gap;
		} else {
			j++;
		}
	}

	dim[r]=dim[r]-ni;
	n_data=n;
	delete [] data;
	data=new_data;

	return NO_PROBLEM;

}

template <class dtype>
errtype dependent<dtype>::preload(dtype *new_data, sub_1d_type n, 
			sub_1d_type offset) {
  if (offset+n > n_data) n=n_data-offset;
  for (sub_1d_type i=0; i<n; i++) {
    data[i+offset]=new_data[i];
  }
  return n;
}

template <class dtype>
errtype dependent<dtype>::read_chunk(FILE *fileptr, sub_1d_type offset, 
			sub_1d_type n) {
  long nread;

//  printf("Read chunk called for regular dependent\n");

  if (offset > n_data) return SUBSCRIPT_OUT_OF_RANGE;
  if (offset + n > n_data) n=n_data-offset;
  nread=fread(&data[offset], sizeof(dtype), n, fileptr);

  return nread;
}

template <class dtype>
errtype dependent<dtype>::iter_get_next(ind_type *ind_cur, rank_type &rank_cur) {
  for (rank_type i=0; i<rank; i++) ind_cur[i]=0;
  rank_cur=0;
}

template <class dtype>
void dependent<dtype>::print_meta() {
  printf("Type: %d\n", type);
  for (rank_type i=0; i<rank-1; i++) printf("%d X ", dim[i]);
  printf("%d = %d\n", dim[rank-1], n_data);
  printf("\n");
  for (rank_type i=0; i<rank; i++) {
    dependencies[i]->print();
  }
}


/*
template <class dtype>
dependent<interpol_index> *dependent<dtype>::intcoeff(stype *s) {
  dependent<interpol_index> *c;
  c=new dependent<interpol_index>(dependencies, rank);
  for (sub_1d_type i=0; i<n_data; i++) c->cel_1d(s->interpol(data[i]), i);
  return c;
}

template <class dtype>
dependent<dtype> *dependet<dtype>::interpolate(dependent<interpol_index> **c) {
  interpol_index ind[rank];
  dependent<dtype> *result;
  result=new dependent<dtype>(c[0]->dependencies, c[0]->rank);

  for (sub_1d_type i=0; i<result->n_data; i++) {
    for (rank_type k=0; k<rank; k++) ind[k]=c[k]->data[i];
    result->data[i]=interpol(ind);
  }

  return result;
}

template <class dtype>
errtype dependent<dtype>::interpolate(dependent<interpol_index> **c, dependent<dtype> *result) {
  interpol_index ind[rank];

  for (sub_1d_type i=0; i<result->n_data; i++) {
    for (rank_type k=0; k<rank; k++) c[k]->get_1d(ind[k], i);
    result->data[i]=interpol(ind);
  }

  return result;
}
*/

template<>
void dependent<float>::set_type() {
  missing=FLOAT_MISSING;
  type=DEPENDENT_FLOAT;
}

template<>
void dependent<double>::set_type() {
  missing=DOUBLE_MISSING;
  type=DEPENDENT_DOUBLE;
}

template<>
void dependent<int32_t>::set_type() {
  missing=INT32_MISSING;
  type=DEPENDENT_INT32;
}

template<>
void dependent<int64_t>::set_type() {
  missing=INT64_MISSING;
  type=DEPENDENT_INT64;
}


template class dependent<float>;
template class dependent<double>;
template class dependent<int32_t>;
template class dependent<int64_t>;

} //end namespace datasets
} //end namespace libpetey
