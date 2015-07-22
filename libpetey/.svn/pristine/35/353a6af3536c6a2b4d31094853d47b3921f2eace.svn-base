class iter2_obj {
  ind_type *ind;
  rank_type rank;

  ind_type *map1;
  ind_type *map2;

  ind_type *dim;
  rank_type n;

  simple_dataset *dep;

  dependent_dataset *ds1;
  dependent_dataset *ds2;

  iter2_obj(dependent<dtype1> *d1, dependent<dtype2> *d2);

  ~iter2_obj;

  void zero();

  errtype next(sub_1d_type *i1, sub_1d_type *i2);

  template <class dtype1, class dtype2, class dtype3>
  dependent<dtype3> * iter2_obj::bi_op(dtype3 &bi_op(dtype1, dtype2));
};

//there's a lot of extra, unnecessary redundancies here
//but they're probably necessary for efficiency reasons...
iter2_obj::iter2_obj(dependent_dataset *d1, dependent_dataset *d2<dtype2>) {

  long *sind1;
  long *sind2;

  long k1, k2;

  ds1=d1;
  ds2=d2;

  sind1=heapsort(ds1->dim, ds1->rank);
  sind2=heapsort(ds1->dim, ds2->rank);

  map1=new ind_type[ds1->rank+ds2->rank];
  map2=new ind_type[ds1->rank+ds2->rank];

  dim=new rank_type[ds1->rank+ds2->rank];

  dep=new simple_dataset[ds1->rank+ds2->rank];

  k1=0;
  k2=0;
  n=0;
  while (k1 < ds1->rank && k2 < ds2->rank) {
    if (s1->dependencies[sind1[k1]] == s2->dependencies[sind2[k2]]) {
      map1[n]=sind[k1];
      map2[n]=sind[k2];
      dim[n]=s1->dim[sind1[k1]];
      dep[n]=s1->dependencies[sind1[k1]];
      k1++;
      k2++;
      n++;
    } else if (s1->dependencies[sind1[k1]] < s2->dependencies[sind2[k2]]) {
      map1[n]=sind[k1];
      dim[n]=s1->dim[sind1[k1]];
      dep[n]=s1->dependencies[sind1[k1]];
      k1++;
      n++;
    } else {
      map2[n]=sind[k2];
      dim[n]=s2->dim[sind2[k2]];
      dep[n]=s2->dependencies[sind2[k2]];
      k2++;
      n++;
    }
  }

  ind=new ind_type[n];

  delete [] sind1;
  delete [] sind2;

}


iter2_obj::~iter2_obj {
  delete [] map1;
  delete [] map2;
  delete [] ind;
  delete [] dim;
}

void iter2_obj::zero() {
  for (ind_type i=0; i<=n; i++) {
    ind[i]=0;
  }
  rank=0;
}

errtype iter2_obj::next(sub_1d_type *i1, sub_1d_type *i2);
  ind_type *ind1;
  ind_type *ind2;

  ind1=new ind_type[ds1->rank];
  ind2=new ind_type[ds2->rank];

  i1=ds1->calc_sub(ind1);
  i2=ds2->calc_sub(ind2);

  errtype err=0;

  ind[rank]++;

  if (ind[rank] >= dim[rank]) {
    ind[rank]=0;
    rank++;
    if (rank >= ndim) rank=0;
    err=1;
  }

  for (ind_type i=0; i<ndim; i++) {
    ind1[map1[i]]=ind[i];
    ind2[map2[i]]=ind[i];
  }

  delete [] ind1;
  delete [] ind2;

}

template <class dtype1, class dtype2, class dtype3>
dependent<dtype3> * iter2_obj::bi_op(dtype3 &bi_op(dtype1, dtype2)) {

  dependent<dtype1> *d1;
  dependent<dtype2> *d2;
  dependent<dtype3> *result;

  sub_1d_type sub1;
  sub_1d_type sub2;

  dtype1 val1;
  dtype2 val2;

  d1=(dependent<dtype1> *) ds1;
  d2=(dependent<dtype2> *) ds2;

  result=new dependent<dtype3>(dep, n);

  zero();

  while (next(&sub1, &sub2) == 0) {
    d1->get_1d(val1, sub1);
    d2->get_1d(val2, sub2);
    result->cel(bi_op(val1, val2), ind);
  }

  return result;

}


