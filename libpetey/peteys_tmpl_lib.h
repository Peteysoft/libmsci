//
// Copyright 2004 Peter Mills.  All rights reserved.
//
// This is the header file for "peteys_tmpl_lib.cc."  Only templates are
// given here.  The desired instantiation may or may not exist--check the
// main source file.

#ifndef TMPL_LIB_INCLUDED
#define TMPL_LIB_INCLUDED

#pragma interface

namespace libpetey {

  template <class dt>
  void heapsort_inplace(dt *data, long n);

  template <class dt>
  long * heapsort(dt *data, long n);

  template <class dt>
  void heapsort(dt *data, long *ind, long n);

  template <class dt>
  long * treesort(dt *data, long n);

  //because we're too lazy to write a heapsort routine that sorts in descending order:
  template <class dt>
  void reverse(dt *data, long n);

  //maps a vector to a set of new positions given by a vector of indices:
  template <class dt>
  dt * map_vector(dt * vector, long * indices, long n);

  //maps a vector to a set of new positions given by a vector of indices:
  template <class dt>
  void map_vector_inplace(dt * vector, long * indices, long n);

  //uses a bisection to search an ordered list:
  //in ascending order:
  template <class dt>
  long bin_search (dt *list, long n, dt value, long last_ind=-1);
  template <class dt>
  //in descending order:
  long bin_search_r (dt *list, long n, dt value, long last_ind=-1);
  //both ascending and descending:
  template <class dt>
  long bin_search_g (dt *list, long n, dt value, long last_ind=-1);
  //performs an interpolation at the same time:
  template <class dt>
  double interpolate (dt *list, long n, dt value, long last_ind=-1);

  //templated Runge-Kutta integrator:
  template <class ind_type, class dep_type>
  dep_type **rk_dumb(ind_type t0, dep_type *xinit, long nx, ind_type dt, long nt,
                 void (* derivs) (ind_type, dep_type *, dep_type *, long) );

  template <class ind_type, class dep_type>
  void rk_dumb(ind_type t0, dep_type *xinit, long nx, ind_type dt, long nt,
		dep_type ** xvec,
		 void (* derivs) (ind_type, dep_type *, dep_type *, long) );

  //reverse a vector:
  template <class dt>
  void reverse(dt *data, long n);

  //maps a vector to a set of new positions given by a vector of indices:
  template <class dt>
  inline void map_vector(dt * vector, long * indices, dt *result, long n) {
    long i;

    for (i=0;i<n;i++) result[i]=vector[indices[i]];
  }

  //invert a mapping
  inline long *invert_mapping(long *ind, long n) {
    long *result;
    result=new long[n];
    for (long i=0; i<n; i++) result[ind[i]]=i;
    return result;
  }

}

#endif
