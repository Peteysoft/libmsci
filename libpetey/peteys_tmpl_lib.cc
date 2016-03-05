//
// Copywrite 2004 Peter Mills.  All rights reserved.
//
// Compiles a number of instantiations of templated subroutines into
// a single file
//
#include <stdint.h>
#include <vector>

#include "time_class.h"
#include "string_petey.h"

#include "bin_search.cc"
#include "heapsort_tmpl.cc"
#include "rk_dumb_tmpl.cc"
#include "treesort.cpp"

using namespace std;

namespace libpetey {

//binary search:
template long bin_search_g<float>(float *, long, float, long);
template long bin_search_g<int32_t>(int32_t *, long, int32_t, long);
template long bin_search_g<int64_t>(int64_t *, long, int64_t, long);
template long bin_search_g<double>(double *, long, double, long);
template long bin_search_g<time_class>(time_class *, long, time_class, long);
template long bin_search_g<string_petey>(string_petey *, long, string_petey, long);
template long bin_search<vector<int> >(vector<int> *, long, vector<int>,
long);
//template long bin_search<int>(int *, long, int, long);

template double interpolate<float>(float *, long, float, long);
//template double interpolate<long>(long *, long, long, long);
template double interpolate<int32_t>(int32_t *, long, int32_t, long);
template double interpolate<int64_t>(int64_t *, long, int64_t, long);
template double interpolate<double>(double *, long, double, long);
template double interpolate<time_class>(time_class *, long, time_class, long);
template double interpolate<string_petey>(string_petey *, long, string_petey,
long);

//sorting etc.:
template void heapsort_inplace<float>(float *, long);
template void heapsort_inplace<int32_t>(int32_t *, long);
template void heapsort_inplace<int64_t>(int64_t *, long);
template void heapsort_inplace<double>(double *, long);
template void heapsort_inplace<time_class>(time_class *, long);
template void heapsort_inplace<string_petey>(string_petey *, long);

template long * heapsort<float>(float *, long);
template long * heapsort<int32_t>(int32_t *, long);
template long * heapsort<int64_t>(int64_t *, long);
template long * heapsort<double>(double *, long);

template long * heapsort<time_class>(time_class *, long);
template long * heapsort<string_petey>(string_petey *, long);

template long * heapsort<vector<int> >(vector<int> *, long);
template long * heapsort<vector<float> >(vector<float> *, long);
template long * heapsort<vector<double> >(vector<double> *, long);

template long * heapsort<void *>(void **, long);

template void heapsort<float>(float *data, long *, long);
template void heapsort<int32_t>(int32_t *data, long *, long);
template void heapsort<int64_t>(int64_t *data, long *, long);
template void heapsort<double>(double *data, long *, long);

template void heapsort<time_class>(time_class *data, long *, long);
template void heapsort<string_petey>(string_petey *data, long *, long);

template void heapsort<vector<int> >(vector<int> *, long *, long);
template void heapsort<vector<float> >(vector<float> *, long *, long);
template void heapsort<vector<double> >(vector<double> *, long *, long);

template void heapsort<void *>(void **data, long *ind, long n);

template long * treesort<float>(float *data, long n);
template long * treesort<int32_t>(int32_t *data, long n);
template long * treesort<int64_t>(int64_t *data, long n);
template long * treesort<double>(double *data, long n);
template long * treesort<time_class>(time_class *data, long n);
template long * treesort<string_petey>(string_petey *data, long n);

template void reverse<float>(float *, long);
template void reverse<int32_t>(int32_t *, long);
template void reverse<int64_t>(int64_t *, long);
template void reverse<double>(double *, long);
template void reverse<time_class>(time_class *, long);
template void reverse<string_petey>(string_petey *, long);

template float * map_vector<float>(float *, long *, long);
template int32_t * map_vector<int32_t>(int32_t *, long *, long);
template int64_t * map_vector<int64_t>(int64_t *, long *, long);
template double * map_vector<double>(double *, long *, long);
template time_class * map_vector<time_class>(time_class *, long *, long);
template string_petey * map_vector<string_petey>(string_petey *, long *, long);
template char ** map_vector<char *>(char **, long *, long);

template void map_vector_inplace<float>(float *, long *, long);
template void map_vector_inplace<double>(double *, long *, long);
template void map_vector_inplace<int32_t>(int32_t *, long *, long);
template void map_vector_inplace<int64_t>(int64_t *, long *, long);
template void map_vector_inplace<time_class>(time_class *, long *, long);
template void map_vector_inplace<string_petey>(string_petey *, long *, long);

//ODEs:
template float **rk_dumb<float, float>(float, float *, long, float, long,
		void (* derivs) (float, float *, float *, long));
template double **rk_dumb<float, double>(float, double *, long, float, long,
		void (* derivs) (float, double *, double *, long));
template float **rk_dumb<double, float>(double, float *, long, double, long,
		void (* derivs) (double, float *, float *, long));
template double **rk_dumb<double, double>(double, double *, long, double, long,
		void (* derivs) (double, double *, double *, long));
template float **rk_dumb<time_class, float>(time_class, float *, long, time_class, long,
		void (* derivs) (time_class, float *, float *, long));
template double **rk_dumb<time_class, double>(time_class, double *, long, time_class, long,
		void (* derivs) (time_class, double *, double *, long));

template void rk_dumb<float, float>(float, float *, long, float, long,
		float **,
		void (* derivs) (float, float *, float *, long));
template void rk_dumb<float, double>(float, double *, long, float, long,
		double **,
		void (* derivs) (float, double *, double *, long));
template void rk_dumb<double, float>(double, float *, long, double, long,
		float **,
		void (* derivs) (double, float *, float *, long));
template void rk_dumb<double, double>(double, double *, long, double, long,
		double **,
		void (* derivs) (double, double *, double *, long));
template void rk_dumb<time_class, float>(time_class, float *, long, 
		time_class, long, float **,
		void (* derivs) (time_class, float *, float *, long));
template void rk_dumb<time_class, double>(time_class, double *, long, 
		time_class, long, double **,
		void (* derivs) (time_class, double *, double *, long));

} //end namespace libpetey

