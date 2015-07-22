
//*FLAG*

// *** BE CAREFULL! ***

// THIS FILE IS AN ATTEMPTED UPDATE TO tracer_anal.h

// DO NOT CONFUSE THE TWO

// CURRENTLY NOT USED ANYWHERE...

//*FLAG*


#include <stdint.h>

#include "sparse.h"
#include "meas_data.h"

#include "ctraj_vfield_base.h"
#include "ctraj_tfield_base.h"

float ** tracer_multiply(sparse_matrix *a, int32_t N, float *v0);

float ** tracer_multiply_renorm(sparse_matrix *a, int32_t N, float *v0);

//assumes standard side-length of 20000 km
void tracer_interp(float **qall, time_class *t, int32_t nt, int32_t np, meas_data *data, long nsamp);

float ** pc_proxy(sparse_matrix *matall, time_class *t, int32_t N, int32_t nall, meas_data *samp, long nsamp, int32_t nev, int32_t ncv, int cflag=0);

//returns the coefficients:
float * proxy_tracer(float **tall, int32_t np, time_class *t, int32_t nall, meas_data *samp, long nsamp, int32_t order);

//new paradigm:
template <class real>
void tracer_map(ctraj_tfield_base<real> * tracer, 
		ctraj_vfield_base<real> *vfield,
		sparse<int32_t, real> *map,
		real t, real dt, int32_t nt);

