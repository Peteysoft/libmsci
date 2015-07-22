#ifndef AGF_METRIC2_H_INCLUDED
#define AGF_METRIC2_H_INCLUDED 1

#include "agf_defs.h"

namespace libagf {

  //constants and global variables:
  extern const enum_a global_nmetric;		//number of metrics
  //default metric:
  extern real_a (* global_metric2) (real_a *, real_a *, dim_ta m);
  //list of metrics:
  extern real_a (* global_metric2_pointer_list[]) (real_a *, real_a *, dim_ta m);
  //names and numbers of metrics (for help screens):
  extern char global_metric_type[];

  //all metrics are squared unless otherwise noted:

  //Cartesian metric (squared):
  template <class real>
  real metric2(real *v1, 		//first coordinate
		real *v2, 		//second coordinate
		dim_ta m);		//number of dimensions

  //Manhattan metric:
  template <class real>
  real manhattan2(real *v1, real *v2, dim_ta m);

  //quartic:
  template <class real>
  real quad_metric2(real *v1, real *v2, dim_ta m);

  //simple maximum:
  template <class real>
  real max_metric2(real *v1, real *v2, dim_ta m);

  //quartic to fourth power:
  template <class real>
  real quad_metric4(real *v1, real *v2, dim_ta m);

  //Manhattan (not squared):
  template <class real>
  real manhattan(real *v1, real *v2, dim_ta m);
}

#endif

