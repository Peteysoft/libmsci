#ifndef ctraj_COORDTRAN_H
#define ctraj_COORDTRAN_H

#include "simple_temp.h"
#include "dependent_temp.h"

#include "az_eq_t.h"

namespace ctraj {
  using namespace libpetey::datasets;

  simple<float> * create_grid(long ngrid, float slen);

  //given a pair of grids, calculate interpolation coefficients
  //to go from lon-lat to azimuthal equidistant...
  void intcoeff(simple<float> *lon,
		simple<float> *lat,
		simple<float> *xgrid,
		simple<float> *ygrid,
		dependent<interpol_index> *c1, 
		dependent<interpol_index> *c2,
		short hemi);

  //opposite direction:
  void intcoeff2(simple<float> *lon,
		simple<float> *lat,
		simple<float> *xgrid,
		simple<float> *ygrid,
		dependent<interpol_index> *c1, 
		dependent<interpol_index> *c2,
		short hemi);

  //interpolate a field:
  void intq(dependent<float> q0,
		dependent<interpol_index> c1,
		dependent<interpol_index> c2,
		dependent<float> q1);

  //interpolate a velocity field (includes metric transformations):


  //create mapping for az. eq. gridding with corners to one without:
  void nocorner_map(simple<float> *xgrid, simple<float> *ygrid, dependent<sub_1d_type> *map, sub_1d_type &nmap, az_eq_t<float> *m=NULL);

  //bloody stupid...
  sub_1d_type calc_nmap(ind_type ngrid);

} //end namespace ctraj

#endif


