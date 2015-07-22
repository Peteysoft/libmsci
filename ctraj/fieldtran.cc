//lon-lat ascii to azimuthal-equidistant binary...
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <getopt.h>

#include "linked.h"

#include "simple_temp.h"
#include "dependent_temp.h"

#include "coordtran.h"

void ll2ae(dependent<float> *q1, 
		simple<float> *lon, 
		simple<float> *lat,
		dependent<float> *q2,
		simple<float> *x,
		simple<float> *y,
		int hemi) {

  dependent<interpol_index> **k;

  k=new dependent<interpol_index> *[2];

  k[0]=new dependent<interpol_index>(x, y);
  k[1]=new dependent<interpol_index>(x, y);

  intcoeff2(lon, lat, x, y, k[0], k[1]);

  q1->interpol(k, q2);
}

//inverse of previous (should be more closely connected...)
void ae2ll(dependent<float> *q1, 
		simple<float> *lon, 
		simple<float> *lat,
		dependent<float> *q2,
		simple<float> *x,
		simple<float> *y,
		int hemi) {

  dependent<interpol_index> **k;

  k=new dependent<interpol_index> *[2];

  k[0]=new dependent<interpol_index>(x, y);
  k[1]=new dependent<interpol_index>(x, y);

  intcoeff(lon, lat, x, y, k[0], k[1]);

  q1->interpol(k, q2);
}

//another inverse (this one is exact...)
void apply_map(float *vec, 
		dependent<sub_1d_type> *map) {
  ind_type nx, ny;

  map->get_dim(nx, ny);

  

