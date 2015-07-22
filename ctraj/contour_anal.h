#ifndef ctraj_CONTOUR_ANAL__H
#define ctraj_CONTOUR_ANAL__H

#include <stdint.h>
#include "time_class.h"

namespace ctraj {

  using namespace libpetey;

  int32_t read_bev_file(char *filename, int32_t index, time_class &t, float *&lon, float *&lat);

  int32_t bev_index(char *filename, time_class *&t, int32_t *&n);

  long read_ascii_contour(FILE *fs, float *&x, float *&y);

  //returns the uncertainty fraction of a boundary, given an error value:
  float uncertainty_fraction(float *lon, float *lat, long npt, float epsilon, 
		int hemi, int min_uncertain, int max_total, int wrap_flag=0);

  long count_boxes_global(float *lon, float *lat, long npt, float epsilon, 
		int hemi, int wrapflag=0);

  long measure_contour(float *lon, float *lat, long npt, float epsilon, int wrap_flag=0);

}

#endif

