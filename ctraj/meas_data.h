#ifndef ctraj_MEAS_DATA_H
#define ctraj_MEAS_DATA_H

#include <stdio.h>

#include "time_class.h"

#include "ctraj_defaults.h"

namespace ctraj {
  using namespace libpetey;

  struct meas_data {
	time_class t;
	float lon;
	float lat;
	double q;
	double qerr;
  };

  meas_data *generate_random_global(time_class t0, time_class tf, long n, int hemi=0);

  meas_data *read_meas(FILE *fs, long *n, int32_t twid=TFIELD_WIDTH);

  inline meas_data *read_meas(char *infile, long *n, int32_t twid=TFIELD_WIDTH) {
    FILE *fs;
    meas_data *data;

    fs=fopen(infile, "r");
    data=read_meas(fs, n, twid);
    fprintf(stderr, "%ld measurements read from file, %s\n", *n, infile);

    return data;

  }

  int write_meas(meas_data *data, long n, FILE *fs);

  meas_data * sort_meas(meas_data *dat, long n);

  meas_data * randomize_meas(meas_data *dat, long n);

  meas_data *select_meas(time_class t0, time_class t1, meas_data *dat, long n1, long *n2, int hemi=0);

} //end namespace ctraj

#endif

