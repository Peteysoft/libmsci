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

  //if flag is set, .qerr is printed as regular data with a wider format
  int write_meas(meas_data *data, long n, FILE *fs, int flag=0);

  inline meas_data *copy_meas(meas_data *dat, long n) {
    meas_data *result=new meas_data[n];
    for (long i=0; i<n; i++) result[i]=dat[i];
    return result;
  }

  meas_data * sort_meas(meas_data *dat, long n);

  meas_data * randomize_meas(meas_data *dat, long n);

  meas_data *select_meas(time_class t0, time_class t1, meas_data *dat, long n1, long *n2, int hemi=0);

  meas_data *select_lat_range(meas_data *dat, long n1, float min, float max, long *n2);

  double correlate_meas(meas_data *samp1, 	//first set of samples
		  meas_data *samp2, 		//second set of samples
		  long n,			//number of samples
		  int flag=0);			//correlate error
} //end namespace ctraj

#endif

