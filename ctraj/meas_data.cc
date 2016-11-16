#include <stdio.h>
#include <sys/timeb.h>
#include <math.h>
#include <string.h>

#include "randomize.h"

#include "time_class.h"
#include "read_ascii_all.h"

#include "linked.h"
#include "peteys_tmpl_lib.h"

#include "meas_data.h"

using namespace libpetey;

namespace ctraj {

meas_data *generate_random_global(time_class t0, time_class tf, long n, int hemi) {

  meas_data *result;
  time_class *t, tdiff;
  double sinlat;
  double td, dt;

  result=new meas_data[n];

  tdiff=tf-t0;
  dt=tdiff;
  t=new time_class[n];

  for (long i=0; i<n; i++) {
    td=dt*ranu();
    t[i]=t0;
    t[i].add(td);
    //t[i].write_string(tstring);
    //printf("%23s\n", tstring);
  }
  //sort the dates::
  heapsort_inplace(t, n);

  for (long i=0; i<n; i++) {
    result[i].t=t[i];
    result[i].q=0;
    result[i].qerr=0;
    result[i].lon=ranu()*360;
    if (hemi!=0) {
      sinlat=hemi*ranu();
    } else {
      sinlat=2.*ranu()-1.;
    }
    result[i].lat=asin(sinlat)*180/M_PI;
  }
  delete [] t;

  return result;
}

//we specify the width of the date field explicitly
//=ugly hack...
meas_data *read_meas(FILE *fs, long *n, int32_t twid) {
  linked_list<char *> lines;
  char *line;
  char **data;
  meas_data *data2;

  char tstring[30];

  while (feof(fs) == 0) {
    if ((line=fget_line(fs))==NULL) {
      delete [] line;
      break;
    }
    lines.add(line);
  }

  data=lines.make_array(*n);
  data2=new meas_data[*n];

  for (long i=0; i<*n; i++) {
    data2[i].qerr=0.;
    sscanf(data[i]+twid, "%f %f %lf %lf", &data2[i].lon, &data2[i].lat, &data2[i].q, &data2[i].qerr);
    strncpy(tstring, data[i], twid+1);
    tstring[twid]='\0';
    data2[i].t.read_string(tstring);
    //printf("%s %f %f %lf %lf\n", tstring, data2[i].lon, data2[i].lat);
  }

  for (long i=0; i<*n; i++) delete [] data[i];
  delete [] data;

  return data2;

}

int write_meas(meas_data *data, long n, FILE *fs) {
  char tstring[30];

  for (long i=0; i<n; i++) {
    data[i].t.write_string(tstring);
    //this is annoying: you can't stick a macro in a string (I don't think...)
    fprintf(fs, "%23s %9.3f %9.3f %14.7lg %8.2lg\n", tstring, data[i].lon, data[i].lat, data[i].q, data[i].qerr);
  }

}

meas_data * sort_meas(meas_data *dat, long n) {
  double t[n];
  long ind[n];
  meas_data *new_data;

  for (int i=0; i<n; i++) t[i]=(double) dat[i].t;

  heapsort(t, ind, n);

  new_data=new meas_data[n];
  for (long i=0; i<n; i++) {
    new_data[i]=dat[ind[i]];
  }

  return new_data;
}

meas_data * randomize_meas(meas_data *dat, long n) {
  double t[n];
  long ind[n];
  meas_data *new_data;
  gsl_rng *rann;
  long seed;

  timeb now;
  ftime(&now);
  seed= (long) now.time + ((long) now.millitm) << 7;

  rann=gsl_rng_alloc(gsl_rng_mt19937);

  gsl_rng_set(rann, seed);

  for (long i=0; i<n; i++) t[i]=gsl_rng_uniform(rann);

  heapsort(t, ind, n);

  new_data=new meas_data[n];
  for (long i=0; i<n; i++) {
    new_data[i]=dat[ind[i]];
  }

  gsl_rng_free(rann);

  return new_data;
}

meas_data * select_meas(time_class t0, time_class t1, meas_data *dat, long n1, long *n2, int hemi) {
  time_class t[n1];
  long l1, l2;
  meas_data *new_data;

  for (long i=0; i<n1; i++) t[i]=dat[i].t;

  l1=ceil(interpolate(t, n1, t0, -1));
  if (l1 < 0) l1=0;
  l2=bin_search(t, n1, t1, -1);
  if (l2 >= n1) l2=n1-1;

  fprintf(stderr, "Selecting data between %d and %d\n", l1, l2);

  *n2=l2-l1+1;

  if (*n2 == 0) return NULL;

  new_data=new meas_data[*n2];

  if (hemi==0) {
    for (long i=l1; i<=l2; i++) {
      new_data[i-l1]=dat[i];		//if this works, it's probably a gnu ext.
    }
  } else {
    *n2=0;
    for (long i=l1; i<=l2; i++) {
      if (hemi*dat[i].lat >=0) {
        new_data[*n2]=dat[i];
	(*n2)++;
      }
    }
  }

  return new_data;
}

  meas_data *select_lat_range(meas_data *dat, long n1, float min, float max, long *n2) {
    meas_data *result;
    result=new meas_data[n1];
    *n2=0;
    if (min > max) {
      for (int i=0; i<n1; i++) {
        if (dat[i].lat >= min || dat[i].lat <= max) {
          result[*n2]=dat[i];
	  (*n2)++;
	}
      }
    } else {
      for (int i=0; i<n1; i++) {
        if (dat[i].lat >= min && dat[i].lat <= max) {
          result[*n2]=dat[i];
	  (*n2)++;
	}
      }
    }
    return result;
  }

  float correlate_meas(meas_data *samp1, 
		  meas_data *samp2, 
		  long nsamp) {
    float ave1, ave2;
    float cov, var1, var2;
    float corr;
    //calculate averages:
    ave1=0;
    ave2=0;
    for (long i=0; i<nsamp; i++) {
      ave1+=samp1[i].q;
      ave2+=samp2[i].q;
    }
    ave1/=nsamp;
    ave2/=nsamp;
    //calculate covariance and standard deviations:
    cov=0;
    var1=0;
    var2=0;
    for (long i=0; i<nsamp; i++) {
      float diff1=samp1[i].q-ave1;
      float diff2=samp2[i].q-ave2;
      cov+=diff1*diff2;
      var1+=diff1*diff1;
      var2+=diff2*diff2;
    }
    corr=cov/sqrt(var1/(nsamp-1))/sqrt(var2/(nsamp-1))/(nsamp-1);

    return corr;
  }

} //end namespace ctraj

