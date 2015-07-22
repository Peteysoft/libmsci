#include <math.h>
#include <stdio.h>
#include <sys/timeb.h>

#include "gbin_ll.h"
//#include "lonlat_binsub.h"
#include "global_metric.h"
#include "azeq_binsub.h"
#include "meas_data.h"
#include "randomize.h"

int main() {
  meas_data *loc;
  meas_data *test;
  gbin_ll<int64_t, lonlat_coord, azeq_binsub> *bins;
  global_bin<int64_t, lonlat_coord> *el;
  global_bin<int64_t, lonlat_coord> *nnel;
  time_class t0, tf;
  long minn=1000;
  long maxn=100000;
  long nn=5;
  long ntrial=100;
  long n;
  lonlat_coord test2;

  timeb t1, t2;
  float ttot1[nn], ttot2[nn];
  float setup_t, bs_t, nbs_t;

  float d, dmin;
  int64_t indm;

  float dmin2;

  lonlat_coord *v;

  t0="2000/1/1";
  tf="2000/1/2";

  ran_init();

  for (long i=0; i<nn; i++) {
    n=minn*pow(pow(1.*maxn/minn, 1./(nn-1)), i);
    loc=generate_random_global(t0, tf, n);
    ftime(&t1);
    bins=new gbin_ll<int64_t, lonlat_coord, azeq_binsub>();
    //v=new lonlat_coord[n];
    for (long j=0; j<n; j++) {
      el=new global_bin<int64_t, lonlat_coord>;
      el->coords.lon=loc[j].lon;
      el->coords.lat=loc[j].lat;
      el->ind=j;
      bins->add_el(el);
      //v[j].lon=loc[j].lon;
      //v[j].lat=loc[j].lat;
    }
    //bins=new gbin_ll<int64_t, lonlat_coord, azeq_binsub>(v, n);
    //delete [] v;
    bins->update();
    ftime(&t2);
    setup_t=t2.time-t1.time+(t2.millitm/1000.-t1.millitm/1000);
    ttot1[i]=setup_t;
    ttot2[i]=0;

    for (long j=0; j<ntrial; j++) {
      test=generate_random_global(t0, tf, ntrial);
      ftime(&t1);
      test2.lon=test[j].lon;
      test2.lat=test[j].lat;
      dmin2=bins->nn(test2, nnel);
      ftime(&t2);
      bs_t=t2.time-t1.time+(t2.millitm/1000.-t1.millitm/1000);
      ttot1[i]+=bs_t;

      ftime(&t1);
      dmin=sdist(loc[0].lon, loc[0].lat, test[j].lon, test[j].lat);
      indm=0;
      for (long k=1; k<n; k++) {
        d=sdist(loc[k].lon, loc[k].lat, test[j].lon, test[j].lat);
	if (d < dmin) {
          dmin=d;
	  indm=k;
	}
      }
      ftime(&t2);
      nbs_t=t2.time-t1.time+(t2.millitm/1000.-t1.millitm/1000);
      ttot2[i]+=nbs_t;

      printf("%d %d %f %f %f %ld (%f, %f)\n", j, n, setup_t, bs_t, sqrt(dmin2), nnel->ind, nnel->coords.lon, nnel->coords.lat);
      printf("%d %d %f %f %f %d (%f, %f)\n", j, n, 0., nbs_t, sqrt(dmin), indm, loc[indm].lon, loc[indm].lat);
      if (nnel->ind != indm) {
        printf("******************************\n");
        printf("* ERROR * ERROR * ERROR!! ?? *\n");
        printf("******************************\n");
      }
      delete [] test;
    }
    delete bins;
    delete [] loc;

  }

  for (int32_t i=0; i<nn; i++) {
    printf("%f %f\n", ttot1[i], ttot2[i]);
  }
  ran_end();

}
        
