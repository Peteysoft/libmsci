#include "kselect.h"
#include "randomize.h"

#include "rk_dumb_ts.h"
#include "solve_lode.h"
#include "bit_array.h"

namespace libpetey{
  template <class real>
  int test_rk_ts(int, real, real, real);
  template <class real>
  int test_rk_ts2(int, real, int);
}

//this is pretty elementary:
template <typename type>
int verify_ascending(type *list, int n) {
  for (int i=1; i<n; i++) if (list[i]<list[i-1]) return 1;
  return 0;
}

template <typename type>
int test_sorting(void (*sort)(type *, long), long n) {
  type list[n];
  //assumes a floating point type:
  for (long i=0; i<n; i++) list[i]=ranu();
  (*sort) (list, n);
  return verify_ascending(list, n);
}

template <typename type>
int test_sorting(void (*sort)(type *, long *, long), long n) {
  type list1[n];
  type list2[n];
  long ind[n];
  //assumes a floating point type:
  for (long i=0; i<n; i++) {
    list1[i]=ranu();
    //agnostic about whether or not the function also sorts the list:
    list2[i]=list[i];
  }
  (*sort) (list2, ind, n);
  for (long i=0; i<n; i++) list1[ind[i]]
  return verify_ascending(list1, n);
}

using namespace libpetey;

int main(int argc, char **argv) {

  //test kselect classes (all four of them):
  const int nk=5;
  int k[nk]={1, 2, 5, 10, 20};
  const int nn=7;
  int n[nn]={1, 2, 5, 10, 20, 50, 100};
  kiselect_base<float> *kisel;
  int err;
  int exit_code=0;

  ran_init();

  for (int i=0; i<nk; i++) {
    for (int j=0; j<nn; j++) {
      if (k[i]>n[j]) continue;
      kisel=new kiselect_naive<float>(k[i]);
      err=kisel->test(n[j]);
      if (err!=0) {
        fprintf(stderr, "kiselect_naive test failed for n=%d, k=%d\n", n[j], k[i]);
	exit_code=1;
      }
      delete kisel;
      kisel=new kiselect_tree<float>(k[i]);
      err=kisel->test(n[j]);
      if (err!=0) {
        fprintf(stderr, "kiselect_tree test failed for n=%d, k=%d\n", n[j], k[i]);
	exit_code=1;
      }
      delete kisel;
      kisel=new kiselect_heap<float>(k[i]);
      err=kisel->test(n[j]);
      if (err!=0) {
        fprintf(stderr, "kiselect_heap test failed for n=%d, k=%d\n", n[j], k[i]);
	exit_code=1;
      }
      delete kisel;
      kisel=new kiselect_quick<float>(k[i]);
      err=kisel->test(n[j]);
      if (err!=0) {
        fprintf(stderr, "kiselect_quick test failed for n=%d, k=%d\n", n[j], k[i]);
	exit_code=1;
      }
      delete kisel;
    }
  }

  //test R-K integrator:
  double h=0.01;
  int nterm=5;
  double t1=0.;
  double t2=1.;
  err=test_rk_ts<double>(nterm, t1, t2, h);
  if (err!=0) {
    fprintf(stderr, "Runge-Kutta integrator failed for n=%d, h=%g; integration limits: [%g, %g]\n", nterm, h);
    exit_code=1;
  }

  h=0.1;
  nterm=8;
  t1=0.;
  t2=100.;
  err=test_rk_ts<double>(nterm, t1, t2, h);
  if (err!=0) {
    fprintf(stderr, "Runge-Kutta integrator failed for n=%d, h=%g; integration limits: [%g, %g]\n", nterm, h);
    exit_code=1;
  }

  err=test_lode(5, 1e-10);
  //err=test_lode(5, 1e-15);
  /*
  err=test_lode(5, 1e-10);
  err=test_lode(100, 1e-10);
  err=test_lode(10000, 1e-10);
  err=test_lode(1000000, 1e-10);
  */

  err=test_rk_ts2<double>(5, 0.001, 100);
  err=test_rk_ts2<double>(5, 0.01, 1000);

  err=test_bit_array(1, 10);
  if (err!=0) exit_code=err;
  err=test_bit_array(10, 10);
  if (err!=0) exit_code=err;
  err=test_bit_array(100, 10);
  if (err!=0) exit_code=err;

  

  ran_end();

  return exit_code;
}

