#include "randomize.h"

#include "agf_lib.h"

using namespace libagf;

int main(int argc, char **argv) {
  const int ntest1=11;
  int n1[ntest1]={1, 1, 1,  2,  5, 1, 2, 10, 20, 50, 100};
  int n2[ntest1]={1, 2, 10, 2, 15, 1, 1, 1, 20, 15, 100};
  const int ntest2=6;
  int n[ntest2]={2, 3, 5, 10, 20, 100};
  int err;
  int exit_code=0;

  ran_init();

  for (int i=0; i<ntest1; i++) {
    err=test_oppositesample<real_a>(n1[i], n2[i]);
    if (err!=0) {
      fprintf(stderr, "test_opposite sample failed for (%d, %d)\n", n1[i], n2[i]);
      exit_code=-1;
    }
  }

  for (int i=0; i<ntest2; i++) {
    err=test_oppositesample_small<real_a>(n[i]);
    if (err!=0) {
      fprintf(stderr, "test_opposite sample small failed for (%d)\n", n[i]);
      exit_code=-1;
    }
  }

  for (int i=0; i<10; i++) {
    err=test_parse_multi_partitions<int>(5, 5);
    if (err!=0) {
      fprintf(stderr, "test_parse_multi_partitions failed\n");
      exit_code=err;
    }
    err=test_parse_multi_partitions<int>(100, 100);
    if (err!=0) {
      fprintf(stderr, "test_parse_multi_partitions failed\n");
      exit_code=err;
    }
  }

  ran_end();

  return exit_code;
}

