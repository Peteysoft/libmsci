#include <stdio.h>

#include "error_codes.h"

#include "parse_command_opts.h"

#include "ctraj_defaults.h"
#include "ctraj_vfield_2d.h"

using namespace libpetey;
using namespace datasets;

namespace ctraj {

template <class real>
ctraj_vfield_2d<real>::ctraj_vfield_2d() {
  fstream=NULL;
}

//initializes the object from a binary file
template <class real>
int ctraj_vfield_2d<real>::init(char *name, 		//data file
			int64_t page_size) {		//page size for file
  long loc, dum;

  printf("Opening: %s\n", name);
  fstream=fopen(name, "r");

  all.read(fstream);
  loc=all.search_var("xgrid", dum);
  x=(simple<float> *) all.get_var(loc);
  loc=all.search_var("ygrid", dum);
  y=(simple<float> *) all.get_var(loc);
  loc=all.search_var("zgrid", dum);
  z=(simple<float> *) all.get_var(loc);
  loc=all.search_var("time", dum);
  t=(simple<time_class> *) all.get_var(loc);
  loc=all.search_var("u", dum);
  U=(dependent_swap<float> *) all.get_var(loc);
  loc=all.search_var("v", dum);
  V=(dependent_swap<float> *) all.get_var(loc);

  if (page_size > 0) {
    U->set_page_size(page_size, 1);
    V->set_page_size(page_size, 1);
  }
  return 0;

}

template <class real>
void ctraj_vfield_2d<real>::help(FILE *fs) {
  fprintf(fs, "Velocity field: 2-D Cartesian on date grid\n");
  fprintf(fs, "arg1 velocity field\n");
  ctraj_optargs(fs, "B", 1);
}

// parses the command arguments and passes the results to the 
//initialization routine
template <class real>
int ctraj_vfield_2d<real>::setup(int argc, char **argv) {
  void *optarg[10];
  int flag[10];
  int ref;                      //reference domain
  time_class t1, t2;
  int ret;
  int err=1;
  int err2;
  char *fname;
  int64_t page_size=-1;         //vfield page size in bytes

  optarg[0]=&page_size;

  ret=parse_command_opts(argc, argv, "B", "%ld", optarg, flag, OPT_WHITESPACE+2);
  if (ret < 0) {
    fprintf(stderr, "vfield_standard: error parsing command line");
    err=-1;
    ret=ret*err;
  }
  if (ret < 2) {
    fprintf(stderr, "vfield_standard: insufficient command line arguments\n");
    return -ret;
  }

  fname=argv[1];

  err2=init(fname, page_size);
  if (err2!=0) err=-1;

  for (int i=1; i<ret; i++) {
    argv[i-1]=argv[i];
  }
  argv[ret-1]=fname;

  return (ret-1)*err;
}

  
template <class real>
ctraj_vfield_2d<real>::~ctraj_vfield_2d() {
  if (fstream!=NULL) fclose(fstream);
}

template <class real>
int ctraj_vfield_2d<real>::v(int32_t domain, double t, real *x1, real *v) {
  double xind, yind;

  xind=x->interp(x1[0]);
  xind=y->interp(x1[1]);

  U->interpol(v[0], xind, yind, 0, t);
  V->interpol(v[1], xind, yind, 0, t);

  return 0;

}

template <class real>
int32_t ctraj_vfield_2d<real>::ndim() {
  return 2;
}

template <class real>
double ctraj_vfield_2d<real>::maxt() {
  return t->nel()-1;
}

template <class real>
double ctraj_vfield_2d<real>::get_tind(char *date) {
  time_class t0;
  t0.read_string(date);
  return t->interp(t0);
}

template <class real>
int ctraj_vfield_2d<real>::get_t(double tind, char *date) {
  time_class t0;
  t->get(t0, tind);
  return (t0.write_string(date)!=NULL);
}

template <class real>
double ctraj_vfield_2d<real>::get_tind(time_class date) {
  return t->interp(date);
}

template <class real>
time_class ctraj_vfield_2d<real>::get_t(double tind) {
  time_class t0;
  t->get(t0, tind);
  return t0;
}

template class ctraj_vfield_2d<float>;

} //end namespace ctraj

