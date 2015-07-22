#include <stdio.h>
#include <stdlib.h>
#include <dlfcn.h>

#include "error_codes.h"
#include "parse_command_opts.h"

#include "ctraj_defaults.h"

#include "ctraj_vfield_anal.h"

using namespace libpetey;

namespace ctraj {

template <class real>
ctraj_vfield_anal<real>::ctraj_vfield_anal(int32_t nd) {
  ndims=nd;
  vfield=NULL;
  handle=NULL;
  nparam=0;
  param=NULL;
}

template <class real>
ctraj_vfield_anal<real>::~ctraj_vfield_anal() {
  if (handle!=NULL) dlclose(handle);
  if (param!=NULL) delete [] (double *) param;
}

template <class real>
int ctraj_vfield_anal<real>::init(char *file, char *sym, void *p) {
  char *err;

  param=p;

  handle=dlopen(file, RTLD_NOW);
  if (handle==NULL) {
    err=dlerror();
    fprintf(stderr, "ctraj_vfield_anal: dlopen returned the following error message:\n");
    fprintf(stderr, "%s\n", err);
    exit(FILE_READ_ERROR);
  }

  vfield=(int (*) (double, float *, float *, void *)) dlsym(handle, sym);
  if (handle==NULL) {
    err=dlerror();
    fprintf(stderr, "ctraj_vfield_anal: dlsym returned the following error message:\n");
    fprintf(stderr, "%s\n", err);
    exit(FILE_READ_ERROR);
  }
  return 0;

}

template <class real>
void ctraj_vfield_anal<real>::help(FILE *fs) {
  fprintf(fs, "Velocity field: analytical\n");
  fprintf(fs, "    The user must write and compile a function of the form:\n");
  fprintf(fs, "    int v(double t, float *x, float *v, void *param);\n");
  fprintf(fs, "    and specify the shareable object library and entry point.\n");
  fprintf(fs, "where:\n");
  ctraj_optargs(fs, "t", 1);
  fprintf(fs, "  arg1 object file\n");
  fprintf(fs, "  arg2 entry point\n");
  fprintf(fs, "  [arg3] first parameter\n");
  fprintf(fs, "  [arg4] second parameter\n");
  fprintf(fs, "  [arg5] ...\n");
}

template <class real>
int ctraj_vfield_anal<real>::setup(int argc, char **argv) {
  char *fname;
  char *sym;
  void *optarg[10];
  int flag[10];
  int ret, err, err2;

  optarg[1]=&nparam;

  ret=parse_command_opts(argc, argv, "t", "%d", optarg, flag, OPT_WHITESPACE+2);
  if (ret < 0) {
    fprintf(stderr, "vfield_anal: error parsing command line");
    err=-1;
    ret=ret*err;
  }
  if (ret < 3) {
    fprintf(stderr, "vfield_standard: insufficient command line arguments\n");
    return -3;
  }

  fname=argv[1];
  sym=argv[2];
  if (nparam>0) param=new double[nparam];

  for (int i=3; i<3+nparam; i++) {
    sscanf(argv[i], "%lg", (double *) param+i-3);
  }

  err2=init(fname, sym, param);
  if (err2!=0) err=-1;

  for (int i=2+nparam; i<ret; i++) {
    argv[i-2+nparam]=argv[i];
  }
  argv[ret-2]=fname;
  argv[ret-1]=sym;

  return (ret-1)*err;
}

template <class real>
int ctraj_vfield_anal<real>::v(int32_t domain, double t, real *x, real *v) {
  
  return (*vfield) (t, x, v, param);

}

template <class real>
int32_t ctraj_vfield_anal<real>::ndim() {
  return ndims;
}

template <class real>
double ctraj_vfield_anal<real>::maxt() {
  return -1;
}

template class ctraj_vfield_anal<float>;

} //end namespace ctraj

