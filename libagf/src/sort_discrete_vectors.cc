#include <math.h>
#include <vector>

#include "agf_lib.h"

#include "peteys_tmpl_lib.h"

using namespace libpetey;
using namespace libagf;
using namespace std;

int main(int argc, char **argv) {
  FILE *fs;
  nel_ta ntrain;
  dim_ta nvar;
  real_a **train;
  cls_ta *cls;
  vector<int> **vec0;
  vector<int> *vec;
  long *sind;

  if (argc>1) fs=fopen(argv[1], "r"); else fs=stdin;
  ntrain=read_lvq(fs, train, cls, nvar);
  if (fs!=stdin) fclose(fs);

  vec=new vector<int>[ntrain];

  for (nel_ta i=0; i<ntrain; i++) {
    //printf("%d: ", i);
    vec[i].resize(nvar+1);
    //vec0[i]=new vector<int>(nvar+1);
    for (dim_ta j=0; j<nvar; j++) vec[i][j]=floor(train[i][j]);
    vec[i][nvar]=cls[i];
    //for (dim_ta j=0; j<=nvar; j++) printf(" %d", vec[i][j]);
    //printf("\n");
  }

  //vec=new vector<int>[ntrain];
  //for (nel_ta i=0; i<ntrain; i++) vec[i]=*vec0[i];

  sind=heapsort(vec, ntrain);

  if (argc>2) fs=fopen(argv[2], "w"); else fs=stdout;
  fprintf(fs, "%d\n", nvar);
  for (dim_ta j=0; j<nvar; j++) fprintf(fs, "%g ", train[sind[0]][j]);
  fprintf(fs, "%d\n", cls[sind[0]]);
  for (nel_ta i=1; i<ntrain; i++) {
    //printf("%d\n", i);
    if (vec[sind[i]]!=vec[sind[i-1]]) fprintf(fs, "\n");
    for (dim_ta j=0; j<nvar; j++) fprintf(fs, "%g ", train[sind[i]][j]);
    fprintf(fs, "%d\n", cls[sind[i]]);
    //for (dim_ta j=0; j<=nvar; j++) fprintf(fs, "%d ", vec[sind[i]][j]);
    //fprintf(fs, "\n");
  }
  if (fs!=stdout) fclose(fs);

  delete [] train;
  delete [] cls;
  delete [] vec;
  delete [] sind;

}

  
