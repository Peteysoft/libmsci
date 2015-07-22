#include <math.h>
#include <assert.h>
#include <string.h>

#include "agf_lib.h"

#include "vector_s.h"
#include "peteys_tmpl_lib.h"
#include "randomize.h"

using namespace libpetey;
using namespace libagf;

int main(int argc, char **argv) {
  char *outfile;
  char *confile;
  FILE *fs;
  nel_ta ntrain;
  dim_ta nvar;
  real_a **train;
  cls_ta *cls;
  vector_s<int> *vec;

  real_a *minvec;
  real_a *maxvec;
  real_a min;
  real_a max;

  real_a **test;
  cls_ta *testcls;
  vector_s<int> *tvec;
  dim_ta nvar2;
  nel_ta ntest;

  long ind;

  real_a *p;
  cls_ta ncls=1;
  cls_ta *result;
  real_a *con;
  int nmatch;

  ran_init();

  fs=fopen(argv[1], "r");
  ntrain=read_lvq(fs, train, cls, nvar);
  if (fs!=stdin) fclose(fs);

  vec=new vector_s<int>[ntrain];
  minvec=new real_a[nvar+1];
  maxvec=new real_a[nvar+1];
  for (dim_ta j=0; j<nvar; j++) {
    minvec[j]=train[0][j];
    maxvec[j]=train[0][j];
  }
  minvec[nvar]=cls[0];
  maxvec[nvar]=cls[0];
  for (nel_ta i=0; i<ntrain; i++) {
    if (cls[i]>=ncls) ncls=cls[i]+1;
    vec[i].resize(nvar+1);
    //vec0[i]=new vector_s<int>(nvar+1);
    for (dim_ta j=0; j<nvar; j++) {
      vec[i][j]=floor(train[i][j]);
      if (train[i][j]<minvec[j]) minvec[j]=train[i][j];
      		else if (train[i][j]>maxvec[j]) maxvec[j]=train[i][j];
    }
    vec[i][nvar]=cls[i];
    if (cls[i]<minvec[nvar]) minvec[nvar]=cls[i];
		else if (cls[i]>maxvec[nvar]) maxvec[nvar]=cls[i];
  }

  fs=fopen(argv[2], "r");
  ntest=read_lvq(fs, test, testcls, nvar2);
  if (fs!=stdout) fclose(fs);
  assert(nvar2<=nvar);

  for (nel_ta i=0; i<ntest; i++) {
    for (dim_ta j=0; j<nvar2; j++) {
      if (test[i][j]<minvec[j]) minvec[j]=test[i][j];
      		else if (test[i][j]>maxvec[j]) maxvec[j]=test[i][j];
    }
  }
  min=minvec[0];
  max=maxvec[0];
  for (dim_ta j=1; j<nvar; j++) {
    if (minvec[j]<min) min=minvec[j]; else if (maxvec[j]>max) max=maxvec[j];
  }

  fprintf(stderr, "min      max\n");
  for (dim_ta j=0; j<nvar; j++) fprintf(stderr, "%g %g\n", minvec[j], maxvec[j]);

  fs=stdout;
  tvec=new vector_s<int>[ntest];
  p=new real_a[ncls];
  result=new cls_ta[ntest];
  con=new real_a[ntest];
  for (nel_ta i=0; i<ntest; i++) {
    long lastind=-1;
    dim_ta vused=nvar;

    tvec[i].resize(nvar+1, min-2);
    for (dim_ta j=0; j<nvar; j++) {
      if (tvec[i][j]!=NAN) tvec[i][j]=floor(test[i][j]);
    }
    do {
      tvec[i][vused]=tvec[i].missing;
      ind=bin_search(vec, ntrain, tvec[i], lastind);
      vused--;
    } while (tvec[i]!=vec[ind] && vused>=0);
    nmatch=0;
    for (dim_ta j=0; j<ncls; j++) p[j]=0;
    for (nel_ta k=ind; k<ntrain && vec[k]==tvec[i]; k++) {
      //for (dim_ta j=0; j<nvar; j++) fprintf(fs, "%g ", train[k][j]);
      //fprintf(fs, "%d\n", cls[k]);
      p[cls[k]]++;
      nmatch++;
    }
    for (nel_ta k=ind-1; k>=0 && vec[k]==tvec[i]; k--) {
      //for (dim_ta j=0; j<nvar; j++) fprintf(fs, "%g ", train[k][j]);
      //fprintf(fs, "%d\n", cls[k]);
      p[cls[k]]++;
      nmatch++;
    }
    result[i]=choose_class(p, ncls);
    for (dim_ta j=1; j<ncls; j++) {
      p[j]=p[j]/nmatch;
      fprintf(fs, "%g ", p[j]);
    }
    fprintf(fs, "%d\n", result[i]);
    con[i]=(ncls*p[result[i]]-1)/(ncls-1);
  }

  outfile=new char[strlen(argv[3])+5];
  sprintf(outfile, "%s.cls", argv[3]);
  fs=fopen(outfile, "w");
  fwrite(result, sizeof(cls_ta), ntest, fs);
  fclose(fs);

  confile=new char[strlen(argv[3])+5];
  sprintf(outfile, "%s.con", argv[3]);
  fs=fopen(outfile, "w");
  fwrite(con, sizeof(real_a), ntest, fs);
  fclose(fs);

  delete [] outfile;
  delete [] confile;

  ran_end();

  delete [] train;
  delete [] cls;
  delete [] vec;

  delete [] test;
  delete [] tvec;
  delete [] testcls;

  delete [] minvec;
  delete [] maxvec;

  delete [] p;
  delete [] result;

  delete [] con;

}

  
