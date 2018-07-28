#include "agf_lib.h"

using namespace libagf;

//dirt simple program, takes a list of joint or conditional probabilities
//sums their logarithms and divides by the number of samples:
//(assumes they are already importance-sampled...)
int main(int argc, char **argv) {

  cls_ta *cls;			//class labels
  real_a **p;			//joint probabilities 
  dim_ta ncls;			//number of class labels
  nel_ta *nc;			//numbers of each class
  nel_ta n;			//number of samples

  real_a hj;			//joint entropy
  real_a hc;			//conditional entropy
  real_a pt;			//total probability
  real_a hx;			//entropy of the dependent variable(s)
  real_a hi;			//entropy of the class variable

  n=read_lvq(argv[1], p, cls, ncls, 1);

  //count the numbers of each class:
  nc=new cls_ta[ncls];
  for (cls_ta i=0; i<ncls; i++) nc[i]=0;
  for (nel_ta i=0; i<n; i++) {
    nc[cls[i]]++;
  }

  hj=0;
  hc=0;
  hx=0;
  for (nel_ta i=0; i<n; i++) {
    pt=0;
    for (dim_ta j=0; j<ncls; j++) {
      pt+=p[i][j];
      if (p[i][j]!=0) hj-=nc[j]*log(p[i][j])/n;
    }
    if (pt!=0) hx-=log(pt);
    for (dim_ta j=0; j<ncls; j++) {
      if (p[i][j]!=0) hc-=nc[j]*log(p[i][j]/pt)/n;
    }
  }

  hi=0;
  for (cls_ta i=0; i<ncls; i++) {
    if (nc[i]!=0) hi-=nc[i]*log((real_a) nc[i]/n);
  }

  printf("joint entropy= %g\n", hj/log(2)/n);
  printf("cond. entropy= %g\n", hc/log(2)/n);
  printf("total entropy= %g\n", hx/log(2)/n);
  printf("class entropy= %g\n", hi/log(2)/n);

  return 0;

}
