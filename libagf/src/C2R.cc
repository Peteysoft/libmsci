#include <assert.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include "parse_command_opts.h"

#include "agf_lib.h"

using namespace std;
using namespace libagf;
using namespace libpetey;

int main(int argc, char **argv) {
  FILE *fs;  

  real_a *con;
  cls_ta *cls;
  nel_ta n1, n2;

  char *agfbase;
  char *rfile;

  char *confile;
  char *clsfile;

  real_a *r;
  real_a r0;

  void *optarg[10];
  int flag[10];

  optarg[0]=&r0;
  argc=parse_command_opts(argc, argv, "r", "%g", optarg, flag, 1);

  if (argc < 3) {
    printf("Converts classes and confidence ratings into a \n");
    printf("a difference in conditional probabilities:\n");
    printf("     R = 2*C*(c - 1);   (C = |R|)\n\n");
    printf("usage: C2R basename outfile\n\n");
    printf("    where:\n");
    printf("basename = base name of binary files:\n");
    printf("           .cls for classes, .con for confidence ratings\n");
    printf("outfile  = name of output file\n");
    printf("\n");
    printf("-r r0    = reverse the operation.  Must supply threshold\n");
    exit(INSUFFICIENT_COMMAND_ARGS);
  }

  if (flag[0]) {
    rfile=argv[1];
    agfbase=argv[2];
  } else {
    agfbase=argv[1];
    rfile=argv[2];
  }

  //generate the agf class file names:
  clsfile=new char[strlen(agfbase)+5];
  strcpy(clsfile, agfbase);
  strcat(clsfile, ".cls");

  confile=new char[strlen(agfbase)+5];
  strcpy(confile, agfbase);
  strcat(confile, ".con");

  //printf("%s; %s\n", confile, clsfile);

  if (flag[0]) {
    r=read_datfile(rfile, n1);
    cls=new cls_ta[n1];
    con=new real_a[n1];
    for (nel_ta i=0; i<n1; i++) {
      if (r[i]<r0) {
        cls[i]=0;
        con[i]=(r0-r[i])/(1+r0);
      } else {
        cls[i]=1;
        con[i]=(r[i]-r0)/(1-r0);
      }
      //printf("%g %d %g\n", r[i], cls[i], con[i]);
    }
    fs=fopen(clsfile, "w");
    fwrite(cls, sizeof(cls_ta), n1, fs);
    fclose(fs);
    fs=fopen(confile, "w");
    fwrite(con, sizeof(real_a), n1, fs);
    fclose(fs);
  } else {
    cls=read_clsfile(clsfile, n1);
    con=read_datfile(confile, n2);

    if (n1!=n2) {
      fprintf(stderr, "C2R: sample count mismatch: %d elements in %s\n", 
		n1, clsfile);
      fprintf(stderr, "                            %d elements in %s\n",
		n2, confile);
      exit(SAMPLE_COUNT_MISMATCH);
    }

    r=new real_a[n1];
    for (nel_ta i=0; i<n1; i++) {
      if (cls[i] < 1) {
        r[i]=-con[i];
      } else {
        r[i]=con[i];
      }
      //printf("%f\n", r[i]);
    }
    //printf("ii\n");
    fs=fopen(rfile, "w");
    fwrite(r, sizeof(real_a), n1, fs);
    fclose(fs);
  }

  //printf("iii\n");
  delete [] con;
  delete [] cls;
  delete [] r;

}

