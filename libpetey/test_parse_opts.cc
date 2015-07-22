#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "parse_command_opts.h"

using namespace std;
using namespace libpetey;

namespace libpetey {extern void get_format(const char *, int, char *);}

int main(int argc, char **argv) {
  void *optarg[20];
  int flags[20];
  int nopts;
  int err;
  int opts=0;

  err=parse_command_opts(argc, argv, "pfbh", "%s%s%%", optarg, flags, 3);
  for (int i=0; i<argc; i++) printf("%s ", argv[i]);
  printf("\n");
  if (flags[2]==1) opts=1;

  if (flags[3]) {
    printf("test_parse_opts [-g<number>] [-i<integer>] [-s<string>] [-c<char>] \n");
    printf("                [-d<double>] [-l<long>] [-h] [-b] [args]\n");
    printf("-g read a floating point number and print it out\n");
    printf("-i read an integer and print it out\n");
    printf("-s read a string and print it out\n");
    printf("-c read a character and print it out\n");
    printf("-d read a double-precision floating point number and print it out\n");
    printf("-l read a longword integer and print it out\n");
    printf("-h print this help screen\n");
    printf("-b accept whitespace between option and argument\n");
    printf("-p custom string of option letters\n");
    printf("-f format codes used in conjunction with -p\n");
    printf("\n");
  }

  if (flags[0] && flags[1]) {
    char *code=(char *) optarg[0];
    char *fstring=(char *) optarg[1];
    char form[3];
    int n=strlen(code);
    void *p2[n];
    int f2[n];

    for (int i=0; i<n; i++) {
      get_format(fstring, i, form);
      if (form[1]=='c') {
        p2[i]=new char;
      } else if (form[1]=='d') {
        p2[i]=new int32_t;
      } else if (form[1]=='g') {
        p2[i]=new float;
      } else if (form[1]=='l') {
        if (form[2]=='d') {
          p2[i]=new int64_t;
        } else if (form[2]=='g') {
          p2[i]=new double;
        }
      }
    }

    err=parse_command_opts(abs(err), argv, code, fstring, p2, f2, opts);

    for (int i=0; i<n; i++) {
      if (f2[i]) {
        get_format(fstring, i, form);
        if (form[1]=='s') {
          printf("-%c (string)=%s\n", code[i], (char *) p2[i]);
	  delete [] (char *) p2[i];
	} else if (form[1]=='c') {
          printf("-%c (char)=%c\n", code[i], * (char *) p2[i]);
	  delete (char *) p2[i];
        } else if (form[1]=='d') {
          printf("-%c (integer)=%d\n", code[i], * (int32_t *) p2[i]);
	  delete (int32_t *) p2[i];
        } else if (form[1]=='g') {
          printf("-%c (float)=%g\n", code[i], * (float *) p2[i]);
	  delete (float *) p2[i];
        } else if (form[1]=='l') {
          if (form[2]=='d') {
            printf("-%c (long int)=%ld\n", code[i], * (int64_t *) p2[i]);
	    delete (int64_t *) p2[i];
          } else if (form[2]=='g') {
            printf("-%c (double)=%lg\n", code[i], * (double *) p2[i]);
	    delete (double *) p2[i];
          }
        }
      }
    }
  } else {
    float g=0;
    int i=0;
    char c;
    int64_t l;
    double d;
    char *s;

    optarg[0]=(void *) &g;
    optarg[1]=(void *) &i;
    optarg[2]=NULL;
    optarg[3]=(void *) &c;
    optarg[4]=(void *) &d;
    optarg[5]=(void *) &l;

    err=parse_command_opts(abs(err), argv, "giscdl", "%g%d%s%c%lg%ld", optarg, flags, opts);

    s=(char *) optarg[2];
  
    if (flags[0]) printf("-g (float)=%g\n", g);
    if (flags[1]) printf("-i (integer)=%d\n", i);
    if (flags[2]) printf("-s (string)=%s\n", s);
    if (flags[3]) printf("-c (char)=%c\n", c);
    if (flags[4]) printf("-d (double)=%lg\n", d);
    if (flags[5]) printf("-l (long)=%ld\n", l);
  }
  
  //printf("%d%d%d%d\n", flags[0], flags[1], flags[2], flags[3]);

  printf("Arguments: ");
  for (int i=0; i<abs(err); i++) printf("%s ", argv[i]);
  printf("\n");

  return err;

}

