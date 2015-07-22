#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <string.h>

#include "parse_command_opts.h"
#include "randomize.h"

#include "contour_anal.h"
#include "ctraj_defaults.h"

namespace ctraj {
  int cd_glob_min_unc;
  int cd_glob_max_ntrial;
  int cd_glob_hemi;
}

using namespace libpetey;
using namespace ctraj;

int main(int argc, char **argv) {
  FILE *fs;
  float *lon, *lat;
  long npt;		//total number of points
  int nf;		//number of files
  int *npf;		//number of points in each file
  int k;

  float eps1, eps2;
  int32_t neps;
  int32_t min_unc, max_trial;
  float *eps;
  float **f;
  int nside;
  int nbox;
  int32_t ntype;	//number of types to calculate
  char *Dtypes;		//the types: 0=box-counting;
			//           u=uncertainty fraction;
			//           l=box-counting using line segments

  void *optarg[20];
  int optflag[20];

  int32_t nvar;

  int logflag;
  int wrapflag=1;
  int hemi=0;

  //set defaults:
  eps1=EPS_MIN;
  eps2=EPS_MAX;
  neps=NEPS;
  min_unc=NUNC_MIN;
  max_trial=UNC_MAXN;

  optarg[0]=&eps1;
  optarg[1]=&eps2;
  optarg[2]=&neps;
  optarg[3]=&min_unc;
  optarg[4]=&max_trial;

  argc=parse_command_opts(argc, argv, "IFeuUDg-+?o", "%g%g%d%d%d%s%%%%%", optarg, optflag, OPT_WHITESPACE);

  if (optflag[9]) {
    fprintf(stdout, "syntax: contour_dimension [-D type1[type2[type3]]] [-I eps1] [-F eps2] [-e neps] [-g] [--] [-+] [file1 [file 2...]]\n");
    ctraj_optargs(stdout, "DIFeguU-+o?");
    fprintf(stdout, "**note: -u and -U options only applicable to u type;\n");
    fprintf(stdout, "        -- and -+ switches not applicable to l type\n");
    return 0;
  }
  if (optflag[5]) {
    Dtypes=(char *) optarg[5];
  } else {
    Dtypes=new char[4];
    strcpy(Dtypes, "u");
  }
  ntype=strlen(Dtypes);

  logflag=optflag[6];
  if (optflag[7]) hemi=-1;
  if (optflag[8]) hemi=1;
  if (optflag[10]) wrapflag=0;

  if (argc==1) {
    npt=read_ascii_contour(stdin, lon, lat);
  } else {
    nf=argc-1;
    npf=new int[nf];
    npt=0;
    for (int i=0; i<nf; i++) {
      fs=fopen(argv[i+1], "r");
      fread(&nvar, sizeof(int32_t), 1, fs);
      assert(nvar==2);
      fseek(fs, 0, SEEK_END);
      npf[i]=(ftell(fs)-sizeof(nvar))/nvar/sizeof(float);
      npt+=npf[i];
      fclose(fs);
    }
    lon=new float[npt];
    lat=new float[npt];

    k=0;
    for (int i=0; i<nf; i++) {
      fs=fopen(argv[i+1], "r");
      fseek(fs, sizeof(nvar), SEEK_SET);
      for (int j=0; j<npf[i]; j++) {
        fread(lon+k+j, sizeof(float), 1, fs);
        fread(lat+k+j, sizeof(float), 1, fs);
      }
      k+=npf[i];
      fclose(fs);
    }
    delete [] npf;

  }

  eps=new float[neps];
  f=new float *[ntype];
  f[0]=new float[ntype*neps];
  for (int i=1; i<ntype; i++) f[i]=f[0]+i*neps;

  ran_init();
  for (int i=0; i<neps; i++) {
    if (logflag) {
      eps[i]=eps1*pow(eps2/eps1, 1.0*i/(neps-1));
    } else {
      eps[i]=eps1+i*(eps2-eps1)/(neps-1);
    }
  }

  for (int i=0; i<neps; i++) {
    printf("%12.6g", eps[i]);
    for (int j=0; j<ntype; j++) {
      switch (Dtypes[j]) {
        case ('u'):
          f[j][i]=uncertainty_fraction(lon, lat, npt, eps[i], hemi, min_unc, max_trial, wrapflag);
          break;
        case ('0'):
          f[j][i]=(float) count_boxes_global(lon, lat, npt, eps[i], hemi, wrapflag);
          break;
        case ('l'):
          f[j][i]=(float) measure_contour(lon, lat, npt, eps[i], wrapflag);
          break;
      }
      //printf("%12.6g %12.6g\n", eps[i], f[i]);
      printf(" %12.6g", f[j][i]);
    }
    printf("\n");
  }
  ran_end();

//  for (int i=0; i<neps; i++) {
//    printf("%12.6g %12.6g\n", eps[i], ftotal[i]);
//  }
  printf("\n");

  delete [] lon;
  delete [] lat;


  delete [] eps;
  delete [] f[0];
  delete [] f;

}

