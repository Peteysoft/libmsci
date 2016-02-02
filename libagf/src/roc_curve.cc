#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <string.h>

#include <getopt.h>

#include "parse_command_opts.h"
#include "agf_lib.h"

#define NTRIAL_DEFAULT 21

using namespace std;
using namespace libpetey;
using namespace libagf;

int main(int argc, char **argv) {
  int err;
  FILE *fs;
  int32_t n=NTRIAL_DEFAULT;

  char *fbase;
  char *fname=NULL;
  int baselen;
  char *outbase0;
  char **outbase;
  char *outname;

  real_a vmin, vmax;
  real_a thresh;
  int gflag;
  int Dflag;
  real_a r_val=0;
  void *optarg[20];
  int flag[20];
  char *command;

  real_a *data=NULL;
  real_a dmin, dmax;
  nel_ta nel;
  cls_ta *truth=NULL;
  cls_ta *ret=NULL;
  cls_ta ct;
  real_a *con=NULL;
  cls_ta *ret0=NULL;
  cls_ta *ret1=NULL;

  nel_ta n0, n1, nhit, nmiss;
  real_a *hitrate, *missrate;

  int32_t algtyp=0;
  int strpos;

  optarg[0]=&vmin;
  optarg[1]=&vmax;
  optarg[4]=&n;
  optarg[5]=&r_val;
  optarg[6]=&algtyp;

  //this whole business is a bit silly:
  argc=parse_command_opts(argc, argv, "01Dgqrc", "%g%g%%%d%g%d", optarg, flag, 3);
  Dflag=flag[2];
  gflag=flag[3];

  while (getopt(argc, argv, "a:d:f:h:i:I:k:l:m:N:s:S:t:v:V:W:Knu")!=-1);

  if (argc<2) {
    printf("Generates the ROC-curve [1] in three possible ways:\n");
    printf("1. by converting floating-point data **\n");
    printf("2. by shifting the discrimination border using the -r option *** (-D option)\n");
    printf("3. by shifting the discrimination border based on calculated\n");
    printf("  conditional probabilities (-D option, -r option must be set ***)\n");
    printf("\n");
    printf("syntax: roc_curve [options] [-D] [-c type] [-q n] [-g] train\n");
    printf("\n");
    printf("where:\n");
    printf("  train   is the base-name for the training data files\n");
    printf("\n");
    printf("options:\n");
    printf("  -0  min  = minimum threshold for continuous data\n");
    printf("  -1  max  = maximum threshold for continuous data\n");
    printf("  -c  type = is the type of algorithm to use (0=AGF borders, 1=AGF, 2=KNN)\n");
    printf("               --default AGF borders\n");
    printf("  -q    nt = is the number trials [%d]\n", NTRIAL_DEFAULT);
    printf("  -g       = use a geometric progression\n");
    printf("  -D       = generate ROC curve by shifting discrimination border\n");
    printf("\n");
    printf("  -a normf = file containing normalization data (input/output)\n");
    printf("  -d ndiv  = number of divisions in data (default=%d)\n", (int32_t) (1/F_DEFAULT));
    printf("  -i maxi1 =  maximum number of iterations when searching for class border (%d)\n", (int32_t) agf_global_borders_maxiter);

    printf("  -I maxi2 = maximum number of iterations when calculating weights (%d, %d)\n", (int32_t) agf_global_weights_maxiter, (int32_t) agf_global_borders_maxiter);
    printf("  -k k     = number of nearest neighbours\n");
    printf("  -K       = keep temporary files\n");
    printf("  -l tol   = tolerance of W (default=%g)\n", (float) agf_global_weights_tol);
    printf("  -n       = option to normalise the data\n");
    printf("  -N maxi3 = maximum number of iterations in supernewton (default=%d, %d)\n", (int32_t) agf_global_weights_maxiter, (int32_t) agf_global_borders_maxiter);
    printf("  -r r0    = location of discrimination border (default=0)\n");
    printf("  -s n     = number of times to sample the border (default=%d)\n", (int32_t) NBORD_DEFAULT);
    printf("  -S nsv   = number of singular values from SVD\n");
    printf("  -t tol   = tolerance of border samples (default=%g)\n", (float) TOL_DEFAULT);
    printf("  -u       = store borders data in un-normalized coordinates\n");
    printf("  -v var1  = lower filter variance bound\n");
    printf("               --default is to use the total variance of the data/n^(2/D)\n");
    printf("  -V var2  = lower filter variance bound\n");
    printf("               --default is to use the total variance of the data\n");
    printf("  -W Wc    = objective total weight (default=%g)\n", (float) W_DEFAULT);
    //printf("  -r thresh= discrimination border for AGF borders\n");
    //command=new char[35];
    //sprintf(command, "nfold | tail -n 17 | head -n 14");
    //system(command);
    printf("\n");
    printf("Notes:\n");
    printf("[1] Jolliffe & Stephenson, eds., 2003. Forecast Verification:\n");
    printf("    a Practitioner's Guide in the Atmospheric Sciences. Wiley\n");
    printf("\n");
    printf("**  overwrites <train>.cls if it exists\n");
    printf("\n");
    printf("*** -r option has no effect in AGF and KNN algorithms (type=1; type=2)\n");
    printf("    --simply set it to 0 for the -D option to take effect\n");
    //delete [] command;
    return 1;
  }

  if (algtyp==-1) algtyp=0;

  fbase=argv[optind];
  baselen=strlen(fbase);
  if (argc-optind > 1) {
    outbase0=argv[optind+1];
  } else {
    outbase0=new char[baselen+20];
    sprintf(outbase0, "%s%11.11ld", fbase, seed_from_clock());
  }

  fname=new char[baselen+5];
  if (Dflag==0) {
    sprintf(fname, "%s.dat", fbase);
    data=read_datfile<real_a>(fname, nel);
    if (data==NULL) {
      fprintf(stderr, "Unable to open file, %s, for reading\n", fname);
      return UNABLE_TO_OPEN_FILE_FOR_READING;
    }
    truth=new cls_ta[nel];
    //for (nel_ta i=0; i<nel; i++) printf("%f ", data[i]);
  } else {
    sprintf(fname, "%s.cls", fbase);
    truth=read_clsfile<cls_ta>(fname, nel);
  }

  if (flag[0]==0 || flag[1]==0) {
    if (Dflag) {
      if (flag[0]==0) vmin=-1.;
      if (flag[1]==0) vmax=1.;
    } else {
      dmin=data[0];
      dmax=data[1];
      for (nel_ta i=1; i<nel; i++) {
        if (data[i]<dmin) dmin=data[i];
        else if (data[i]>dmax) dmax=data[i];
      }
      if (flag[0]==0) vmin=dmin;
      if (flag[1]==0) vmax=dmax;
    }
    if (flag[0]==0 && flag[5]==0) dmin=vmin+(vmax-vmin)/(n+1); else dmin=vmin;
    if (flag[1]==0 && flag[5]==0) dmax=vmax-(vmax-vmin)/(n+1); else dmax=vmax;
    vmin=dmin;
    vmax=dmax;
  }

  printf("%d\n", n);

  command=new char[10000];
  sprintf(command, "nfold ");
  strpos=6;
  for (int i=1; i<optind; i++) {
    strcat(command, argv[i]);
    strcat(command, " ");
    strpos+=strlen(argv[i])+1;
  }

  outbase=new char*[n];
  outname=new char[strlen(outbase0)+20];
  for (int32_t i=0; i<n; i++) {
    outbase[i]=new char[strlen(outbase0)+12];
    sprintf(outbase[i], "%s-%3.3d", outbase0, i);
  }

  if (Dflag && flag[5]) {
    //we only need to run the validation routine once:
    sprintf(command+strpos, " -r %g -c %d %s %s", r_val, algtyp, fbase, outbase0);
    printf("%s\n", command);
    err=system(command);
    if (err!=0) exit(err);
    sprintf(outname, "%s.cls", outbase0);
    ret0=read_clsfile<cls_ta>(outname, nel);
    sprintf(outname, "%s.con", outbase0);
    con=read_datfile<real_a>(outname, nel);
    ret1=new cls_ta[nel];
  }

  for (int32_t i=0; i<n; i++) {
    if (gflag) {
      thresh=vmin*pow((vmax/vmin), 1.0*i/(n-1));
    } else {
      thresh=vmin+(vmax-vmin)*i/(n-1);
    }
    if (Dflag) {
      r_val=thresh;
      if (flag[5]) {
        for (nel_ta j=0; j<nel; j++) {
          if (ret0[j]>0) r_val=con[j]; else r_val=-con[j];
          if (r_val>thresh) ret1[j]=1; else ret1[j]=0;
          //ret1[j]=(cls_ta) (r_val>thresh);
          //printf("%d %g %d %g\n", ret[j], con[j], ret1[j], r_val);
        }
        sprintf(outname, "%s.cls", outbase[i]);
        fs=fopen(outname, "w");
        fwrite(ret1, sizeof(cls_ta), nel, fs);
        fclose(fs);
      }
    } else {
      ct=0;
      for (nel_ta j=0; j<nel; j++) {
        truth[j]=(cls_ta) (data[j] > thresh);
        //printf("%g %g %d\n", data[j], thresh, truth[j]);
        ct+=truth[j];
      }
      if (ct==0 || ct==nel) continue;
      sprintf(fname, "%s.cls", fbase);
      fs=fopen(fname, "w");
      fwrite(truth, sizeof(cls_ta), nel, fs);
      fclose(fs);
    }

    if (Dflag==0 || flag[5]==0) {
      sprintf(command+strpos, " -r %g -c %d %s %s", r_val, algtyp, fbase, outbase[i]);
      printf("%s\n", command);
      err=system(command);
      if (err!=0) exit(err);
    }
  }

  hitrate=new real_a[n];
  missrate=new real_a[n]; 
  for (int32_t i=0; i<n; i++) {

    if (Dflag==0) {
      ct=0;
      for (nel_ta j=0; j<nel; j++) {
        truth[j]=(cls_ta) (data[j] > thresh);
        ct+=truth[j];
      }
      if (ct==0 || ct==nel) continue;
    }

    sprintf(outname, "%s.cls", outbase[i]);
    ret=read_clsfile<cls_ta>(outname, nel);
    
    n0=0;
    n1=0;
    nhit=0;
    nmiss=0;
    for (nel_ta j=0; j<nel; j++) {
      if (truth[j]>0) {
        n1++;
        if (ret[j]>0) nhit++;
      } else {
        n0++;
        if (ret[j]>0) nmiss++;
      }
    }
    hitrate[i] =1.*nhit/n1;
    missrate[i]=1.*nmiss/n0;
    //printf("%d %d %d %d\n", n0, n1, nhit, nmiss);
  }

  if (argc-optind == 1) {
    for (int32_t i=0; i<n; i++) {
      sprintf(command, "rm %s.cls", outbase[i]);
      printf("%s\n", command);
      system(command);
      if (Dflag==0 || flag[5]==0) {
        sprintf(command, "rm %s.con", outbase[i]);
        printf("%s\n", command);
        system(command);
      }
    }
    if (Dflag && flag[5]) {
      sprintf(command, "rm %s.cls", outbase0);
      printf("%s\n", command);
      system(command);
      sprintf(command, "rm %s.con", outbase0);
      printf("%s\n", command);
      system(command);
    }
    delete [] outbase0;
  }


  printf("\n");
  for (int32_t i=0; i<n; i++) {
    if (gflag) {
      thresh=vmin*pow((vmax/vmin), 1.0*i/(n-1));
    } else {
      thresh=vmin+(vmax-vmin)*i/(n-1);
    }

    printf("%g %g %g\n", thresh, missrate[i], hitrate[i]);

    delete [] outbase[i];
  }

  delete [] outbase;
  delete [] hitrate;
  delete [] missrate;
  delete [] command;
  delete [] outname;
  if (fname!=NULL) delete [] fname;
  if (ret!=NULL) delete [] ret;
  if (ret0!=NULL) delete [] ret0;
  if (ret1!=NULL) delete [] ret1;
  if (truth!=NULL) delete [] truth;
  if (con!=NULL) delete [] con;
  if (data!=NULL) delete [] data;

  return 0;
}
     

