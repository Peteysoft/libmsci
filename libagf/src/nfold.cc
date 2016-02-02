#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include <sys/timeb.h>
#include <unistd.h>

#include "agf_lib.h"

using namespace std;
using namespace libagf;

int main(int argc, char *argv[]) {
  FILE *fs, *fs2, *fs3;
  FILE *diagfs;

  char *basename;		//base name of input files
  char *controlfile;		//name of control input file

  char *tempbase;		//base name for temporary files
  char *trainbase;		//base name of separated training data
  char *outbase;		//base name for optional output file
  char *trainvecfile;		//for finding the class borders

  char *trainclsfile;		//	"
  char *testbase;
  char *testfile;		//test data
  char *outfile;		//output classes
  char *confile;		//output confidence ratings
  char *brdfile;		//contains class borders

  nel_ta nsamp;		//total number of samples
  dim_ta nvar;		//number of variables
  cls_ta nclass;		//number of classes
  nel_ta ntest;		//number of test data points
  dim_ta nvar1;

  real_a **train;	//training data vectors
  void *ord;		//training data ordinates
  size_t ordsize;		//size of ordinates
  real_a **test;		//test data vectors
  void *result;		//results of classification
  real_a *con;

  real_a *all;		//for memory allocation

  nel_ta ntrue;
  nel_ta nfalse;
  real_a corr;

  nel_ta nfold;
  nel_ta k;

  char *command;		//the system command
  int commandlen=0;		//maximum length of command

  timeb date;
  pid_t pid;
  int64_t session_id;		//for naming temporary files

  agf_command_opts opt_args;	//option parameters
  //pass to the machine learning commands:
  char normflag[3];
  char normout[7];
  char flagstr[20];		//-n -S svd -u
  char *normfile;		//normalization file name plus option

  int err_code;

  diagfs=stdout;

  //set defaults and parse command line options:
  opt_args.k=-1;
  opt_args.W2=W_DEFAULT_BORDERS;
  opt_args.n=NBORD_DEFAULT;
  opt_args.tol=TOL_DEFAULT;
  opt_args.ftest=F_DEFAULT;
  opt_args.algtype=0;
  
  err_code=agf_parse_command_opts(argc, argv, "a:c:d:f:i:I:k:l:m:N:r:s:S:t:v:V:W:Knu", &opt_args);
  if (err_code==FATAL_COMMAND_OPTION_PARSE_ERROR) return err_code;

  if (argc < 1) {
    printf("\n");
    printf("syntax:  nfold [-n] [-d div] [-c type] [-k k] [-W Wc] [-s n] [-t tol] \n");
    printf("               [control] train [outbase]\n");
    printf("\n");
    printf("arguments:\n");
    printf("  control  = control file for multi-borders training (type=5)\n");
    printf("  train    = binary files containing input training data:\n");
    printf("              .vec for vectors;\n");
    printf("              .cls for classes;\n");
    printf("              .dat for continuum ordinates\n");
    printf("  outbase  = optional binary output files:\n");
    printf("              .cls for class data;\n");
    printf("              .dat for interpolates;\n");
    printf("              .con for confidence ratings;\n");
    printf("              .err for error tolerances\n");
    printf("\n");
    printf("options:\n");
    printf("  -a normf = file containing normalization data (input/output)\n");
    printf("  -c type  = algorithm \n");
    printf("          (0=AGF borders, 1=AGF, 2=KNN, 3=AGF int., 4=KNN int., 5=multi-borders)\n");
    printf("               --default AGF borders\n");
    printf("  -d ndiv  = number of divisions in data (default=%d)\n", (int32_t) (1/opt_args.ftest));
    printf("  -i maxi1 =  maximum number of iterations when searching for class border (%d)\n", (int32_t) agf_global_borders_maxiter);

    printf("  -I maxi2 = maximum number of iterations when calculating weights (%d, %d)\n", (int32_t) agf_global_weights_maxiter, (int32_t) agf_global_borders_maxiter);
    printf("  -k k     = number of nearest neighbours\n");
    printf("  -K       = keep temporary files\n");
    printf("  -l tol   = tolerance of W (default=%g)\n", (float) agf_global_weights_tol);
    printf("  -n       = option to normalise the data\n");
    printf("  -N maxi3 = maximum number of iterations in supernewton (default=%d, %d)\n", (int32_t) agf_global_weights_maxiter, (int32_t) agf_global_borders_maxiter);
    printf("  -r r0    = location of discrimination border (default=0)\n");
    printf("  -s n     = number of times to sample the border (default=%d)\n", (int32_t) opt_args.n);
    printf("  -S nsv   = number of singular values from SVD\n");
    printf("  -t tol   = tolerance of border samples (default=%g)\n", (float) opt_args.tol);
    printf("  -u       = store borders data in un-normalized coordinates\n");
    printf("  -v var1  = lower filter variance bound\n");
    printf("               --default is to use the total variance of the data/n^(2/D)\n");
    printf("  -V var2  = lower filter variance bound\n");
    printf("               --default is to use the total variance of the data\n");
    printf("  -W Wc    = objective total weight (default=%g)\n", (float) opt_args.W2);
    printf("\n");
    return INSUFFICIENT_COMMAND_ARGS;
  }

  //normalization options:
  //default
  strcpy(flagstr, "");
  strcpy(normout, "");

  if (opt_args.svd>0) {
    sprintf(flagstr, "-S %6d ", opt_args.svd);
    strcpy(normout, "-n ");
    commandlen+=10;
  }

  //create string versions of flags:
  if (opt_args.normflag) {
    strcat(flagstr, "-n ");
    strcpy(normout, "-n ");
    commandlen+=3;
  }

  //create string versions of flags:
  if (opt_args.uflag) {
    strcat(flagstr, "-u");
    strcat(normout, "-u");
    commandlen+=3;
  }
  if (opt_args.normfile!=NULL) {
    normfile=new char[strlen(opt_args.normfile)+4];
    sprintf(normfile, "-a %s", opt_args.normfile);
    commandlen+=strlen(normfile)+3;
  } else {
    normfile=new char[1];
    normfile[0]='\0';
  }

  //figure out all the various file names:
  //create a unique session id for naming the temporary files:
  ftime(&date);
  pid=getpid();
  session_id=pid+date.time*100000;        //how big do PID's get??

  if (opt_args.algtype == 5) {
    controlfile=argv[0];
    argc--;
    argv++;
    commandlen+=strlen(controlfile)+1;
  }

  basename=argv[0];
  if (argc>1) {
    tempbase=new char[strlen(argv[1])+21];
    sprintf(tempbase, "%s.%ld", argv[1], session_id);
    outbase=argv[1];
  } else {
    tempbase=new char[strlen(argv[0])+21];
    sprintf(tempbase, "%s.%ld", argv[0], session_id);
    outbase=NULL;
  }

  if (opt_args.algtype >= 3 && opt_args.algtype <= 4) {
    //interpolation/regression algorithm:
    //doesn't have a bundled read routine, so we have to generate the file names:
    char vecname[strlen(basename)+5];
    char ordname[strlen(basename)+5];
    nel_ta n1;

    strcpy(vecname, basename);
    strcat(vecname, ".vec");
    strcpy(ordname, basename);
    strcat(ordname, ".dat");

    train=read_vecfile<real_a>(vecname, nsamp, nvar);
    ord=read_datfile<real_a>(ordname, n1);

    if (nsamp!=n1) {
      fprintf(stderr, "Sample count mismatch: %s, n=%d; %s, n=%d\n", 
		vecname, nsamp, ordname, n1);
      return SAMPLE_COUNT_MISMATCH;
    }

    ordsize=sizeof(real_a);
  } else {
    cls_ta *dum;
    //classification algorithm:
    err_code=agf_read_train(basename, train, dum, nsamp, nvar);
    ord=(cls_ta *) dum;
    if (err_code!=0) return err_code;
    ordsize=sizeof(cls_ta);
  }

  all=train[0];		//must do this because of the way the data is allocated...

  //count the number of classes:
  if (opt_args.algtype<3 || opt_args.algtype>4) {
    nclass=1;
    for (nel_ta i=0; i<nsamp; i++) if (((cls_ta *) ord)[i]>=nclass) nclass++;
  }
  
  //file name for the separated data (training and test):
  trainbase=new char[strlen(tempbase)+5];
  sprintf(trainbase, "%s.trn", tempbase);
  commandlen+=strlen(trainbase)+1;
  
  trainvecfile=new char[strlen(trainbase)+5];
  strcpy(trainvecfile, trainbase);
  strcat(trainvecfile, ".vec");
  
  trainclsfile=new char[strlen(trainbase)+5];
  strcpy(trainclsfile, trainbase);
  
  testbase=new char[strlen(tempbase)+6];
  sprintf(testbase, "%s.test", tempbase);
  
  testfile=new char[strlen(testbase)+5];
  strcpy(testfile, testbase);
  strcat(testfile, ".vec");
  commandlen+=strlen(testfile)+1;
  
  //files containing class borders:
  brdfile=new char[strlen(tempbase)+5];
  sprintf(brdfile, "%s.brd", tempbase);
  commandlen+=strlen(brdfile)+1;
  
  //output data (the stuff we need!):
  outfile=new char[strlen(testbase)+5];
  strcpy(outfile, testbase);
  
  confile=new char[strlen(testbase)+5];
  strcpy(confile, testbase);
  if (opt_args.algtype > 2 && opt_args.algtype <=4) {
    //file extensions for interpolation:
    strcat(trainclsfile, ".dat");
    strcat(outfile, ".dat");
    strcat(confile, ".err");
  } else {
    //file extensions for classification:
    strcat(trainclsfile, ".cls");
    strcat(outfile, ".cls");
    strcat(confile, ".con");
  }

  fprintf(diagfs, "%d training vectors found: %s\n", nsamp, basename);

/*
  for (long i=0; i<nsamp; i++) {
    for (long j=0; j<nvar; j++) printf("%f ", train[i][j]);
    printf("%ld\n", cls[i]);
  }
*/

  //randomize the training data:
  fprintf(diagfs, "Randomizing the data...\n");

  //use a dummy array to randomize the ordinates
  //(which could be real_aing point or longword...)
  cls_ta *ind;
  void *ord1;

  ind=new cls_ta[nsamp];
  ord1=malloc(ordsize*nsamp);

  for (nel_ta i=0; i<nsamp; i++) ind[i]=i;

  randomize_vec(train, nvar, nsamp, ind);

  //this is starting to look a little messy:
  for (nel_ta i=0; i<nsamp; i++) {
    for (size_t j=0; j<ordsize; j++) {
      //((char *) ord1)[ind[i]*ordsize+j]=((char *) ord)[i*ordsize+j];
      ((char *) ord1)[i*ordsize+j]=((char *) ord)[ind[i]*ordsize+j];
    }
  }

  ord=ord1;		//***memory leak***

  //calculate the number of test vectors:
  if (opt_args.div==-1) nfold=1./opt_args.ftest; else nfold=opt_args.div;
  opt_args.ftest=1./(real_a) nfold;

  if (opt_args.algtype != 2 && opt_args.algtype != 4) {
    if (opt_args.k <= opt_args.W2 || opt_args.k >= nsamp) {
      if (opt_args.k != -1) {
        fprintf(stderr, "Parameter k=%d out of range.  Using all the training data.\n", opt_args.k);
        opt_args.k=-1;
        err_code=PARAMETER_OUT_OF_RANGE;
      }
    }
  } else if (opt_args.k <= 0 || opt_args.k >= nsamp) {
    if (opt_args.k!=-1) {
      fprintf(stderr, "Parameter k=%d out of range.  Using k=11.\n", opt_args.k);
      err_code=PARAMETER_OUT_OF_RANGE;
    }
    opt_args.k=11;
  }

  //read in the results of our efforts:
  //result=new cls_ta[nsamp];
  result=malloc(nsamp*ordsize);
  con=new real_a[nsamp];

  k=0;

  //I'm sure this should be enough: (I know, it's a stupid way of doing it...)
  commandlen+=2*strlen(tempbase)+200;
  command=new char[commandlen];
  
  for (nel_ta fold=0; fold<nfold; fold++) {

    //now we need to output the training and test data to separate files:
    fs=fopen(trainvecfile, "w");
    fwrite(&nvar, 1, sizeof(nvar), fs);
    fs2=fopen(trainclsfile, "w");
    fs3=fopen(testfile, "w");
    fwrite(&nvar, 1, sizeof(nvar), fs3);
    //should be able to calculate from first principles,
    //but we're lazy...
    ntest=0;
    for (nel_ta i=0; i<nsamp; i++) {
      if ((nel_ta) (i/opt_args.ftest/nsamp) != fold) {
        fwrite(train[i], sizeof(real_a), nvar, fs);
        fwrite(((char *) ord)+i*ordsize, ordsize, 1, fs2);
      } else {
        ntest++;
        fwrite(train[i], sizeof(real_a), nvar, fs3);
      }

    }
    fclose(fs);
    fclose(fs2);
    fclose(fs3);

    fprintf(diagfs, "Using %d vectors for testing\n", ntest);
  
    //begin the classification scheme:
    fprintf(diagfs, "Beginning classification...\n");

    switch (opt_args.algtype) {
      case (0):
        //now we run the program to find the class borders:
        //generate the command:
        sprintf(command, "time %sclass_borders%s %s %s -r %g -I %ld -l %g -i %ld -v %g -V %g -k %d -W %f -s %d -t %g %s %s", 
		AGF_COMMAND_PREFIX, AGF_OPT_VER, flagstr, normfile, 
		opt_args.rthresh, 
		agf_global_weights_maxiter, agf_global_weights_tol,
		agf_global_borders_maxiter,
		opt_args.var[0], opt_args.var[1], opt_args.k, opt_args.W2, opt_args.n, opt_args.tol, 
		trainbase, brdfile);
        fprintf(diagfs, "%s\n", command);
        err_code=system(command);
        if (err_code!=0) {
          fprintf(stderr, "Command: \"%s\"\n exited with code, %d\n", command, err_code);
          goto fail;
        }
        //now we use the class borders so generated to classify the test data:
        sprintf(command, "%stime classify_b%s %s %s %s %s %s", 
		AGF_COMMAND_PREFIX, AGF_OPT_VER, normout, normfile, brdfile, 
		  testfile, testbase);
        fprintf(diagfs, "%s\n", command);
        err_code=system(command);
        if (err_code!=0) {
          fprintf(stderr, "Command: \"%s\"\n exited with code, %d\n", command, err_code);
          goto fail;
        }
	break;
      case (1):
        sprintf(command, "time %sagf%s %s %s -I %ld -l %g -v %g -V %g -k %d -W %f classify %s %s %s", 
		AGF_COMMAND_PREFIX, AGF_OPT_VER, normfile, flagstr,
		agf_global_weights_maxiter, agf_global_weights_tol,
		opt_args.var[0], opt_args.var[1], opt_args.k, opt_args.W2,  
		trainbase, testfile, testbase);
        fprintf(diagfs, "%s\n", command);
        err_code=system(command);
        if (err_code!=0) {
          fprintf(stderr, "Command: \"%s\"\n exited with code, %d\n", command, err_code);
          goto fail;
        }
        break;
      case (2):
        sprintf(command, "time %sknn%s %s %s -k %d classify %s %s %s", 
		AGF_COMMAND_PREFIX, AGF_OPT_VER, normfile, flagstr,
		opt_args.k, trainbase, testfile, testbase);
        fprintf(diagfs, "%s\n", command);
        err_code=system(command);
        if (err_code!=0) {
          fprintf(stderr, "Command: \"%s\"\n exited with code, %d\n", command, err_code);
          goto fail;
        }
        break;
      case (3):
        sprintf(command, "time %sagf%s %s %s -I %ld -l %g -v %g -V %g -k %d -W %f interp %s %s %s", 
		AGF_COMMAND_PREFIX, AGF_OPT_VER, normfile, flagstr,
		agf_global_weights_maxiter, agf_global_weights_tol,
		opt_args.var[0], opt_args.var[1], opt_args.k, opt_args.W2,  
		trainbase, testfile, testbase);
        fprintf(diagfs, "%s\n", command);
        err_code=system(command);
        if (err_code!=0) {
          fprintf(stderr, "Command: \"%s\"\n exited with code, %d\n", command, err_code);
          goto fail;
        }
        break;
      case (4):
        sprintf(command, "time %sknn%s %s %s -k %d interp %s %s %s", 
		AGF_COMMAND_PREFIX, AGF_OPT_VER, flagstr, normfile, 
		opt_args.k, trainbase, testfile, testbase);
        fprintf(diagfs, "%s\n", command);
        err_code=system(command);
        if (err_code!=0) {
          fprintf(stderr, "Command: \"%s\"\n exited with code, %d\n", command, err_code);
          goto fail;
        }
        break;
      case (5):
        sprintf(command, "time %smulti_borders%s %s %s -r %g -I %ld -l %g -i %ld -v %g -V %g -k %d -W %f -s %d -t %g %s %s %s %s.txt > %s.sh", 
		AGF_COMMAND_PREFIX, AGF_OPT_VER, flagstr, normfile, 
		opt_args.rthresh,
		agf_global_weights_maxiter, agf_global_weights_tol,
		agf_global_borders_maxiter,
		opt_args.var[0], opt_args.var[1], opt_args.k, opt_args.W2, opt_args.n, opt_args.tol, 
		controlfile, trainbase, brdfile, tempbase, tempbase);
        fprintf(diagfs, "%s\n", command);
        err_code=system(command);
        if (err_code!=0) {
          fprintf(stderr, "Command: \"%s\"\n exited with code, %d\n", command, err_code);
          goto fail;
        }
        sprintf(command, "chmod u+x %s.sh", tempbase);
        fprintf(diagfs, "%s\n", command);
        sprintf(command, "bash -e %s.sh", tempbase);
        fprintf(diagfs, "%s\n", command);
        err_code=system(command);
        if (err_code!=0) {
          fprintf(stderr, "Command: \"%s\"\n exited with code, %d\n", command, err_code);
          if (opt_args.Kflag==0) {
            sprintf(command, "rm %s.*.vec", brdfile);
            fprintf(stderr, "%s\n", command);
            system(command);
            sprintf(command, "rm %s.*.cls", brdfile);
            fprintf(stderr, "%s\n", command);
            system(command);
          }
          goto fail;
        }
        //now we use the class borders so generated to classify the test data:
        if (opt_args.normfile==NULL && (opt_args.svd > 0 || opt_args.normflag)) {
          sprintf(command, "time %sclassify_m%s -a %s.std %s %s.txt %s %s", 
  		AGF_COMMAND_PREFIX, AGF_OPT_VER, brdfile, normout, 
		tempbase, testfile, testbase);
        } else {
          sprintf(command, "%stime classify_m%s %s %s %s.txt %s %s", 
  		AGF_COMMAND_PREFIX, AGF_OPT_VER, normfile, normout, 
		tempbase, testfile, testbase);
        }
        fprintf(diagfs, "%s\n", command);
        err_code=system(command);
        if (err_code!=0) {
          fprintf(stderr, "Command: \"%s\"\n exited with code, %d\n", command, err_code);
          goto fail;
        }
        break;
      default:
        fprintf(stderr, "nfold: classification algorithm number %d not recognized", opt_args.algtype);
	exit(COMMAND_OPTION_PARSE_ERROR);
    }

    if (err_code!=0) exit(err_code);

    //read in the results of our efforts:
    fs=fopen(outfile, "r");
    fread((char *) result+k*ordsize, ordsize, ntest, fs);
    fclose(fs);  

    //no confidence rating yet!!  
    fs=fopen(confile, "r");
    fread(con+k, sizeof(real_a), ntest, fs);
    fclose(fs);  

    k+=ntest;

  }

  if (opt_args.algtype > 2 && opt_args.algtype <=4) {
    real_a ave1, ave2;
    real_a var1, var2;
    real_a diff1, diff2;
    real_a cov, r;

    ave1=0;
    ave2=0;
    for (nel_ta i=0; i<nsamp; i++) {
      ave1+=((real_a *) result)[i];
      ave2+=((real_a *) ord)[i];
    }
    ave1/=nsamp;
    ave2/=nsamp;

    var1=0;
    var2=0;
    cov=0;
    for (nel_ta i=0; i<nsamp; i++) {
      diff1=((real_a *) result)[i]-ave1;
      diff2=((real_a *) ord)[i]-ave2;
      var1+=diff1*diff1;
      var2+=diff2*diff2;
      cov+=diff1*diff2;
    }

    r=cov/sqrt(var1*var2);

    printf("Correlation: %g\n", r);

  } else {
    ntrue=0;
    for (nel_ta i=0; i<nsamp; i++) if (((cls_ta *) result)[i]==((cls_ta *) ord)[i]) ntrue++;

    class_eval((cls_ta *) ord, (cls_ta *) result, nsamp);
    check_confidence(((cls_ta *) ord), (cls_ta *) result, con, nsamp, NCONHIST);
  }

  //write the optional output files:
  if (outbase!=NULL) {
    char *outfile1;
    char *confile1;
    void * result1;

    outfile1=new char[strlen(outbase)+10];
    strcpy(outfile1, outbase);
  
    confile1=new char[strlen(outbase)+10];
    strcpy(confile1, outbase);

    result1=malloc(nsamp*ordsize);

    if (opt_args.algtype > 2 && opt_args.algtype <=4) {
      strcat(outfile1, ".dat");
      strcat(confile1, ".err");
      //result1=new cls_ta[nsamp];
      for (nel_ta i=0; i<nsamp; i++) 
		((cls_ta *) result1)[ind[i]]=((cls_ta *) result)[i];
    } else {
      strcat(outfile1, ".cls");
      strcat(confile1, ".con");
      //result1=new real_a[nsamp];
      for (nel_ta i=0; i<nsamp; i++) 
		((real_a *) result1)[ind[i]]=((real_a *) result)[i];
    }

    fs=fopen(outfile1, "w");
    fwrite(result1, ordsize, nsamp, fs);
    fclose(fs);

    real_a con1[nsamp];

    for (nel_ta i=0; i<nsamp; i++) con1[ind[i]]=con[i];

    fs=fopen(confile1, "w");
    fwrite(con1, sizeof(real_a), nsamp, fs);
    fclose(fs);

    delete [] outfile1;
    delete [] confile1;

    free(result1);
  }

  //clean up:


  fail:
    if (opt_args.Kflag) exit(err_code);

    //delete the temporary files so they don't clutter the directory:
    sprintf(command, "rm %s", testfile);
    printf("%s\n", command);
    system(command);
    sprintf(command, "rm %s", outfile);
    printf("%s\n", command);
    system(command);
    sprintf(command, "rm %s", confile);
    printf("%s\n", command);
    system(command);
    sprintf(command, "rm %s", trainvecfile);
    printf("%s\n", command);
    system(command);
    sprintf(command, "rm %s", trainclsfile);
    printf("%s\n", command);
    system(command);
    if (opt_args.algtype==0) {
      sprintf(command, "rm %s.brd", brdfile);
      printf("%s\n", command);
      system(command);
      sprintf(command, "rm %s.bgd", brdfile);
      printf("%s\n", command);
      system(command);
      if ((opt_args.normflag || opt_args.svd>0) && opt_args.normfile==NULL) {
        sprintf(command, "rm %s.std", brdfile);
        system(command);
      }
    } else if (opt_args.algtype==5) {
      sprintf(command, "rm %s*.brd", brdfile);
      printf("%s\n", command);
      system(command);
      sprintf(command, "rm %s*.bgd", brdfile);
      printf("%s\n", command);
      system(command);
      sprintf(command, "rm %s.txt", tempbase);
      printf("%s\n", command);
      system(command);
      sprintf(command, "rm %s.sh", tempbase);
      printf("%s\n", command);
      system(command);
      if ((opt_args.normflag || opt_args.svd>0) && opt_args.normfile==NULL) {
        sprintf(command, "rm %s.std", brdfile);
        system(command);
      }
    } else {
      if ((opt_args.normflag || opt_args.svd>0) && opt_args.normfile==NULL) {
        sprintf(command, "rm %s.std", testbase);
        system(command);
      }
    }

    //delete [] result;
    delete [] all;
    delete [] train;

    delete [] command;
    delete [] normfile;
    delete [] tempbase;
  
    //delete various file names:
    delete [] trainbase;
    delete [] trainvecfile;
    delete [] trainclsfile;
  
    delete [] testbase;
    delete [] testfile;
    delete [] outfile;
    delete [] confile;
  
    delete [] brdfile;
  
    delete [] ind;
  //  delete [] clind;
    //delete [] acc_mat[0];
    //delete [] acc_mat;

    free(ord1);
    free(result);
    delete [] con;

  return err_code;

}


