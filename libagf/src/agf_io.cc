
//Copyright (C) 2007 Peter Mills.  All rights reserved.

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>

#include "error_codes.h"

#include "peteys_tmpl_lib.h"
#include "read_ascii_all.h"
#include "full_util.h"
#include "linked.h"

#include "agf_util.h"
#include "agf_io.h"

using namespace std;
using namespace libpetey;

namespace libagf {

#define STATS_HEADER "dim    average  std. dev.\n"

//-returns null pointer on failure
//-dimension is -1 if this is a bad value
//-sample number is -1 if there is an allocation failure
real_a ** read_vecfile(const char *filename, nel_ta &m, dim_ta &n) {
  FILE *fs;
  real_a **data;
  nel_ta n1;			//dimensions have to be same type for matrix routines

  fs=fopen(filename, "r");
  if (fs == NULL) {
    m=0; n=0;
    return NULL;
  }

  //just piggy-back off the matrix routines:
  data=read_matrix<real_a, nel_ta>(fs, m, n1);
  n=n1;

  fclose(fs);

  return data;
}

cls_ta * read_clsfile(const char *filename, nel_ta &n) {
  FILE *fs;
  cls_ta *data;

  fs=fopen(filename, "r");
  if (fs == NULL) {
    n=0;
    return NULL;
  }

  fseek(fs, 0, SEEK_END);
  n=ftell(fs);
  if (n % sizeof(cls_ta) != 0) { //check for consistency
    n=-1;
    fclose(fs);
    return NULL;
  }
  n=n/sizeof(cls_ta);
  fseek(fs, 0, SEEK_SET);

  data=new cls_ta [n];

  fread(data, sizeof(cls_ta), n, fs);

  fclose(fs);

  return data;
}

real_a * read_datfile(const char *filename, nel_ta &n) {
  FILE *fs;
  real_a *data;

  fs=fopen(filename, "r");
  if (fs == NULL) {
    n=0;
    return NULL;
  }

  fseek(fs, 0, SEEK_END);
  n=ftell(fs);
  if (n % sizeof(real_a) != 0) { //check for consistency
    n=-1;
    fclose(fs);
    return NULL;
  }
    
  n=n/sizeof(real_a);
  fseek(fs, 0, SEEK_SET);

  data=new real_a [n];
  if (data==NULL) {
    //doesn't c++ throw bad allocs anyway??
    n=-1;
    fclose(fs);
    return data;
  }

  fread(data, sizeof(real_a), n, fs);

  fclose(fs);

  return data;
}

#define MAXLL 200

int read_stats(const char *filename, real_a *ave, real_a *std, dim_ta ndim) {
  FILE *fs;
  dim_ta ivar;
  char header[MAXLL];

  fs=NULL;
  fs=fopen(filename, "r");
  if (fs==NULL) return UNABLE_TO_OPEN_FILE_FOR_READING;

  fgets(header, MAXLL, fs);
  for (dim_ta i=0; i<ndim; i++) {
    fscanf(fs, "%d %g %g", &ivar, &ave[i], &std[i]);
  }

  fclose(fs);

  return 0;
}

real_a ** read_stats2(const char *filename, real_a *&ave, dim_ta &m, dim_ta &n) {
  FILE *fs;
  real_a **mat;
  dim_ta ivar;
  char *header=NULL;
  char **line=NULL;
  long nline;
  int ncon;

  ave=NULL;

  fs=fopen(filename, "r");
  if (fs==NULL) return NULL;

  header=fget_line(fs, 1);
  if (header==NULL) goto fail;
  if (strcmp(header, STATS_HEADER)==0) {
    line=read_ascii_all(fs, &nline, 1);
    if (line==NULL || nline<=0) goto fail;
    n=nline;
    ave=new real_a[n];
    mat=allocate_matrix<real_a, nel_ta>(n, n);
    for (dim_ta i=0; i<n; i++) {
      if (strlen(line[i])==0) {
        n=i+1;
        break;
      }
      ncon=sscanf(line[i], "%d %g %g", &ivar, ave+i, mat[i]+i);
      if (ncon!=3) {
        goto fail;
      }
      delete [] line[i];
    }
    m=n;
    delete [] line;
  } else if (header==NULL) {
    goto fail;
  } else {
    fseek(fs, 0, SEEK_SET);
    mat=read_matrix<real_a, nel_ta>(fs, m, n);
    if (mat==NULL || m<=0 || n<=0) goto fail;
    ave=new real_a[m];
    //printf("read_stats2: constant term:\n");
    n--;
    for (dim_ta i=0; i<m; i++) {
      ave[i]=mat[i][n];
      //fprintf(stderr, "%g\n", ave[i]);
    }
    //printf("read_stats2: transformation matrix:\n");
    //print_matrix(stdout, mat, m, n);
  }

  delete [] header;

  fclose(fs);

  return mat;

  fail:			//clean up in event of read failure and set error indicators
    fprintf(stderr, "read_stats2: an error occurred reading data file, %s\n", filename);
    if (line!=NULL) {
      for (long i=0; i<nline; i++) delete [] line[i];
      delete [] line;
    }
    if (ave!=NULL) delete [] ave;
    if (mat!=NULL) delete_matrix(mat);
    if (fs!=NULL) fclose(fs);
    if (header!=NULL) delete [] header;
    n=-1;
    m=-1;
  return NULL;
}

int print_stats(FILE *fs, real_a *ave, real_a *std, dim_ta ndim) {
  fprintf(fs, STATS_HEADER);
  for (dim_ta i=0; i<ndim; i++) {
    fprintf(fs, "%3d %10.6g %10.6g\n", i, ave[i], std[i]);
  }

  return 0;
}

int agf_read_train(const char *fbase, real_a **&train, cls_ta *&cls, nel_ta &n, dim_ta &nvar) {
  char *vecfile=NULL;
  char *classfile=NULL;
  nel_ta n1;
  int err=0;			//return code

  train=NULL;
  cls=NULL;

  vecfile=new char[strlen(fbase)+5];
  sprintf(vecfile, "%s.vec", fbase);

  classfile=new char[strlen(fbase)+5];
  sprintf(classfile, "%s.cls", fbase);

  //read in the training data:
  train=read_vecfile(vecfile, n, nvar);
  if (nvar <= 0 || n <= 0) {
    fprintf(stderr, "Error reading file: %s\n", vecfile);
    err=FILE_READ_ERROR;
    goto fail;
  }
  if (train==NULL) {
    fprintf(stderr, "Unable to open file, %s, for reading.\n", vecfile);
    err=UNABLE_TO_OPEN_FILE_FOR_READING;
    goto fail;
  }

  cls=read_clsfile(classfile, n1);

  if (n1 <= 0) {
    fprintf(stderr, "Error reading file: %s\n", classfile);
    err=FILE_READ_ERROR;
    goto fail;
  }
  if (cls == NULL) {
    fprintf(stderr, "Unable to open file, %s, for reading.\n", classfile);
    err=UNABLE_TO_OPEN_FILE_FOR_READING;
    goto fail;
  }
  if (n1!=n) {
    fprintf(stderr, "Sample count mismatch: %d in %s, %d in %s.\n", n, vecfile, n1, classfile);
    err=SAMPLE_COUNT_MISMATCH;
    goto fail;
  }

  fail:

    delete [] classfile;
    delete [] vecfile;

  return err;
}	     

int agf_read_borders(const char *fbase, real_a **&brd, real_a **&grd, nel_ta &n, dim_ta &nvar) {
  char *brdfile;
  char *grdfile;
  nel_ta n1;
  dim_ta nvar1;
  int err=0;

  brd=NULL;
  grd=NULL;

  brdfile=new char[strlen(fbase)+5];
  strcpy(brdfile, fbase);
  strcat(brdfile, ".brd");

  grdfile=new char[strlen(fbase)+5];
  strcpy(grdfile, fbase);
  strcat(grdfile, ".bgd");

  //read in the decision surface data:
  brd=read_vecfile(brdfile, n, nvar);
  if (nvar <= 0 || n <= 0) {
    fprintf(stderr, "Error reading file: %s\n", brdfile);
    err=FILE_READ_ERROR;
    goto fail;
  }
  if (brd == NULL) {
    fprintf(stderr, "Unable to open input file: %s\n", brdfile);
    err=UNABLE_TO_OPEN_FILE_FOR_READING;
    goto fail;
  }

  grd=read_vecfile(grdfile, n1, nvar1);

  if (nvar1 <= 0 || n1 <= 0) {
    fprintf(stderr, "Error reading file: %s\n", grdfile);
    err=FILE_READ_ERROR;
    goto fail;
  }
  if (grd == NULL) {
    fprintf(stderr, "Unable to open input file: %s\n", grdfile);
    err=UNABLE_TO_OPEN_FILE_FOR_READING;
    goto fail;
  }
  if (nvar1 != nvar) {
    fprintf(stderr, "Error: dimensions of border and gradient vectors do not agree:\n");
    fprintf(stderr, "       %s: D=%d vs. %s: D=%d\n", brdfile, nvar, grdfile, nvar1);
    err=DIMENSION_MISMATCH;
    goto fail;
  }
  if (n1 != n) {
    fprintf(stderr, "Error: number of samples in border and gradient files do not agree:\n");
    fprintf(stderr, "       %d in %s vs. %d in %s\n", n, brdfile, n1, grdfile);
    err=SAMPLE_COUNT_MISMATCH;
    goto fail;
  }

  fail:
  
    delete [] grdfile;
    delete [] brdfile;

  return err;
}

char *compile_precommand(const char *fname, agf_command_opts *optargs) {
  char *command=NULL;

  if (optargs->normflag || optargs->svd>0 || optargs->normfile!=NULL) {
    //if the user wants some pre-processing we farm it out to "agf_precondition"

    //if (optargs.normfile == NULL) {
    //  optargs.normfile=new char[strlen(argv[3])+5];
    //  sprintf(optargs.normfile, "%s.std", argv[3]);
    //}
    command=new char[strlen(optargs->normfile)+strlen(fname)+100];
    sprintf(command, "%s%s%s -a %s", AGF_COMMAND_PREFIX, 
		AGF_LTRAN_COM, AGF_OPT_VER, optargs->normfile);
    if (optargs->normflag) strcat(command, " -n");
    if (optargs->asciiflag) strcat(command, " -A");
    if (optargs->Mflag) strcat(command, " -M");
    if (optargs->svd>0) {
      sprintf(command+strlen(command), " -S %d", optargs->svd);
    }
    sprintf(command+strlen(command), " %s", fname);
  }
  
  return command;
  
}

//increasingly I'm finding this bit really brain-dead:
real_a **agf_get_features(const char *fbase, agf_command_opts *opt_args, dim_ta &nvar, nel_ta &n, flag_a sufflag)
{
  real_a **train;
  char *vecfile;

  vecfile=new char[strlen(fbase)+5];
  if (sufflag) {
    sprintf(vecfile, "%s", fbase);
  } else {
    sprintf(vecfile, "%s.vec", fbase);
  }

  if (opt_args->normflag || opt_args->svd>0 || opt_args->normfile != NULL) {
    //if the user wants some pre-processing we farm it out to "agf_precondition"
    FILE *fs;
    char *command;
    nel_ta nvar1;

    //if (opt_args.normfile == NULL) {
    //  opt_args.normfile=new char[strlen(argv[3])+5];
    //  sprintf(opt_args.normfile, "%s.std", argv[3]);
    //}
    command=new char[strlen(opt_args->normfile)+strlen(fbase)+50];
    sprintf(command, "%s%s%s -a %s", AGF_COMMAND_PREFIX, 
		AGF_LTRAN_COM, AGF_OPT_VER, opt_args->normfile);
    if (opt_args->normflag) strcat(command, " -n");
    //if (opt_args->asciiflag) strcat(command, " -A");
    //if (opt_args->Mflag) strcat(command, " -M");
    if (opt_args->svd>0) {
      sprintf(command+strlen(command), " -S %d", opt_args->svd);
    }
    sprintf(command+strlen(command), " %s", vecfile);
    fprintf(stderr, "%s\n", command);
    fs=popen(command, "r");
    train=read_matrix<real_a, nel_ta>(fs, n, nvar1);
    nvar=nvar1;
    pclose(fs);
    delete [] command;
  } else {
    //read in the training data:
    train=read_vecfile(vecfile, n, nvar);
  }
  if (nvar <= 0 || n <= 0) {
    fprintf(stderr, "Error reading file: %s\n", vecfile);
    exit(FILE_READ_ERROR);
  }
  if (train==NULL) {
    fprintf(stderr, "Unable to open file, %s, for reading.\n", vecfile);
    exit(UNABLE_TO_OPEN_FILE_FOR_READING);
  }

  delete [] vecfile;
  return train;

}

template <>
void get_format_code<float>(char *code) {
  code[0]='g';
  code[1]='\0';
}

template <>
void get_format_code<double>(char *code) {
  code[0]='l';
  code[1]='g';
  code[2]='\0';
}

//reads ASCII files in the same format as used by Kohonen's LVQ package
//returns number of lines read in
//if there is an error, returns a negative number
//flag: 1st bit=no header; 2nd bit=no class data; 3rd bit: omit class data
nel_ta read_lvq(FILE *fs, real_a **&train, cls_ta *&cls, dim_ta &nvar, int flags) {

  char **line;			//all the lines in the file
  int nread;			//for counting the number of characters read in a line
  int pos;			//string pointer
  int err;			//counts number of items read
  long n;			//number of lines
  real_a dum;			//for counting number of features
  nel_ta count;			//number of lines actually read
  int hflag, cflag, oflag;	//no header, omit class data, no class data
  int hflag2;			//not hflag
  linked_list<real_a> first_line;		//first line data

  char fcode[4];
  char format1[10];

  //separate out each of the flags:
  hflag=flags & 1;		//no header
  cflag=(flags & 2) >> 1;	//there is no class data
  oflag=(flags & 4) >> 2;	//omit class data

  //printf("hflag=%d; cflag=%d; oflag=%d\n", hflag, cflag, oflag);

  //determine the format codes:
  get_format_code<real_a>(fcode);
  sprintf(format1, "%%%s%%n", fcode);

  train=NULL;
  cls=NULL;

  if (hflag) hflag2=0; else hflag2=1;

  line=read_ascii_all(fs, &n, 1);
  //won't use explicit error messages since we are reading a text file:
  //users can examine it for themselves
  //however we will return the line number where the error occurred...
  if (line==NULL || n<=1) return -1;

  if (hflag) {
    pos=0;
    nvar=0;
    while (sscanf(line[0]+pos, format1, &dum, &nread)==1) {
      first_line.add(dum);
      pos+=nread;
      nvar++;
    }
    if (cflag==0) nvar--;
    //printf("read_lvq: found %d features\n", nvar);
  } else {
    err=sscanf(line[0], "%d", &nvar);
    if (err!=1) return -1;
    n--;
  }

  train=new real_a*[n];
  train[0]=new real_a[n*nvar];
  if (cflag==0 && oflag==0) {
    cls=new cls_ta[n];
  }

  //so we don't have to re-read the first line:
  if (hflag) {
    for (dim_ta i=0; i<nvar; i++) {
      first_line.pop(dum);
      train[0][i]=dum;
    }
    if (cflag==0 && oflag==0) {
      first_line.pop(dum);
      cls[0]=(cls_ta) dum;
    }
  }

  count=hflag;
  for (nel_ta i=hflag; i<n; i++) {
    if (strlen(line[i+hflag2])==0) continue;
    train[count]=train[0]+count*nvar;
    pos=0;
    for (dim_ta j=0; j<nvar; j++) {
      err=sscanf(line[i+hflag2]+pos, format1, train[count]+j, &nread);
      if (err != 1) {
        if (j>0) count=-count-1-hflag;	//not an error if blank line
        break;
      }
      pos+=nread;
    }
    if (err != 1) break;
    if (cflag || oflag) {
      count++;
      continue;
    }
    err=sscanf(line[i+hflag2]+pos, "%d", cls+count);
    if (err != 1) {
      count=-count-1-hflag;
      break;
    }
    count++;
  }

  return count;

}

nel_ta read_lvq(const char *fname, real_a **&train, cls_ta *&cls, dim_ta &nvar, int flags) {
  nel_ta result;

  FILE *fs;
  fs=fopen(fname, "r");
  if (fs==NULL) return -1;
  result=read_lvq(fs, train, cls, nvar, flags);
  fclose(fs);

  return result;
}

//pretty stupid, but kind of fun to write:
cls_ta * read_lvq_classes(FILE *fs, nel_ta &n, int hflag) {
  int c1;
  char c2;
  char cold;
  int tranloc;			//last position where there is a transition
				//between space and number
  int count;			//position in line
  linked_list<cls_ta> data;
  cls_ta val;
  long n1;
  cls_ta *cls;

  if (hflag==0) {
    //get header, don't care what's in it...
    do {
      c1=fgetc(fs);
      if (c1==EOF) break;
      c2=(char) c1;
    } while (c2!='\n');
  }

  while (feof(fs)==0) {
    c1=fgetc(fs);
    if (c1==EOF) break;
    cold=(char) c1;
    count=1;
    if (isdigit(cold)) tranloc=0; else tranloc=-1;
    do {
      c1=fgetc(fs);
      if (c1==EOF) break;
      c2=(char) c1;
      if (isspace(cold) and isdigit(c2)) tranloc=count;
      count++;
    } while (c2!='\n');
    if (tranloc<0) break;
    fseek(fs, tranloc-count, SEEK_CUR);
    fscanf(fs, "%d", &val);
    data.add(val);
  }
      
  cls=data.make_array(n1);
  n=n1;
  return cls;

}

nel_ta read_svm(FILE *ifs, real_a **&train, cls_ta *&cls, dim_ta &nvar, real_a missing, int Uflag) {

  long n;
  real_a **vec;			//vectors as read in (may be too short)
  real_a **vec2;		//final vectors
  int *cls_lst;			//list of classes
  cls_ta ncls;			//number of classes
  dim_ta *nfeat;		//total number of features in each vector in vec
  dim_ta ind;
  int err;

  char **line;			//raw ascii lines read from file
  int pos;			//position in line
  int rel;			//relative position in line

  linked_list<real_a> feat;	//raw features read in 
  linked_list<dim_ta> sub;	//raw indices

  line=read_ascii_all(ifs, &n);

  vec=new real_a *[n];
  cls=new cls_ta[n];
  nfeat=new dim_ta[n];

  char fcode[4];
  char format1[10];

  get_format_code<real_a>(fcode);
  sprintf(format1, "%%%s%%d%%n", fcode);

  nvar=0;
  ncls=0;

  //read in the raw data (no preset maximums, lines can have any number of features):
  for (nel_ta i=0; i<n; i++) {
    real_a val;			//value of feature
    real_a *raw;		//array of raw features data
    dim_ta *ind2;		//array of feature indices
    long cnt1, cnt2;		//# of indices/ features read in
    int len=strlen(line[i]);
    err=sscanf(line[i], "%d%d%n", cls+i, &ind, &pos);
    if (err==0) {
      n=i;
      break;
    }
    if (cls[i]>=ncls) ncls=cls[i]+1;
    if (err==1) {
      fprintf(stderr, "read_svm: error; missing features data on line %d\n", i);
      exit(FILE_READ_ERROR);
    } 

    nfeat[i]=0;

    while (pos<len) {
      sub.add(ind);
      if (ind>nfeat[i]) {
        nfeat[i]=ind;
      }
      for (pos=pos; line[i][pos]!=':' && pos<len; pos++);
      if (pos>=len) {
        fprintf(stderr, "read_svm: error; feature %d missing value on line %d\n", ind, i);
        exit(FILE_READ_ERROR);
      } 
      pos++;
      err=sscanf(line[i]+pos, format1, &val, &ind, &rel);
      if (err<=0) {
        fprintf(stderr, "read_svm: error; feature %d missing value on line %d\n", ind, i);
        exit(FILE_READ_ERROR);
      }
      feat.add(val);
      if (err < 2) break;
      //printf("2: nscan=%d\n", err);
      pos+=rel;
    }
    ind2=sub.make_array(cnt1);
    raw=feat.make_array(cnt2);
    assert(cnt1==cnt2);
    vec[i]=new real_a[nfeat[i]];
    for (dim_ta j=0; j<nfeat[i]-1; j++) vec[i][j]=missing;
    //printf("line %d; %d features ...", i, nfeat[i]);
    for (dim_ta j=0; j<cnt1; j++) {
      vec[i][ind2[j]-1]=raw[j];
      //printf(" %d:%g", ind2[j], raw[j]);
    }
    //printf("\n");
    if (nfeat[i]>nvar) nvar=nfeat[i];

    delete [] raw;
    delete [] ind2;
    delete [] line[i];
    sub.reset();
    feat.reset();
  }
  delete [] line;

  //printf("features data: %d X %d\n", n, nvar);

  vec2=new real_a *[n];
  vec2[0]=new real_a[n*nvar];
  cls_lst=new cls_ta[ncls];
  for (cls_ta j=0; j<ncls; j++) cls_lst[j]=0;

  //fill in missing values:
  for (nel_ta i=0; i<n; i++) {
    //printf("line %d; %d features ...", i, nfeat[i]);
    cls_lst[cls[i]]=1;
    vec2[i]=vec2[0]+i*nvar;
    for (dim_ta j=0; j<nvar; j++) {
      if (j<nfeat[i]) {
        vec2[i][j]=vec[i][j];
      } else {
        vec2[i][j]=missing;
      }
      //printf("%g ", vec2[i][j]);
    }
    delete [] vec[i];
    //printf("%d\n", cls[i]);
  }
  delete [] vec;

  //convert class labels to go from 0..nc-1:
  if (Uflag) {
    cls_ta ncls2=0;
    for (cls_ta j=0; j<ncls; j++) {
      if (cls_lst[j]!=0) {
        cls_lst[j]=ncls2;
        ncls2++;
      }
    }
    for (nel_ta i=0; i<n; i++) cls[i]=cls_lst[cls[i]];
  }

  delete [] cls_lst;
  delete [] nfeat;
  train=vec2;

  return n;

}

nel_ta read_svm(const char *fname, real_a **&train, cls_ta *&cls, dim_ta &nvar, real_a missing, int Uflag) {
  nel_ta result;

  FILE *fs;
  fs=fopen(fname, "r");
  if (fs==NULL) return -1;
  result=read_svm(fs, train, cls, nvar, missing, Uflag);
  fclose(fs);

  return result;
}

template <class cls_t, class real>
void print_lvq_svm(FILE *fs, real **vec, cls_t *cls, nel_ta n, dim_ta nvar, int svmflag, int nhflag) {
  char fcode[4];
  char format[10];
  real dum;

  get_format_code<real>(fcode);
  if (svmflag) {
    sprintf(format, " %%d:%%%s", fcode);
    for (nel_ta i=0; i<n; i++) {
      fprintf(fs, "%d", cls[i]);
      for (dim_ta j=0; j<nvar; j++) fprintf(fs, format, j+1, vec[i][j]);
      fprintf(fs, "\n");
    }
  } else {
    sprintf(format, "%%%s ", fcode);
    if (nhflag==0) fprintf(fs, "%d\n", nvar);
    for (nel_ta i=0; i<n; i++) {
      for (dim_ta j=0; j<nvar; j++) fprintf(fs, format, vec[i][j]);
      if (cls==NULL) fprintf(fs, "\n"); else fprintf(fs, "%d\n", cls[i]);
    }
  }
}

template <class cls_t, class real>
nel_ta read_svmout(FILE *ifs, cls_t *&cls, real **&p, cls_t &ncls, nel_ta n) {
  char **line;
  long nline;                   //number of lines

  linked_list<cls_ta> label;    //list of labels
  cls_t *label_list;           //list of labels

  long ncls1;                   //number of class labels
  cls_t ncls2;                 //max label+1
  real pdf1;			//cond. prob. read in
  int err;

  char format1[10];
  char fcode[4];

  int strptr;
  int32_t stradv;

  //determine the format codes:
  get_format_code<real>(fcode);
  sprintf(format1, "%%%s%%n", fcode);

  line=read_ascii_all(ifs, &nline);

  for (strptr=0; line[0][strptr]!='\0' && isblank(line[0][strptr]); strptr++);

  if (line[0][strptr]=='l') {
    cls_ta clab;

    //read in the labels:
    strptr+=6;
    err=1;
    ncls2=0;
    for (;;) {
      err=sscanf(line[0]+strptr, "%d%n", &clab, &stradv);
      //printf("%d %d\n", clab, err);
      if (err != 1) break;
      if (clab >= ncls2) ncls2=clab+1;
      label.add(clab);
      strptr+=stradv;
    }

    //index the labels:
    label_list=label.make_array(ncls1);

    /*
    printf("class labels:\n");
    for (int i=0; i<ncls1; i++) printf("%d ", label_list[i]);
    printf("\n");
    */

    //only allocate the arrays if n is not set:
    if (n<=0) {
      n=nline-1;
      cls=new cls_ta[n];
      p=new real *[n];
      p[0]=new real[n*ncls2];
    }

    delete [] line[0];
    for (nel_ta i=0; i<n; i++) {
      err=sscanf(line[i+1], "%d%n", cls+i, &strptr);
      if (err != 1) {
        n=i;
        break;
      }
      for (long j=0; j<ncls2; j++) {
        err=sscanf(line[i+1]+strptr, format1, &pdf1, &stradv);
        if (err != 1) {
          fprintf(stderr, "svmout2agf: read error--missing data element %d on line %d\n", j, i+2);
          exit(FILE_READ_ERROR);
        }
        strptr+=stradv;
        p[i][label_list[j]]=pdf1;
      }

      //printf("%d %f\n", cls, con);
      //for (cls_ta j=0; j<ncls2; j++) printf("%10.6f ", pdf2[j]);
      //printf("%5d\n", cls[i]);
      delete [] line [i+1];
    }

    for (nel_ta i=n+1; i<nline; i++) delete [] line[i];
    delete [] label_list;
    ncls=ncls2;
  } else if (isdigit(line[0][strptr])) {
    ncls=1;
    if (n<=0) {
      cls=new cls_t[nline];
      p=NULL;
    }
    n=nline;
    for (nel_ta i=0; i<n; i++) {
      err=sscanf(line[i], "%d", cls+i);
      if (err!=1) {
        n=i;
        break;
      }
      if (cls[i]>=ncls) ncls=cls[i]+1;
      delete [] line[i];
    }
    for (nel_ta i=n; i<nline; i++) delete [] line[i];
    if (p!=NULL) {
      for (nel_ta i=0; i<n; i++) {
        for (cls_ta j=0; j<ncls; j++) p[i][j]=0;
        if (cls[i]>=0 && cls[i]<=ncls-1) p[i][cls[i]]=1;
      }
    }
  } else {
    fprintf(stderr, "svmout2agf: format error on line 1\n");
    exit(FILE_READ_ERROR);
  }

  delete [] line;

  return n;
}


template void print_lvq_svm<cls_ta, real_a>(FILE *fs, real_a **, cls_ta *, nel_ta, dim_ta, int, int);
template nel_ta read_svmout<cls_ta, real_a>(FILE *ifs, cls_ta *&, real_a **&, cls_ta &, nel_ta);


} //end namespace libagf

