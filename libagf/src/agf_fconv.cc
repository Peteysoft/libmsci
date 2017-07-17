//
// Input/output routines for libAGF: chiefly binary dumps of vectors and
// matrices.
//

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
#include "agf_fconv.h"

using namespace std;
using namespace libpetey;

namespace libagf {

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

template <>
void get_format_code<int32_t>(char *code) {
  code[0]='d';
  code[1]='\0';
}

template <>
void get_format_code<int64_t>(char *code) {
  code[0]='l';
  code[1]='d';
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

  //scans a single line for features data in LIBSVM format:
  template <class real>
  dim_ta scan_svm_features(char *line, dim_ta *&ind2, real *&raw) {
    int err;			//number of items scanned
    int pos;			//position in the line
    int len;			//length of the line
    int rel;			//relative position in remaining fraction of line
    dim_ta ind;			//for reading in the dimension index
    real val;			//for reading in the feature data
    linked_list<real> feat;	//linked list of features data
    linked_list<dim_ta> sub;	//linked list of indices
    long cnt1;			//number of indices collected
    long cnt2;			//number of features collected
    char format[10];		//format code
    char fcode[4];		//format code for features data

    //not particularly efficient for something that occurs every line:
    get_format_code<real>(fcode);
    sprintf(format, "%%%s%%d%%n", fcode);

    len=strlen(line);
    err=sscanf(line, "%d%n", &ind, &pos);
    if (err!=1) {
      fprintf(stderr, "scan_svm_features: error; missing features data\n");
      return 0;
      //exit(FILE_READ_ERROR);
    } 

    while (pos<len) {
      sub.add(ind);
      for (pos=pos; line[pos]!=':' && pos<len; pos++);
      if (pos>=len) {
        fprintf(stderr, "scan_svm_features: error; feature %d missing\n", ind);
	return 0;
        //exit(FILE_READ_ERROR);
      } 
      pos++;
      err=sscanf(line+pos, format, &val, &ind, &rel);
      if (err<=0) {
        fprintf(stderr, "read_svm: error; feature %d missing value\n", ind);
	return 0;
        //exit(FILE_READ_ERROR);
      }
      feat.add(val);
      if (err < 2) break;
      //printf("2: nscan=%d\n", err);
      pos+=rel;
    }
    ind2=sub.make_array(cnt1);
    raw=feat.make_array(cnt2);
    assert(cnt1==cnt2);		//should always be true?

    return cnt1;
  }


  //should work for integer as well as floating-point ordinates:
  template <class real, class cls_t>
  nel_ta read_svm(FILE *ifs, real **&train, cls_t *&cls, dim_ta &nvar, real missing) {

    long n;
    real **vec;			//vectors as read in (may be too short)
    dim_ta *nfeat;		//total number of features in each vector in vec
    int err=0;

    char **line;		//raw ascii lines read from file
    int pos;			//position in line

    dim_ta **ind;		//indices of features
    real **raw;			//raw features data

    line=read_ascii_all(ifs, &n);

    cls=new cls_t[n];
    nfeat=new dim_ta[n];

    ind=new dim_ta*[n];
    raw=new real*[n];

    char fcode0[4];
    char format0[10];

    get_format_code<cls_t>(fcode0);
    sprintf(format0, "%%%s%%n", fcode0);

    nvar=0;

    //read in the raw data (no preset maximums, lines can have any number of features):
    for (nel_ta i=0; i<n; i++) {
      long cnt1, cnt2;		//# of indices/ features read in
      //err=sscanf(line[i], "%d%d%n", cls+i, &ind, &pos);
      cnt1=sscanf(line[i], format0, cls+i, &pos);
      if (cnt1==0) {
        n=i;
        break;
      }
      nfeat[i]=scan_svm_features(line[i]+pos, ind[i], raw[i]);
      if (nfeat[i]<=0) {
        fprintf(stderr, "read_svm: error on line %d; skipping\n", i);
        err=FILE_READ_ERROR;
        //probably not good form but I don't feel like rewriting the control structure:
        n--;
        i--;
        continue;
        //exit(FILE_READ_ERROR);
      } 
      for (dim_ta j=0; j<nfeat[i]; j++) if (ind[i][j]>nvar) nvar=ind[i][j];
      delete [] line[i];
    }
    delete [] line;

    //printf("features data: %d X %d\n", n, nvar);

    vec=new real *[n];
    vec[0]=new real[n*nvar];

    //fill in missing values:
    for (nel_ta i=0; i<n; i++) {
      //printf("line %d; %d features ...", i, nfeat[i]);
      vec[i]=vec[0]+i*nvar;
      for (dim_ta j=0; j<nvar; j++) vec[i][j]=missing;
      for (dim_ta j=0; j<nfeat[i]; j++) vec[i][ind[i][j]-1]=raw[i][j];
      //printf("%g ", vec[i][j]);
      delete [] ind[i];
      delete [] raw[i];
      //printf("%d\n", cls[i]);
    }
    delete [] ind;
    delete [] raw;
    delete [] nfeat;

    train=vec;
    if (err!=0) n=-n;
    return n;
  }

template <class real, class cls_t>
nel_ta read_svm(FILE *fs, real **&train, cls_t *&cls, dim_ta &nvar, real missing, int Uflag) {
  cls_t ncls=0;			//number of classes
  nel_ta result;
  //real *ord;

  result=read_svm(fs, train, cls, nvar, missing);

  //cls=new cls_t[result];
  for (nel_ta i=0; i<result; i++) {
    //cls[i]=ord[i];
    if (cls[i]>=ncls) ncls=cls[i]+1;
  }
  //delete [] ord;

  //convert class labels to go from 0..nc-1:
  if (Uflag) compress_labels(cls, result);

  return result;
}

template <class real, class cls_t>
nel_ta read_svm(const char *fname, real **&train, cls_t *&cls, dim_ta &nvar, real missing, int Uflag) {
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

  linked_list<cls_t> label;    //list of labels
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
    cls_t clab;

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
      cls=new cls_t[n];
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
        for (cls_t j=0; j<ncls; j++) p[i][j]=0;
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

template nel_ta read_svm<double, double>(FILE *fs, double **&, double *&, dim_ta &, double);
template nel_ta read_svm<float, float>(FILE *fs, float **&, float *&, dim_ta &, float);
template nel_ta read_svm<float, cls_ta>(FILE *fs, float **&, cls_ta *&, dim_ta &, float);
template nel_ta read_svm<double, cls_ta>(FILE *fs, double **&, cls_ta *&, dim_ta &, double);
template nel_ta read_svm<float, cls_ta>(FILE *fs, float **&, cls_ta *&, dim_ta&, float, int);
template nel_ta read_svm<double, cls_ta>(FILE *fs, double **&, cls_ta *&, dim_ta&, double, int);
template nel_ta read_svm<float, cls_ta>(const char *, float **&, cls_ta *&, dim_ta &, float, int);
template nel_ta read_svm<double, cls_ta>(const char *, double **&, cls_ta *&, dim_ta &, double, int);

template void print_lvq_svm<cls_ta, float>(FILE *fs, float **, cls_ta *, nel_ta, dim_ta, int, int);
template void print_lvq_svm<cls_ta, double>(FILE *fs, double **, cls_ta *, nel_ta, dim_ta, int, int);
template nel_ta read_svmout<cls_ta, float>(FILE *ifs, cls_ta *&, float **&, cls_ta &, nel_ta);
template nel_ta read_svmout<cls_ta, double>(FILE *ifs, cls_ta *&, double **&, cls_ta &, nel_ta);


} //end namespace libagf

