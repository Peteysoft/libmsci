#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "error_codes.h"
#include "linked.h"
#include "read_ascii_all.h"
#include "agf_defs.h"

using namespace std;
using namespace libagf;
using namespace libpetey;

int main(int argc, char **argv) {
  char *clsfile;
  char *confile;

  FILE *ifs;
  FILE *ofs1, *ofs2;

  char **line;
  long nline;			//number of lines
  nel_ta n;			//number of results

  linked_list<cls_ta> label;	//list of labels
  cls_ta *label_list;		//list of labels

  long ncls1;			//number of class labels
  cls_ta ncls2;			//max label+1
  real_a pdf1, *pdf2;		//cond. prob. read in; rearranged
  cls_ta *cls;			//output classes
  real_a *con;			//output confidence ratings
  int err;

  int strptr;
  int32_t stradv;

  if (argc < 3) {
    printf("Usage: svmout2agf SVMclass AGFbase\n");
    printf("\nConverts classification output from LIBSVM to libAGF format\n");
    printf("\nwhere:\n\n");
    printf("SVMclass is the output from LIBSVM\n");
    printf("AGFbase  is the base name for libAGF class output:\n");
    printf("            .cls for classes, .con for confidence ratings\n");
    exit(FATAL_COMMAND_OPTION_PARSE_ERROR);
  }

  //generate the output file names:
  clsfile=new char[strlen(argv[2])+5];
  strcpy(clsfile, argv[2]);
  strcat(clsfile, ".cls");

  confile=new char[strlen(argv[2])+5];
  strcpy(confile, argv[2]);
  strcat(confile, ".con");

  ifs=fopen(argv[1], "r");

  line=read_ascii_all(ifs, &nline);
  fclose(ifs);

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

    n=nline-1;
    cls=new cls_ta[n];
    con=new real_a[n];
    pdf2=new real_a[ncls2];
    for (cls_ta j=0; j<ncls2; j++) pdf2[j]=0;
    delete [] line[0];

    for (nel_ta i=0; i<n; i++) {
      err=sscanf(line[i+1], "%d%n", cls+i, &strptr);
      if (err != 1) {
        n=i;
        break;
      }
      for (long j=0; j<ncls1; j++) {
        err=sscanf(line[i+1]+strptr, "%g%n", &pdf1, &stradv);
        if (err != 1) {
          fprintf(stderr, "svmout2agf: read error--missing data element %d on line %d\n", j, i+2); 
          exit(FILE_READ_ERROR);
        }
        strptr+=stradv;
        pdf2[label_list[j]]=pdf1;
      }
      
      con[i]=(ncls1*pdf2[cls[i]]-1)/(ncls1-1);
      //printf("%d %f\n", cls, con);
      for (cls_ta j=0; j<ncls2; j++) printf("%10.6f ", pdf2[j]);
      printf("%5d\n", cls[i]);
      delete [] line [i+1];
    }

    ofs2=fopen(confile, "w");
    fwrite(con, sizeof(real_a), n, ofs2);
    fclose(ofs2);

    delete [] con;
    delete [] pdf2;

    for (nel_ta i=n+1; i<nline; i++) delete [] line[i];
  } else if (isdigit(line[0][strptr])) {
    n=nline;
    cls=new cls_ta[n];
    con=new real_a[n];
    for (nel_ta i=0; i<n; i++) {
      err=sscanf(line[i], "%d %g", cls+i, con+i);
      if (err==1) {
        con[i]=-1;
      } else if (err!=2) {
        n=i;
        break;
      }
      delete [] line[i];
    }
    for (nel_ta i=n; i<nline; i++) delete [] line[i];
    if (con[0]!=-1) {
      ncls2=1;
      for (nel_ta i=0; i<n; i++) if (cls[i]>=ncls2) ncls2=cls[i]+1;
      for (nel_ta i=0; i<n; i++) con[i]=(ncls2*con[i]-1)/(ncls2-1);
      ofs2=fopen(confile, "w");
      fwrite(con, sizeof(real_a), n, ofs2);
      fclose(ofs2);
    }
    delete [] con;
  } else {
    fprintf(stderr, "svmout2agf: format error on line 1\n"); 
    exit(FILE_READ_ERROR);
  }

  ofs1=fopen(clsfile, "w");
  fwrite(cls, sizeof(cls_ta), n, ofs1);
  fclose(ofs1);

  //printf("%d %d\n", n, nline);

  delete [] line;
  delete [] cls;

  delete [] clsfile;
  delete [] confile;
    
}

