#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>

//#include "full_util.h"
//#include "bit_array.h"
#include <bitset>
#include "gsl_util.h"
#include "error_codes.h"
#include "randomize.h"
#include "tree_lg.h"

#include "agf_lib.h"

using namespace libpetey;
using namespace std;

namespace libagf {

  //print out common control files:
  void one_against_all(FILE *fs, int ncls, const char *options) {
    for (int i=0; i<ncls; i++) {
      if (options!=NULL) {
        if (i==0) fprintf(fs, "\"%s\"", options); else fprintf(fs, "\".\"");
      } else {
        fprintf(fs, "\"\"");
      }
      fprintf(fs, " %d %c", i, PARTITION_SYMBOL);
      for (int j=0; j<i; j++) fprintf(fs, " %d", j);
      for (int j=i+1; j<ncls; j++) fprintf(fs, " %d", j);
      fprintf(fs, ";\n");
    }
  }

  void one_against_one(FILE *fs, int ncls, const char *options) {
    for (int i=0; i<ncls; i++) {
      for (int j=i+1; j<ncls; j++) {
        if (options!=NULL) {
          if (i==0) fprintf(fs, "\"%s\"", options); else fprintf(fs, "\".\"");
        } else {
          fprintf(fs, "\"\"");
        }
        fprintf(fs, " %d %c %d;\n", i, PARTITION_SYMBOL, j);
      }
    }
  }

  void partition_adjacent(FILE *fs, int ncls, const char *options) {
    for (int i=1; i<ncls; i++) {
      if (options!=NULL) {
        if (i==0) fprintf(fs, "\"%s\"", options); else fprintf(fs, "\".\"");
      } else {
        fprintf(fs, "\"\"");
      }
      for (int j=0; j<i; j++) fprintf(fs, " %d", j);
      fprintf(fs, " %c", PARTITION_SYMBOL);
      for (int j=i; j<ncls; j++) fprintf(fs, " %d", j);
      fprintf(fs, ";\n");
    }
  }

  //need to design a version of this for larger n:
  void random_coding_matrix(FILE *fs, int ncls, int ntrial, int strictflag) {
    tree_lg<int64_t> used;
    int row[ncls];
    int64_t cur, curi;
    int mult;
    int list1[ncls];
    int list2[ncls];
    int n1, n2;

    for (int i=0; i<ntrial; i++) {
      do {
        n1=0;
        n2=0;
        mult=1;
	cur=0;
	curi=0;
        for (int j=0; j<ncls; j++) {
          if (strictflag) {
            row[j]=2*ranu();
            cur+=mult*row[j];
            curi+=mult*(1-row[j]);
            row[j]=2*row[j]-1;
            mult*=2;
          } else {
            row[j]=3*ranu();
            cur+=mult*row[j];
            curi+=mult*(2-row[j]);
            row[j]=row[j]-1;
            mult*=3;
	  }
          if (row[j]<0) {
            list1[n1]=j;
            n1++;
          } else if (row[j]>0) {
            list2[n2]=j;
            n2++;
          }
	}
      } while(used.add_member(cur)<0 || used.add_member(curi)<0 || n1<1 || n2<1);
      fprintf(fs, "\"\"");
      for (int j=0; j<n1; j++) fprintf(fs, " %d", list1[j]);
      fprintf(fs, " %c", PARTITION_SYMBOL);
      for (int j=0; j<n2; j++) fprintf(fs, " %d", list2[j]);
      fprintf(fs, ";\n");
    }
  }

  void exhaustive_coding_matrix(FILE *fs, int ncls) {
    int64_t nrow;

    nrow=pow(2, ncls-1)-1;

    for (int64_t i=0; i<nrow; i++) {
      //bit_array *tobits=new bit_array((word *) &i, (sizeof(long)+1)/sizeof(word), ncls-1);
      bitset<sizeof(i)*8> *tobits=new bitset<sizeof(i)*8>(i);
      fprintf(fs, "\"\" ");
      for (long j=0; j<ncls-1; j++) {
        if ((*tobits)[j]==0) fprintf(fs, "%d ", j+1);
      }
      fprintf(fs, "%c ", PARTITION_SYMBOL);
      for (long j=0; j<ncls-1; j++) {
        if ((*tobits)[j]) fprintf(fs, "%d ", j+1);
      }
      fprintf(fs, "0;\n");
      delete tobits;
    }
  }

  template <class vector_t>
  void random_coding_row(vector_t &coding_row, int n, int strictflag) {
    for (int i=0; i<n; i++) {
      coding_row[i]=(1+strictflag)*int(ranu()*(3-strictflag))-1;
    }
  }

  int check_coding_row(int *coding_row, int n) {
    int c1flag=0;
    int c2flag=0;
    for (int i=1; i<n; i++) {
      if (coding_row[i]==-1) c1flag=1;
      if (coding_row[i]==1) c2flag=1;
    }
    if (c1flag && c2flag) return 1;
    return 0;
  }

  //generate orthogonal coding matrix:
  void ortho_coding_matrix_nqbf(FILE *fs, int n, int strictflag) {
    int **coding_matrix;
    double eps=1e-12;

    //to wor properly this flag must be strictly 0 or 1:
    assert(strictflag==0 || strictflag==1);

    //allocate the coding matrix:
    coding_matrix=new int *[n];
    coding_matrix[0]=new int[n*n];
    for (int i=1; i<n; i++) coding_matrix[i]=coding_matrix[0]+n*i;

    //create a random first row:
    do {
      random_coding_row(coding_matrix[0], n, strictflag);
    } while (check_coding_row(coding_matrix[0], n)==0);

    for (int i=1; i<n; i++) {
      //list of partial vectors already tried:
      tree_lg<vector_s<int> > *list=new tree_lg<vector_s<int> >;
      //allocate GSL data structures for linear system:
      gsl_matrix *a=gsl_matrix_alloc(i, i);
      gsl_vector *b=gsl_vector_alloc(i);
      gsl_vector *x=gsl_vector_alloc(i);
      //allocate space for singular value decomposition:
      gsl_matrix *vt=gsl_matrix_alloc(i, i);
      gsl_vector *s=gsl_vector_alloc(i);
      gsl_vector *work=gsl_vector_alloc(i);
      //allocate partial trial vector:
      vector_s<int> *trial=new vector_s<int>(n-i);
      //maximum number of permutations:
      double nperm=pow(3-strictflag, n);

      try_again:	//more spaghetti code... yah!
        //create linear system to solve:
        for (int j=0; j<i; j++) {
          for (int k=0; k<i; k++) {
            gsl_matrix_set(j, k, coding_matrix[j][k]);
          }
        }
        //find the singular value decomposition:
        gsl_linalg_SV_decomp(a, vt, s, work);

        //create a random partial trial vector:
        do {
          //lets write some spaghetti code:
          if (list->nel() >= nperm) goto done;
	  random_coding_row(trial, n, strictflag);
        } while (list->add_member(trial) == list->nel());

        //calculate solution vector:
        for (int j=0; j<i; j++) {
          double temp=0;
          for (int k=0; k<n-i; k++) {
            temp+=coding_matrix[j][k+i]*trial[k];
          }
          gsl_vector_set(b, -temp);
        }

        //solve the linear system to get the rest of the potential new row
        //for the coding matrix:
        gsl_linalg_SV_solve(a, vt, s, b, x);

        //check that the rest of the vector has values of -1 or 1 (or 0):
	for (int j=0; j<i; j++) {
          double val=gsl_vector_get(b, j);
          if (fabs(val-1) > eps && fabs(val+1) > eps && 
			  (fabs(val) > eps || strictflag)) goto try_again;
	  coding_matrix[i][j]=val;
	}

        for (int j=0; j<n-i; j++) {
          coding_matrix[i][i+j]=trial[j];
        }

      if (check_coding_row(coding_matrix[i])==0) goto try_again;

      //clean up:
      gsl_matrix_free(a);
      gsl_vector_free(b);
      gsl_vector_free(x);
      gsl_matrix_free(vt);
      gsl_vector_free(s);
      gsl_vector_free(work);
      delete list;
      delete trial;
    }

    for (int i=0; i<nfilled; i++) {
      fprintf(fs, "\"\" ");
      for (int j=0; j<n; j++) {
        if (coding_matrix[i][j]<0) fprintf(fs, "%d ", j);
      }
      fprintf(fs, "%c", PARTITION_SYMBOL);
      for (int j=0; j<n; j++) {
        if (coding_matrix[i][j]>0) fprintf(fs, " %d", j);
      }
      fprintf(fs, ";\n");
    }

    delete [] coding_matrix[0];
    delete [] coding_matrix;

  }

  //need to design a more efficient version of this...
  void ortho_coding_matrix_brute_force(FILE *fs, int n) {
    long *trial;
    //bit_array *tobits;
    bitset<sizeof(int)*8> *tobits;
    int **coding_matrix;
    int nfilled;
    int nperm;

    nperm=pow(2, n);

    coding_matrix=new int *[n];
    coding_matrix[0]=new int[n*n];
    for (int i=1; i<n; i++) coding_matrix[i]=coding_matrix[0]+n*i;

    if (n>8*sizeof(long)-1) {
      fprintf(stderr, "Size must be %d or less\n", 8*sizeof(long)-1);
      exit(PARAMETER_OUT_OF_RANGE);
    }

    trial=randomize(nperm);

    nfilled=0;
    for (int i=0; i<nperm; i++) {
      int dprod;
      //printf("%d\n", trial[i]);
      //tobits=new bit_array((word *) (trial+i), (sizeof(long)+1)/sizeof(word), n);
      //all zeros or all ones not allowed:
      if (trial[i]==0 || trial[i]==nperm-1) continue;
      tobits=new bitset<sizeof(i)*8>(trial[i]);
      for (long j=0; j<n; j++) {
        coding_matrix[nfilled][j]=2*(long) (*tobits)[j]-1;
        //printf("%d ", (int32_t) (*tobits)[j]);
      }
      //printf("\n");
      //for (int j=0; j<n; j++) printf("%3d", coding_matrix[nfilled][j]);
      //printf("\n");
      dprod=0;
      for (int j=0; j<nfilled; j++) {
        dprod=0;
        for (int k=0; k<n; k++) {
          dprod+=coding_matrix[nfilled][k]*coding_matrix[j][k];
        }
        if (dprod!=0) break;
      }
      if (dprod==0) {
        //for (int j=0; j<n; j++) printf("%3d", coding_matrix[nfilled][j]);
        //printf("\n");
        nfilled++;
      }
      if (nfilled>=n) break;
      delete tobits;
    }

    for (int i=0; i<nfilled; i++) {
      fprintf(fs, "\"\" ");
      for (int j=0; j<n; j++) {
        if (coding_matrix[i][j]<0) fprintf(fs, "%d ", j);
      }
      fprintf(fs, "%c", PARTITION_SYMBOL);
      for (int j=0; j<n; j++) {
        if (coding_matrix[i][j]>0) fprintf(fs, " %d", j);
      }
      fprintf(fs, ";\n");
    }

    delete [] trial;
    delete [] coding_matrix[0];
    delete [] coding_matrix;
  }

  void print_control_hier(FILE *fs, int ncls, int c0, int depth) {
    for (int i=0; i<depth; i++) fprintf(fs, "  ");
    if (ncls==1) {
      fprintf(fs, "%d\n", c0);
    } else if (ncls>1) {
      fprintf(fs, "\"\" {\n");
      print_control_hier(fs, ncls/2, c0, depth+1);
      print_control_hier(fs, ncls-ncls/2, c0+ncls/2, depth+1);
      fprintf(fs, "\n");
      for (int i=0; i<depth; i++) fprintf(fs, "  ");
      fprintf(fs, "}\n");
    }
  }

  void print_control_1vsall(FILE *fs, int ncls, const char *opt) {
    one_against_all(fs, ncls, opt);
    fprintf(fs, "{");
    for (int i=0; i<ncls; i++) fprintf(fs, " %d", i);
    fprintf(fs, "}\n");
  }

  void print_control_1vs1(FILE *fs, int ncls, const char *opt) {
    one_against_one(fs, ncls, opt);
    fprintf(fs, "{");
    for (int i=0; i<ncls; i++) fprintf(fs, " %d", i);
    fprintf(fs, "}\n");
  }

  void print_control_adj(FILE *fs, int ncls, const char *opt) {
    partition_adjacent(fs, ncls, opt);
    fprintf(fs, "{");
    for (int i=0; i<ncls; i++) fprintf(fs, " %d", i);
    fprintf(fs, "}\n");
  }

  void print_control_random(FILE *fs, int ncls, int nrow, int strictflag) {
    random_coding_matrix(fs, ncls, nrow, strictflag);
    fprintf(fs, "{");
    for (int i=0; i<ncls; i++) fprintf(fs, " %d", i);
    fprintf(fs, "}\n");
  }

  void print_control_exhaustive(FILE *fs, int ncls) {
    exhaustive_coding_matrix(fs, ncls);
    fprintf(fs, "{");
    for (int i=0; i<ncls; i++) fprintf(fs, " %d", i);
    fprintf(fs, "}\n");
  }

  void print_control_ortho(FILE *fs, int ncls) {
    ortho_coding_matrix_brute_force(fs, ncls);
    fprintf(fs, "{");
    for (int i=0; i<ncls; i++) fprintf(fs, " %d", i);
    fprintf(fs, "}\n");
  }

}

