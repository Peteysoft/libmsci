#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include <vector>

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

  template <class real, class cls_t>
  int select_classes(real **x, cls_t *cls, int n, cls_t *select, real **x1, cls_t *cls1) {
    real result;
    int ncls1=0;
    for (int i=0; i<n; i++) {
      if (select[cls[i]]) {
        x1[ncls1]=x[i];
	cls1[ncls1]=cls[i];
        ncls1++;
      }
    }
    return ncls1;
  }

  template <class real>
  real hausdorff_metric(real **x0, real **x1, int ncls0, int ncls1, int D) {
    real dmin0[ncls0];
    real dmin1[ncls1];
    real d, result;
    for (int i=0; i<ncls0; i++) dmin0[i]=metric2(x0[i], x1[0], D);
    for (int i=0; i<ncls1; i++) dmin1[i]=metric2(x0[0], x1[i], D);
    for (int i=0; i<ncls0; i++) {
      for (int j=0; j<ncls1; j++) {
        d=metric2(x0[i], x1[j], D);
	if (d<dmin0[i]) dmin0[i]=d;
	if (d<dmin1[j]) dmin1[j]=d;
      }
    }
    result=0;
    for (int i=0; i<ncls0; i++) if (dmin0[i]>result) result=dmin0[i];
    for (int i=0; i<ncls1; i++) if (dmin1[i]>result) result=dmin1[i];
    return result;
  }

  template <class real>
  real hausdorff_width(real **x, int n, int D) {
    real dmin[n];
    real d, result;
    result=metric2(x[0], x[1], D);
    for (int i=0; i<n; i++) {
      for (int j=0; j<i; j++) {
        d=metric2(x[i], x[j], D);
	if (d<dmin[i]) dmin[i]=d;
	if (d<dmin[j]) dmin[j]=d;
      }
    }
    result=0;
    for (int i=0; i<n; i++) if (dmin[i]>result) result=dmin[i];
    return result;
  }

  template <class real, class cls_t>
  real * class_triangle(real **x, cls_t *cls, int n, int D) {
    int ncls=0;
    real *d;
    int k=0;
    cls_t *map;
    real *x0[n], *x1[n];
    cls_t cls0[n], cls1[n];
    int ncls0, ncls1;
    for (int i=0; i<n; i++) if (cls[i]>=ncls) ncls=cls[i]+1;
    map=new cls_t[ncls];
    d=new real[ncls*(ncls-1)/2];
    for (int i=0; i<ncls; i++) map[i]=0;
    for (int i=0; i<ncls; i++) {
      for (int j=0; j<i; j++) {
        map[i]=1;
	ncls0=select_classes(x, cls, n, map, x0, cls0);
	map[i]=0;
	map[j]=1;
	ncls1=select_classes(x, cls, n, map, x1, cls1);
	map[j]=0;
	d[k]=sqrt(hausdorff_metric(x0, x1, ncls0, ncls1, D));
	printf("%12.4g ", d[k]);
	k++;
      }
      map[i]=1;
      ncls0=select_classes(x, cls, n, map, x0, cls0);
      map[i]=0;
      printf("%12.4g\n", sqrt(hausdorff_width(x0, ncls0, D)));
    }
    delete [] map;
    return d;
  }

  template float * class_triangle<float, int32_t>(float **, int32_t *, int, int);

  //a couple of helper functions:
  template <class vector_t>
  void random_coding_row(vector_t &coding_row, int n, int strictflag) {
    for (int i=0; i<n; i++) {
      coding_row[i]=(1+strictflag)*int(ranu()*(3-strictflag))-1;
    }
  }

  template <class vector_t>
  int check_coding_row(vector_t &coding_row, int n, int nt=-1) {
    int c1flag=0;
    int c2flag=0;
    for (int i=0; i<n; i++) {
      if (coding_row[i]==-1) c1flag++;
      if (coding_row[i]==1) c2flag++;
    }
    if (c1flag > 0 && c2flag > 0 && (nt<0 || c1flag+c2flag==nt)) return 1;
    return 0;
  }

  //print out common control files:
  int ** one_against_all(int ncls) {
    int **coding_matrix;
    coding_matrix=new int *[ncls];
    coding_matrix[0]=new int[ncls*ncls];
    for (int i=0; i<ncls; i++) {
      coding_matrix[i]=coding_matrix[0]+i*ncls;
      for (int j=0; j<i; j++) coding_matrix[i][j]=-1;
      for (int j=i+1; j<ncls; j++) coding_matrix[i][j]=-1;
      coding_matrix[i][i]=1;
    }
    return coding_matrix;
  }

  int ** one_against_one(int ncls) {
    int **coding_matrix;
    int k=0;
    int nrow=(ncls-1)*ncls/2;
    coding_matrix=new int *[nrow];
    coding_matrix[0]=new int[nrow*ncls];
    for (int i=0; i<ncls; i++) {
      for (int j=i+1; j<ncls; j++) {
        coding_matrix[k]=coding_matrix[0]+k*ncls;
	for (int m=0; m<ncls; m++) coding_matrix[k][m]=0;
	coding_matrix[k][i]=-1;
	coding_matrix[k][j]=1;
	k++;
      }
    }
    return coding_matrix;
  }

  int ** partition_adjacent(int ncls) {
    int **coding_matrix;
    coding_matrix=new int *[ncls-1];
    coding_matrix[0]=new int[ncls*(ncls-1)];
    for (int i=1; i<ncls; i++) {
      coding_matrix[i-1]=coding_matrix[0]+(i-1)*ncls;
      for (int j=0; j<i; j++) coding_matrix[i-1][j]=-1;
      for (int j=i; j<ncls; j++) coding_matrix[i-1][j]=1;
    }
    return coding_matrix;
  }

  //need to design a version of this for larger n:
  int ** random_coding_matrix(int ncls, int &ntrial, int strictflag) {
    tree_lg<vector<int> > used;
    vector<int> row(ncls);
    int **coding_matrix;
    int err1, err2;
    double nperm=pow(3-strictflag, ncls);
    coding_matrix=new int *[ntrial];
    coding_matrix[0]=new int[ncls*ntrial];

    for (int i=0; i<ntrial; i++) {
      coding_matrix[i]=coding_matrix[0]+i*ncls;
      do {
        random_coding_row(row, ncls, strictflag);
	err1=used.add_member(row);
	//remember that the same row with signs reversed is equivalent:
	for (int j=0; j<ncls; j++) row[j]=-row[j];
	err2=used.add_member(row);
        if (used.nel() >= nperm) {
          ntrial=i;
          coding_matrix[ntrial]=NULL;
	  goto done;
        }
      } while(err1<0 || err2<0 || check_coding_row(row, ncls)==0);
      for (int j=0; j<ncls; j++) coding_matrix[i][j]=row[j];
    }
    done: return coding_matrix;
  }

  int ** exhaustive_coding_matrix(int ncls) {
    int64_t nrow;
    int **coding_matrix;

    nrow=pow(2, ncls-1)-1;

    coding_matrix=new int *[nrow];
    coding_matrix[0]=new int[nrow*ncls];

    for (int64_t i=0; i<nrow; i++) {
      //bit_array *tobits=new bit_array((word *) &i, (sizeof(long)+1)/sizeof(word), ncls-1);
      bitset<sizeof(i)*8> *tobits=new bitset<sizeof(i)*8>(i);
      coding_matrix[i]=coding_matrix[0]+i*ncls;
      for (int j=0; j<ncls-1; j++) {
        coding_matrix[i][j+1]=2*(*tobits)[j]-1;
      }
      coding_matrix[i][0]=1;
      delete tobits;
    }
    return coding_matrix;
  }

  //generate orthogonal coding matrix:
  int ** ortho_coding_matrix_nqbf(int n, int strictflag, int toprow1) {
    tree_lg<vector<int> > list0;
    vector<int> trial0(n);
    double nperm0=pow(3-strictflag, n);
    int **coding_matrix;
    double eps=1e-12;
    int nfilled;

    //to wor properly this flag must be strictly 0 or 1:
    assert(strictflag==0 || strictflag==1);

    //allocate the coding matrix:
    coding_matrix=new int *[n];
    coding_matrix[0]=new int[n*n];
    for (int i=1; i<n; i++) coding_matrix[i]=coding_matrix[0]+n*i;

    do {
      int i;

      if (toprow1) {
        //top row is all 1's:
        for (int i=0; i<n; i++) {
          trial0[i]=1;
          coding_matrix[0][i]=1;
        }
      } else {
        //create a random first row:
        do {
          random_coding_row(coding_matrix[0], n, strictflag);
          //keep track of all initial rows:
          for (int i=0; i<n; i++) trial0[i]=coding_matrix[0][i];
        } while (check_coding_row(coding_matrix[0], n)==0 || 
		      list0.add_member(trial0) < 0);
      }

      for (i=1; i<n; i++) {
        //list of partial vectors already tried:
        tree_lg<vector<int> > *list=new tree_lg<vector<int> >;
        long ntried;		//number of vectors tried
        //allocate GSL data structures for linear system:
        gsl_matrix *a=gsl_matrix_alloc(i, i);
        gsl_vector *b=gsl_vector_alloc(i);
        gsl_vector *x=gsl_vector_alloc(i);
        //allocate space for singular value decomposition:
        gsl_matrix *vt=gsl_matrix_alloc(i, i);
        gsl_vector *s=gsl_vector_alloc(i);
        gsl_vector *work=gsl_vector_alloc(i);
        //allocate partial trial vector:
        vector<int> *trial=new vector<int>(n-i);
        //maximum number of permutations:
        double nperm=pow(3-strictflag, n-i);
	//have to decide which columns to select so that matrix is
	//non-singular:
	double det;		//determinant
	//permutations of the classes:
	long *perm=NULL;
	long inv[n];			//inverse mapping
	vector<int> selvec;		//selection vector
	//list of selection vectors to avoid repeats:
	tree_lg<vector<int> > *selist=new tree_lg<vector<int> >;

        //create linear system to solve:
        do {
          if (perm!=NULL) delete [] perm;
          perm=randomize(n);
          for (int j=0; j<i; j++) selvec[perm[j]]=0;
          for (int j=i; j<n; j++) selvec[perm[j]]=1;
          det=0;
          if (selist->add_member(selvec)<0) continue;

          for (int j=0; j<i; j++) {
            for (int k=0; k<i; k++) {
              gsl_matrix_set(a, j, k, coding_matrix[j][perm[k]]);
            }
          }
          //find the singular value decomposition:
	  printf("matrix:\n");
	  print_gsl_matrix(stdout, a);
          gsl_linalg_SV_decomp(a, vt, s, work);
          det=1;
	  printf("singular values:\n");
	  for (int j=0; j<i; j++) {
            printf("%10.4g ", gsl_vector_get(s, j));
            det*=gsl_vector_get(s, j);
	  }
	  printf("\n");
        } while (det<eps);

	for (int j=0; j<n; j++) {
          printf("%3d ", perm[j]);
          inv[perm[j]]=j;
	}
	printf("\n");
	printf("\n");

        try_again:	//more spaghetti code... yah!
          //create a random partial trial vector:
          do {
            //lets write some spaghetti code:
            if (list->nel() >= nperm) goto done;
            random_coding_row(*trial, n-i, strictflag);
          } while (list->add_member(*trial) < 0);

          //calculate solution vector:
          for (int j=0; j<i; j++) {
            double temp=0;
            for (int k=0; k<n-i; k++) {
              temp+=(*trial)[k]*coding_matrix[j][perm[k+i]];
            }
            gsl_vector_set(b, j, -temp);
          }

          //solve the linear system to get the rest of the potential new row
          //for the coding matrix:
          gsl_linalg_SV_solve(a, vt, s, b, x);

	  for (int j=0; j<n-i; j++) printf("%2d ", (*trial)[j]);
	  for (int j=0; j<i; j++) printf("%10.4g ", gsl_vector_get(x, j));
          printf("\n");

          //check that the rest of the vector has values of -1 or 1 (or 0):
          for (int j=0; j<i; j++) {
            double val=gsl_vector_get(x, j);
            if (fabs(val-1) > eps && fabs(val+1) > eps && 
			  (fabs(val) > eps || strictflag)) goto try_again;
            coding_matrix[i][perm[j]]=round(val);
          }

          for (int j=0; j<n-i; j++) {
            coding_matrix[i][perm[i+j]]=(*trial)[j];
          }

        if (check_coding_row(coding_matrix[i], n)==0) goto try_again;
        printf("\n");

        //clean up:
        gsl_matrix_free(a);
        gsl_vector_free(b);
        gsl_vector_free(x);
        gsl_matrix_free(vt);
        gsl_vector_free(s);
        gsl_vector_free(work);
        delete list;
	delete [] perm;
	delete selist;
      }
      done: nfilled=i;

      for (int i=0; i<nfilled; i++) {
        for (int j=0; j<n; j++) printf("%2d ", coding_matrix[i][j]);
        printf("\n");
      }
      printf("\n");

      for (int i=0; i<nfilled; i++) {
        for (int j=0; j<nfilled; j++) {
          int temp=0;
          for (int k=0; k<n; k++) temp+=coding_matrix[i][k]*coding_matrix[j][k];
          printf("%3d ", temp);
        }
        printf("\n");
      }
      printf("\n");
    //} while (nfilled<n && list0.nel()<nperm0);
    } while (0);

    if (nfilled<n) coding_matrix[nfilled]=NULL;

    return coding_matrix;
  }

  //need to design a more efficient version of this...
  int ** ortho_coding_matrix_brute_force(int n, int toprow1) {
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

    if (toprow1) {
      for (int i=0; i<n; i++) coding_matrix[0][i]=1;
      nfilled=1;
    } else {
      nfilled=0;
    }
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

    if (nfilled<n) coding_matrix[nfilled]=NULL;

    for (int i=0; i<nfilled; i++) {
      for (int j=0; j<n; j++) printf("%2d ", coding_matrix[i][j]);
      printf("\n");
    }
    printf("\n");

    for (int i=0; i<nfilled; i++) {
      for (int j=0; j<nfilled; j++) {
        int temp=0;
        for (int k=0; k<n; k++) temp+=coding_matrix[i][k]*coding_matrix[j][k];
        printf("%3d ", temp);
      }
      printf("\n");
    }
    printf("\n");

    delete [] trial;
    return coding_matrix;
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

  void print_control_nonhier(FILE *fs, int **coding_matrix, int n, int ncls, const char *options) {
    for (int i=0; i<n; i++) {
      if (coding_matrix[i]==NULL) break;
      if (options!=NULL) fprintf(fs, "\"%s\" ", options); else fprintf(fs, "\"\" ");
      for (int j=0; j<ncls; j++) {
        if (coding_matrix[i][j]<0) fprintf(fs, "%d ", j);
      }
      fprintf(fs, "%c", PARTITION_SYMBOL);
      for (int j=0; j<ncls; j++) {
        if (coding_matrix[i][j]>0) fprintf(fs, " %d", j);
      }
      fprintf(fs, ";\n");
    }
    fprintf(fs, "{");
    for (int i=0; i<ncls; i++) fprintf(fs, " %d", i);
    fprintf(fs, "}\n");
  }

}

