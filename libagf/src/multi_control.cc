#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include <vector>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>

#include "full_util.h"
//#include "bit_array.h"
#include <bitset>
#include "gsl_util.h"
#include "error_codes.h"
#include "randomize.h"
#include "tree_lg.h"
#include "kselect.h"

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

  //calculate the Hausdorff distance between two finite sets:
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

  //calculate the Hausdorf metric between the two partitions in a row of 
  //a coding matrix
  template <typename real, typename scalar>
  real partition_distance(real *d, int n, scalar *coding_row) {
    int n1, n2;		//number in each partition
    int ind1[n], ind2[n];		//indices of each partition
    int k;				//1-d subscript

    real *dmin1;
    real *dmin2;
    real result;
    n1=0; n2=0;
    for (int i=0; i<n; i++) {
      if (coding_row[i]<0) {
        ind1[n1]=i;
        n1++;
      } else if (coding_row[i]>0) {
        ind2[n2]=i;
        n2++;
      }
    }
    dmin1=new real[n1];
    dmin2=new real[n2];

    for (int i=0; i<n1; i++) dmin1[i]=-1;
    for (int i=0; i<n2-1; i++) dmin2[i]=-1;

    for (int i=0; i<n1; i++) {
      for (int j=0; j<n2; j++) {
        if (ind1[i]>ind2[j]) {
          k=ind1[i]*(ind1[i]-1)/2+ind2[j];
	} else {
          k=ind2[j]*(ind2[j]-1)/2+ind1[i];
	}
	if (dmin1[i]<0 || d[k]<dmin1[i]) dmin1[i]=d[k];
	if (dmin2[j]<0 || d[k]<dmin2[j]) dmin2[j]=d[k];
      }
    }

    result=0;
    for (int i=0; i<n1; i++) if (dmin1[i]>result) result=dmin1[i];
    for (int i=0; i<n2; i++) if (dmin2[i]>result) result=dmin2[i];

    delete [] dmin1;
    delete [] dmin2;

    return result;

  }

  //calculate the Hausdorf metric between the two partitions in a row of 
  //a coding matrix
  template <typename real, typename scalar, typename cls_t>
  real partition_distance(real *d, 		//distance triangle (all points)
		  cls_t *cls,			//classes
		  nel_ta n,			//number of points
		  scalar *coding_row) {		//row of coding matrix
    int n1, n2;		//number in each partition
    int ind1[n], ind2[n];		//indices of each partition
    int k;				//1-d subscript

    real *dmin1;
    real *dmin2;
    real result;
    n1=0; n2=0;
    for (int i=0; i<n; i++) {
      if (coding_row[cls[i]]<0) {
        ind1[n1]=i;
        n1++;
      } else if (coding_row[cls[i]]>0) {
        ind2[n2]=i;
        n2++;
      }
    }
    dmin1=new real[n1];
    dmin2=new real[n2];

    for (int i=0; i<n1; i++) dmin1[i]=-1;
    for (int i=0; i<n2-1; i++) dmin2[i]=-1;

    for (int i=0; i<n1; i++) {
      for (int j=0; j<n2; j++) {
        if (ind1[i]>ind2[j]) {
          k=ind1[i]*(ind1[i]-1)/2+ind2[j];
	} else {
          k=ind2[j]*(ind2[j]-1)/2+ind1[i];
	}
	if (dmin1[i]<0 || d[k]<dmin1[i]) dmin1[i]=d[k];
	if (dmin2[j]<0 || d[k]<dmin2[j]) dmin2[j]=d[k];
      }
    }

    result=0;
    for (int i=0; i<n1; i++) if (dmin1[i]>result) result=dmin1[i];
    for (int i=0; i<n2; i++) if (dmin2[i]>result) result=dmin2[i];

    delete [] dmin1;
    delete [] dmin2;

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

  //builds a distance triangle for the classes based on the Hausdorff
  //metric
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
	//printf("%12.4g ", d[k]);
	k++;
      }
      map[i]=1;
      ncls0=select_classes(x, cls, n, map, x0, cls0);
      map[i]=0;
      //printf("%12.4g\n", sqrt(hausdorff_width(x0, ncls0, D)));
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

  //generate common coding matrices:
  template <typename scalar>
  scalar ** one_against_all(int ncls) {
    scalar **coding_matrix;
    coding_matrix=new scalar *[ncls];
    coding_matrix[0]=new scalar[ncls*ncls];
    for (int i=0; i<ncls; i++) {
      coding_matrix[i]=coding_matrix[0]+i*ncls;
      for (int j=0; j<i; j++) coding_matrix[i][j]=-1;
      for (int j=i+1; j<ncls; j++) coding_matrix[i][j]=-1;
      coding_matrix[i][i]=1;
    }
    return coding_matrix;
  }

  template <typename scalar>
  scalar ** one_against_one(int ncls) {
    scalar **coding_matrix;
    int k=0;
    int nrow=(ncls-1)*ncls/2;
    coding_matrix=new scalar *[nrow];
    coding_matrix[0]=new scalar[nrow*ncls];
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

  template <typename scalar>
  scalar ** partition_adjacent(int ncls) {
    scalar **coding_matrix;
    coding_matrix=new scalar *[ncls-1];
    coding_matrix[0]=new scalar[ncls*(ncls-1)];
    for (int i=1; i<ncls; i++) {
      coding_matrix[i-1]=coding_matrix[0]+(i-1)*ncls;
      for (int j=0; j<i; j++) coding_matrix[i-1][j]=-1;
      for (int j=i; j<ncls; j++) coding_matrix[i-1][j]=1;
    }
    return coding_matrix;
  }

  //need to design a version of this for larger n:
  template <typename scalar>
  scalar ** random_coding_matrix(int ncls, int &ntrial, int strictflag) {
    tree_lg<vector<int> > used;
    vector<int> row(ncls);
    scalar **coding_matrix;
    int err1, err2;
    double nperm=pow(3-strictflag, ncls);
    coding_matrix=new scalar *[ntrial];
    coding_matrix[0]=new scalar[ncls*ntrial];

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

  //basic brute force:
  //- no attempt to keep track of trials
  //- not halting
  template <typename scalar>
  scalar **orthogonal_coding_matrix(int ncls, int nrow) {
    scalar **result_t;
    scalar **result;
    scalar dprod;
    int notdone;

    result_t=allocate_matrix<scalar>(ncls, nrow);

    for (int i=0; i<ncls; i++) {
      do {
        notdone=0;
        random_coding_row(result_t[i], nrow, 1);
	for (int j=0; j<i; j++) {
	  dprod=0;
	  for (int k=0; k<nrow; k++) dprod+=result_t[i][k]*result_t[j][k];
	  if (dprod!=0) {
            notdone=1;
	    break;
	  }
	}
      } while (notdone);
      print_matrix(stdout, result_t, ncls, nrow);
    }
    result=matrix_transpose(result_t, ncls, nrow);
    delete_matrix(result_t);
    return result;
  }



  template <typename scalar>
  scalar ** exhaustive_coding_matrix(int ncls) {
    int64_t nrow;
    scalar **coding_matrix;

    nrow=pow(2, ncls-1)-1;

    coding_matrix=new scalar *[nrow];
    coding_matrix[0]=new scalar[nrow*ncls];

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
  template <typename scalar>
  scalar ** ortho_coding_matrix_nqbf(int n, int strictflag, int toprow1) {
    tree_lg<vector<int> > list0;
    vector<int> trial0(n);
    double nperm0=pow(3-strictflag, n);
    scalar **coding_matrix;
    double eps=1e-12;
    int nfilled;

    //to wor properly this flag must be strictly 0 or 1:
    assert(strictflag==0 || strictflag==1);

    //allocate the coding matrix:
    coding_matrix=new scalar *[n];
    coding_matrix[0]=new scalar[n*n];
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
  //(greedy, not brute force... brute force down below...)
  template <typename scalar>
  scalar ** ortho_coding_matrix_greedy(int m, int n, int toprow1) {
    long *trial;
    //bit_array *tobits;
    bitset<sizeof(int)*8> *tobits;
    scalar **coding_matrix;
    int nfilled;
    int nperm;

    nperm=pow(2, n);

    coding_matrix=new scalar *[m];
    coding_matrix[0]=new scalar[m*n];
    for (int i=1; i<m; i++) coding_matrix[i]=coding_matrix[0]+n*i;

    if (n>8*sizeof(long)-1) {
      fprintf(stderr, "Size must be %d or less\n", 8*sizeof(long)-1);
      throw(PARAMETER_OUT_OF_RANGE);
    }

    do {
      if (toprow1) {
        for (int i=0; i<n; i++) coding_matrix[0][i]=1;
        nfilled=1;
      } else {
        nfilled=0;
      }
      trial=randomize(nperm);
      for (int i=0; i<nperm; i++) {
        int dprod;
	//printf("%d\n", i);
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
        delete tobits;
        if (nfilled>=m) break;
      }
    } while (nfilled<m);

    //if (nfilled<n) coding_matrix[nfilled]=NULL;

    delete [] trial;
    return coding_matrix;

    //debugging output:
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

  }

  //need to design a more efficient version of this...
  template <typename scalar>
  scalar ** ortho_coding_matrix_brute_force(int n) {
    long *trial;
    //bit_array *tobits;
    bitset<sizeof(int)*8> *tobits;
    scalar **coding_matrix;
    int nfilled;
    int nperm;
    int loc[n];
    int i;
    int dprod;

    nperm=pow(2, n);

    coding_matrix=new scalar *[n];
    coding_matrix[0]=new scalar[n*n];
    for (int i=1; i<n; i++) coding_matrix[i]=coding_matrix[0]+n*i;

    if (n>8*sizeof(long)-1) {
      fprintf(stderr, "Size must be %d or less\n", 8*sizeof(long)-1);
      throw(PARAMETER_OUT_OF_RANGE);
    }

    trial=randomize(nperm);
    loc[0]=0;

    for (i=0; i<n; i++) {
      if (i>0) loc[i]=loc[i-1]+1;
      for (loc[i]=0; loc[i]<nperm; loc[i]++) {
        if (trial[loc[i]]==0 || trial[loc[i]]==nperm-1) continue;
        tobits=new bitset<sizeof(i)*8>(trial[loc[i]]);
        for (long j=0; j<n; j++) {
          coding_matrix[i][j]=2*(long) (*tobits)[j]-1;
          //printf("%d ", (int32_t) (*tobits)[j]);
        }
        delete tobits;
        dprod=0;
        for (int j=0; j<i; j++) {
          dprod=0;
          for (int k=0; k<n; k++) {
            dprod+=coding_matrix[i][k]*coding_matrix[j][k];
          }
          if (dprod!=0) break;
        }
        if (dprod==0) {
          break;
          //for (int j=0; j<n; j++) printf("%3d", coding_matrix[nfilled][j]);
          //printf("\n");
        } 
      }
      if (dprod!=0 && loc[i]==nperm) {
        if (i==0) goto finish;
        i--;
      }
    }

    finish:
	
      if (i<n) coding_matrix[i]=NULL;
      delete [] trial;

    return coding_matrix;

    //debugging output:
    nfilled=i;
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

  }

  template <typename scalar, typename real>
  scalar **optimal_coding_matrix(int n, 	//number of classes
		  int nrow, 			//desired number of rows
		  real *d) {			//distance triangle btw. classes
    long nperm=pow(2, n-1);		//number of possible rows
    scalar **result;			//returned coding matrix
    int ndone;				//number of rows found
    scalar testrow[n];			//test row
    kiselect_base<real> *all;		//rows stored as integers
    long rows[nrow];			//rows as integers (2)
    real dmax[nrow];			//largest distances (not needed)
    long dum;				//for pulling out bits

    all=new kiselect_heap<real>(nrow);
    for (long i=1; i<nperm; i++) {
      //convert integer to an array of {-1, 1} values:
      dum=i;
      for (int j=0; j<n; j++) {
        testrow[j]=dum & 1;
        if (testrow[j]==0) testrow[j]=-1;
	dum=dum >> 1;
      }
      //measure partition distance and add to heap:
      all->add(-partition_distance(d, n, testrow), i);
    }
    result=allocate_matrix<scalar>(nrow, n);

    //extract the resultant coding matrix:
    all->get(dmax, rows);
    for (int i=0; i<nrow; i++) {
      dum=rows[i];
      for (int j=0; j<n; j++) {
        result[i][j]=dum & 1;
        if (result[i][j]==0) result[i][j]=-1;
	dum=dum >> 1;
      }
    }

    delete all;

    return result;
  }

  template <typename scalar, typename real, typename cls_t>
  scalar **optimal_coding_matrix(real *d,	//distance triangle (all points)
		  cls_t *cls,			//classes
		  nel_ta n, 			//number of points
		  int nrow) {			//desired number of rows
    long nperm;				//number of possible rows
    cls_t ncls;				//nuumber of classes
    scalar **result;			//returned coding matrix
    int ndone;				//number of rows found
    scalar testrow[n];			//test row
    kiselect_base<real> *all;		//rows stored as integers
    long rows[nrow];			//rows as integers (2)
    real dmax[nrow];			//largest distances (not needed)
    long dum;				//for pulling out bits

    ncls=1;
    for (nel_ta i=0; i<n; i++) if (cls[i]>=ncls) ncls=cls[i]+1;
    nperm=pow(2, ncls-1);

    all=new kiselect_heap<real>(nrow);
    for (long i=1; i<nperm; i++) {
      //convert integer to an array of {-1, 1} values:
      dum=i;
      for (int j=0; j<n; j++) {
        testrow[j]=dum & 1;
        if (testrow[j]==0) testrow[j]=-1;
	dum=dum >> 1;
      }
      //measure partition distance and add to heap:
      all->add(-partition_distance(d, cls, n, testrow), i);
    }
    result=allocate_matrix<scalar>(nrow, n);

    //extract the resultant coding matrix:
    all->get(dmax, rows);
    for (int i=0; i<nrow; i++) {
      dum=rows[i];
      for (int j=0; j<n; j++) {
        result[i][j]=dum & 1;
        if (result[i][j]==0) result[i][j]=-1;
	dum=dum >> 1;
      }
    }

    delete all;

    return result;
  }

  template <typename scalar>
  scalar ** hierarchical_nonhierarchical(int n) {
    scalar **result;
    if (n==1) return NULL;
    result=allocate_matrix<scalar>(n-1, n);
    if (n==2) {
      result[0][0]=-1;
      result[0][1]=1;
    } else {
      scalar **sub;
      for (int i=0; i<n/2; i++) {
        result[0][i]=-1;
      }
      for (int i=n/2; i<n; i++) {
        result[0][i]=1;
      }
      sub=hierarchical_nonhierarchical<scalar>(n/2);
      for (int i=0; i<n/2-1; i++) {
        for (int j=0; j<n/2; j++) result[i+1][j]=sub[i][j];
        for (int j=0; j<n-n/2; j++) result[i+1][j+n/2]=0;
      }
      if (sub!=NULL) delete_matrix(sub);
      sub=hierarchical_nonhierarchical<scalar>(n-n/2);
      for (int i=0; i<n-n/2-1; i++) {
        for (int j=0; j<n/2; j++) result[i+n/2][j]=0;
        for (int j=0; j<n-n/2; j++) result[i+n/2][j+n/2]=sub[i][j];
      }
      if (sub!=NULL) delete_matrix(sub);
    }
    return result;
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

  template <typename scalar>
  void print_control_nonhier(FILE *fs, scalar **coding_matrix, int n, int ncls, char **options, cls_ta *label) {
    for (int i=0; i<n; i++) {
      if (coding_matrix[i]==NULL) break;
      if (options!=NULL) fprintf(fs, "%s ", options[i]); else fprintf(fs, "\"\" ");
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
    if (label==NULL) {
      for (int i=0; i<ncls; i++) fprintf(fs, " %d", i);
    } else {
      for (int i=0; i<ncls; i++) fprintf(fs, " %d", label[i]);
    }
    fprintf(fs, "}\n");
  }

  template int **optimal_coding_matrix<int, float>(int, int, float *);
  template int **optimal_coding_matrix<int, double>(int, int, double *);
  template float **optimal_coding_matrix<float, float>(int, int, float *);
  template float **optimal_coding_matrix<float, double>(int, int, double *);
  template double **optimal_coding_matrix<double, float>(int, int, float *);
  template double **optimal_coding_matrix<double, double>(int, int, double *);

  template int **optimal_coding_matrix<int, float, cls_ta>(float *, cls_ta *, nel_ta, int);
  template int **optimal_coding_matrix<int, double, cls_ta>(double *, cls_ta *, nel_ta, int);

  template int ** one_against_all<int>(int);
  template float ** one_against_all<float>(int);
  template double ** one_against_all<double>(int);

  template int ** one_against_one<int>(int);
  template float ** one_against_one<float>(int);
  template double ** one_against_one<double>(int);

  template int ** partition_adjacent<int>(int);
  template float ** partition_adjacent<float>(int);
  template double ** partition_adjacent<double>(int);

  template int ** random_coding_matrix<int>(int, int &, int);
  template float ** random_coding_matrix<float>(int, int &, int);
  template double ** random_coding_matrix<double>(int, int &, int);

  template int ** exhaustive_coding_matrix<int>(int);
  template float ** exhaustive_coding_matrix<float>(int);
  template double ** exhaustive_coding_matrix<double>(int);

  template int ** ortho_coding_matrix_nqbf<int>(int,int,int);
  template float ** ortho_coding_matrix_nqbf<float>(int,int,int);
  template double ** ortho_coding_matrix_nqbf<double>(int,int,int);

  template int ** ortho_coding_matrix_greedy<int>(int, int, int);
  template float ** ortho_coding_matrix_greedy<float>(int, int, int);
  template double ** ortho_coding_matrix_greedy<double>(int, int, int);

  template int ** ortho_coding_matrix_brute_force<int>(int);
  template float ** ortho_coding_matrix_brute_force<float>(int);
  template double ** ortho_coding_matrix_brute_force<double>(int);

  template int ** orthogonal_coding_matrix<int>(int, int);
  template float ** orthogonal_coding_matrix<float>(int, int);
  template double ** orthogonal_coding_matrix<double>(int, int);

  template int ** hierarchical_nonhierarchical<int>(int);
  template float ** hierarchical_nonhierarchical<float>(int);
  template double ** hierarchical_nonhierarchical<double>(int);

  template void print_control_nonhier<int>(FILE *fs, int **, int, int, 
		  char **, cls_ta *);
  template void print_control_nonhier<float>(FILE *fs, float **, int, int, 
		  char **, cls_ta *);
  template void print_control_nonhier<double>(FILE *fs, double **, int, int, 
		  char **, cls_ta *);

}

