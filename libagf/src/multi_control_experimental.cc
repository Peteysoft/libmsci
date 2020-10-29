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

#include "multi_control.h"

using namespace libpetey;
using namespace std;

namespace libagf {

  //calculate the Hausdorf metric between the two partitions in a row of 
  //a coding matrix
  template <typename real, typename scalar>
  real partition_distance(real *d, int n, scalar *coding_row) {
    int n1, n2;         //number in each partition
    int ind1[n], ind2[n];               //indices of each partition
    int k;                              //1-d subscript

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
  real partition_distance(real *d,              //distance triangle (all points)
                  cls_t *cls,                   //classes
                  nel_ta n,                     //number of points
                  scalar *coding_row) {         //row of coding matrix
    int n1, n2;         //number in each partition
    int ind1[n], ind2[n];               //indices of each partition
    int k;                              //1-d subscript

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

  template int **optimal_coding_matrix<int, float>(int, int, float *);
  template int **optimal_coding_matrix<int, double>(int, int, double *);
  template float **optimal_coding_matrix<float, float>(int, int, float *);
  template float **optimal_coding_matrix<float, double>(int, int, double *);
  template double **optimal_coding_matrix<double, float>(int, int, float *);
  template double **optimal_coding_matrix<double, double>(int, int, double *);

  template int **optimal_coding_matrix<int, float, cls_ta>(float *, cls_ta *, nel_ta, int);
  template int **optimal_coding_matrix<int, double, cls_ta>(double *, cls_ta *, nel_ta, int);

  template int ** ortho_coding_matrix_nqbf<int>(int,int,int);
  template float ** ortho_coding_matrix_nqbf<float>(int,int,int);
  template double ** ortho_coding_matrix_nqbf<double>(int,int,int);

}

