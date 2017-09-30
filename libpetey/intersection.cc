
#include "peteys_tmpl_lib.h"
#include "full_util.h"
#include "randomize.h"
#include "peteys_no_templates.h"

namespace libpetey {

  template <typename scalar>
  int compare_vectors(scalar *v1, scalar *v2, int n) {
    for (int i=0; i<n; i++) {
      if (v1[i]<v2[i]) return -1;
      if (v1[i]>v2[i]) return 1;
    }
    return 0;
  }

  template <typename scalar>
  int unify_vectors(scalar ***v, 		//list of lists of vectors
		  int *n, 			//number of elements in each list
		  int n0, 			//number of lists
		  int D,			//dimension of each vector
		  scalar **result,		//intersection
		  int **ind) {			//indexes into result
    scalar **all;
    long *sind;
    int k;
    long zone;
    long nt[n0+1];

    //total number of vectors:
    nt[0]=0;
    for (int i=0; i<n0; i++) nt[i+1]=nt[i]+n[i];

    //stick them all in a single array:
    //(and initialize indices while we're at it...)
    all=new scalar *[nt[n0]];
    for (int i=0; i<n0; i++) {
      for (int j=0; j<n[i]; j++) all[j+nt[i]]=v[i][j];
    }

    //sort the array:
    sind=heapsort((void **) all, nt[n0], (void *) &compare_vectors<scalar>);

    //fill result array only with non-duplicates;
    //fill index arrays with location of these
    k=0;
    for (int i=0; i<nt[n0]-1; i++) {
      zone=bin_search(nt, (long) (n0+1), sind[i], (long) -1);
      ind[zone][sind[i]-nt[zone]]=k;
      if (compare_vectors(all[sind[i]], all[sind[i+1]], D)!=0) {
        result[k]=all[sind[i]];
	k++;
      }
    }
    //last one is always in there:
    result[k]=all[sind[nt[n0]-1]];
    zone=bin_search(nt, n0+1, sind[nt[n0]-1], -1);
    ind[zone][sind[nt[n0]-1]-nt[zone]]=k;

    return k+1;

  }

  template <typename tp>
  int32_t intersection(tp *l1, int32_t n1, tp *l2, int32_t n2, tp *lr) {
    long *s1, *s2;
    int32_t i, j, k;

    s1=heapsort(l1, n1);
    s2=heapsort(l2, n2);

    i=0;
    j=0;
    k=0;
    do {
      if (l1[s1[i]] < l2[s2[i]]) {
        i++;
      } else if (l1[s1[i]] > l2[s2[i]]) {
        j++;
      } else {
        lr[k]=l1[i];
	i++;
	j++;
	k++;
      }
    } while (i<n1 && j<n2);

    return k;
  }

  //returns 0 on success:
  template <typename scalar>
  int test_unify_vectors(int nv,		//number of vectors
		  int D,			//dimension of each vector
		  int nss,			//number of sub sets
		  float frac) {			//approx. fraction for each ss
    scalar **v;			//all the vectors
    scalar **ss[nss];		//subsets
    int n[nss];			//number in each subset
    int nt=0;			//number in all subsets

    int **ind;			//list of indices
    scalar **result;		//intersection
    int nres;
    int foundflag;

    //ran_init();

    //create test matrices:
    //printf("Creating test matrices:\n");
    //printf("Filling random matrix:\n");
    v=allocate_matrix<scalar>(nv, D);
    for (int i=0; i<nv; i++) {
      for (int j=0; j<D; j++) {
        v[i][j]=ranu();
      }
    }
    //print_matrix(stdout, v, nv, D);
    //printf("\n");

    //printf("Picking subsets:\n");
    for (int i=0; i<nss; i++) {
      //printf("nv=%d, D=%d\n", nv, D);
      ss[i]=allocate_matrix<scalar>(nv, D);
      n[i]=0;
      for (int j=0; j<nv; j++) {
        //printf("[%d][%d]\n", i, j);
        if (ranu() < frac) {
          for (int k=0; k<D; k++) {
            //printf("ss[%d][%d][%d]=%g\n", i, n[i], k, v[i][k]);
            ss[i][n[i]][k]=v[j][k];
	  }
	  n[i]++;
	}
      }
      //print_matrix(stdout, ss[i], n[i], D);
      nt+=n[i];
      //printf("\n");
    }

    //ran_end();
    delete_matrix(v);

    //find the union:
    //printf("Finding the union:\n");
    ind=allocate_matrix<int>(nss, nv);
    result=new scalar *[nt];
    nres=unify_vectors(ss, n, nss, D, result, ind);

    //indices must be found in the sub-sets:
    //printf("Testing indices:\n");
    for (int i=0; i<nss; i++) {
      for (int j=0; j<n[i]; j++) {
        foundflag=0;
	//printf("ind[%d][%d]=%d\n", i, j, ind[i][j]);
        if (compare_vectors(result[ind[i][j]], ss[i][j], D)!=0) {
          fprintf(stderr, "test_intersect_vectors: incorrect index, ind[%d][%d]=%d\n", i, j, ind[i][j]);
          goto fail;
	}
      }
    }

    //result must be sorted:
    //printf("Checking ordering of results\n");
    for (int i=1; i<nres; i++) {
      if (compare_vectors(result[i-1], result[i], D)>=0) {
        fprintf(stderr, "test_intersect_vectors: result vector at %d duplicated/not sorted\n", i);
        goto fail;
      }
    }

    return 0;
    delete_matrix(ind);
    delete [] result;
    for (int i=0; i<nss; i++) delete_matrix(ss[i]);

    return 0;

    fail:
      print_matrix(stdout, result, nres, D);
      printf("\n");
      for (int i=0; i<nss; i++) {
        print_matrix(stdout, ss[i], n[i], D);
	printf("\n");
	for (int j=0; j<n[i]; j++) {
          printf("%d ", ind[i][j]);
	}
	printf("\n");
	printf("\n");
      }


      delete_matrix(ind);
      delete [] result;
      for (int i=0; i<nss; i++) delete_matrix(ss[i]);

    return -1;
  }

  template int compare_vectors<float>(float *, float *, int);
  template int compare_vectors<double>(double *, double *, int);

  template int unify_vectors<float>(float ***, int *, int, int, float **, int **);
  template int unify_vectors<double>(double ***, int *, int, int, double **, int **);

  template int32_t intersection<int32_t>(int32_t *, int32_t, int32_t *, int32_t, int32_t *);
  template int32_t intersection<int64_t>(int64_t *, int32_t, int64_t *, int32_t, int64_t *);
  template int32_t intersection<float>(float *, int32_t, float *, int32_t, float *);
  template int32_t intersection<double>(double *, int32_t, double *, int32_t, double *);

  template int test_unify_vectors<float>(int, int, int, float);
  template int test_unify_vectors<double>(int, int, int, float);

}
