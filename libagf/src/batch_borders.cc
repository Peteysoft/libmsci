#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "error_codes.h"
#include "full_util.h"
#include "randomize.h"
#include "agf_lib.h"

using namespace libpetey;

namespace libagf {
  template <class real, class cls_t>
  struct bf_batch_param {
    real **x0;			//starting vector
    real **v;			//direction vector
    real **brd;			//current border guess
    general2class<real, cls_t> *twoclass;		//for calculating probabilities
    int n;			//number of elements to minimize
    int D;			//number of dimensions
    real r0;			//location of border

    //"workspace" variables:
    real **x;			//"compressed" independent (feature) variables written to input file
    real *r;			//results
  };

  template <class real, class cls_t>
  int bf_batch_param_init(bf_batch_param<real, cls_t> *bf_param, 	//parameters to pass to bfind_batch
		char *command, 			//command to calculate probabilities
		char *model, 			//model file
		int n, 				//number of variables to zero
		int D,				//number of dimensions
		int (*sample_func) (void *, real *, real *), 	//for sampling feature space
		void *sample_param,		//parameters to pass to sampling function
		real **brd,			//results of root-finding
		real *r1, 			//starting probabilities
		real *r2,			//opposite sign from r2
		real r0,			//threshold value for R
		int Mflag, 			//external command uses LIBSVM format
		int Kflag) {			//keep temporary files
    int nleft, nfilled;
    int err;

    //basic stuff:
    bf_param->twoclass=new general2class<real, cls_t>(model, command, Mflag, Kflag);
    bf_param->brd=brd;
    bf_param->D=D;
    bf_param->r0=r0;

    //allocate space for multiple parameters:
    bf_param->x0=allocate_matrix<real, int>(n, D);
    bf_param->v=allocate_matrix<real, int>(n, D);
    //bf_param->brd=allocate_matrix(real, int>(n, D);

    //allocate workspace variables:
    bf_param->x=new real *[n];
    bf_param->r=new real[n];

    //initialize the parameters for each element:
    nleft=n;
    nfilled=0;
    do {
      for (int i=nfilled; i<n; i++) {
        err=(* sample_func) (sample_param, bf_param->x0[i], bf_param->v[i]);
        if (err!=0) {
          n=i;
          break;
        }
      }
      //farm out probability estimates to external command:
      bf_param->twoclass->batchR(bf_param->x0+nfilled, r1+nfilled, nleft, D);
      bf_param->twoclass->batchR(bf_param->v+nfilled, r2+nfilled, nleft, D);

      //check that probability differences have opposite sign
      //if so, set line parameters:
      for (int i=nfilled; i<n; i++) {
        r1[i]=r1[i]-bf_param->r0;
        r2[i]=r2[i]-bf_param->r0;
        if (r1[i]*r2[i]<=0) {
          for (int j=0; j<D; j++) {
            //looks like they might interfere, but they won't:
            bf_param->v[n-nleft][j]=bf_param->v[i][j]-bf_param->x0[i][j];
            bf_param->x0[n-nleft][j]=bf_param->x0[i][j];
          }
          r1[n-nleft]=r1[i];
          r2[n-nleft]=r2[i];
          nleft--;
        }
      }
      nfilled=n-nleft;

    } while (nleft>0 && err==0);

    bf_param->n=nfilled;

    return nfilled;

    //leave the debug codes down here in case we need them later:
      for (int i=0; i<n; i++) {
        for (int j=0; j<D; j++) printf("%g ", bf_param->x0[i][j]);
        printf("%g\n", r1[i]);
        for (int j=0; j<D; j++) printf("%g ", bf_param->v[i][j]);
        printf("%g\n", r2[i]);
        printf("\n");
      }

    for (int i=0; i<nfilled; i++) {
      printf("%g %g\n", r1[i], r2[i]);
    }

  }

  template <class real, class cls_t>
  void bf_batch_param_free(bf_batch_param<real, cls_t> *parm) {
    delete [] parm->x0[0];
    delete [] parm->x0;
    delete [] parm->v[0];
    delete [] parm->v;
    delete [] parm->r;
    //delete [] parm->brd[0];
    //delete [] parm->brd;
    delete [] parm->x;
    delete parm->twoclass;
  }

  
  template <class real, class cls_t>
  void bfind_batch(real *t, 		//current value for independent variable (line parameter)
		int *ind, 		//indexes into t and r: zero only these variables
		int n, 			//current number to zero: number of values in ind
		real *r0, 		//dependent variable calculate from external ("batch") command
		void *bf_param) {	//parameters

    bf_batch_param<real, cls_t> *p2=(bf_batch_param<real, cls_t> *) bf_param;
    int k;

    //using the line parameter (independent variables), calculate the D-dimensional feature variable
    //only for elements specified by ind:
    for (int i=0; i<n; i++) {
      k=ind[i];
      for (int j=0; j<p2->D; j++) {
        p2->brd[k][j]=p2->x0[k][j]+p2->v[k][j]*t[k];
      }
      //compress samples we're still working on into a contiguous array:
      p2->x[i]=p2->brd[k];
    }

    p2->twoclass->batchR(p2->x, p2->r, n, p2->D);

    //stick the difference in probabilities in the appropriate location
    //in the results (dependent) variables:
    for (int i=0; i<n; i++) {
      k=ind[i];
      r0[k]=p2->r[i]-p2->r0;
    }
  }


  template <class real>
  int batch_rootfind(real *t1, real *t2, int n, 
		void (* func) (real *, int *, int, real *, void *), 
		real tol, int maxiter,
		void *param, 
		real *x1, real *x2) {

    int ind[n];
    real err;
    real t0[n];
    real x0[n];
    int ncur=n;
    int nnew;
    int k;

    //indexes of elements still being minimized:
    for (int i=0; i<n; i++) {
      ind[i]=i;
    }

    //repeat over maximum iterations:
    for (int i=0; i<maxiter; i++) {
      /*
      for (int j=0; j<n; j++) printf("%12g ", t1[j]);
      printf("\n");
      for (int j=0; j<n; j++) printf("%12g ", x1[j]);
      printf("\n");
      for (int j=0; j<n; j++) printf("%12g ", t2[j]);
      printf("\n");
      for (int j=0; j<n; j++) printf("%12g ", x2[j]);
      printf("\n\n");
      */

      //first order estimate--find intercept between two brackets:
      for (int j=0; j<ncur; j++) {
        real m;
        k=ind[j];
        m=(x2[k]-x1[k])/(t2[k]-t1[k]);
        t0[k]=t1[k]-x1[k]/m;			//next estimate for independent variable
      }

      //call "batch" function to get dependent variable:
      (*func) (t0, ind, ncur, x0, param);

      //see which side of the root the new estimates lie at and rebracket the root:
      for (int j=0; j<ncur; j++) {
        k=ind[j];
        if (x0[k]*x1[k]>0) {
          t1[k]=t0[k];
          x1[k]=x0[k];
        } else if (x0[k]*x2[k]>0) {
          t2[k]=t0[k];
          x2[k]=x0[k];
        } else if (x0[k]!=0) {
          fprintf(stderr, "batch_rootfind: root hasn't been bracketed properly in element %d\n", k);
          fprintf(stderr, "        f(%g)=%g; f(%g)=%g\n", t1[k], x1[k], t2[k], x2[k]);
	  assert(0);
        }
      }

      //update the indices into the list:
      nnew=0;
      for (int j=0; j<ncur; j++) {
        k=ind[j];
        //check tolerance:
        err=fabs(2*(t2[k]-t1[k])/(t1[k]+t2[k]));
        //err=fabs(x0[k]);
        if (err>=tol && x0[k]!=0) {
          ind[nnew]=k;
          nnew++;
        }
      }

      //if the tolerance has been satisfied for all the variables, we exit:
      if (nnew==0) break; else ncur=nnew;
    }

    return nnew;

    //these are debug codes:
    for (int j=0; j<n; j++) {
      printf("%12g ", t0[j]);
      //ind[j]=j;
    }
    printf("\n");
    //(*func) (t0, ind, ncur, x0, param);
    for (int j=0; j<n; j++) {
      printf("%12g ", x0[j]);
    }
    printf("\n\n");

  
  }

  template <class real, class cls_t>
  nel_ta batch_borders(char *command, char *model,
		int (*fsamp) (void *, real *, real *),
		void *s_param,
		nel_ta n, dim_ta D, 
		real tol, iter_ta maxit, real ht,
		real **border, real **gradient,
		real r0,
		int Mflag, int Kflag,
		iter_ta diter) {

    bf_batch_param<real, cls_t> bf_param;
    real t1[n], t2[n];
    real r1[n], r2[n];

    //parameters to pass to the root-finding routine:    
    n=bf_batch_param_init(&bf_param, command, model, n, D, fsamp, 
		s_param, border, r1, r2, r0, Mflag, Kflag);

    //initialize line parameters:
    for (nel_ta i=0; i<n; i++) {
      t1[i]=0;
      t2[i]=1;
    }
    //find the borders with simple root-finding:
    batch_rootfind(t1, t2, n, &bfind_batch<real, cls_t>, tol, maxit, (void *) &bf_param, r1, r2);

    //printf("batch_borders: Root-finding complete\n");

    //calculate the gradients using simple finite differencing:
    //(same basic algorith as root-finder--should make a general version of it)
    int nnew[D];
    int ncur[D];
    int flag;
    int **ind;
    int k;

    //initialize all the flags, indices and other bookkeeping:
    real hrel=ht;
    ind=new int *[n];
    ind[0]=new int[n*D];
    for (nel_ta j=0; j<D; j++) ncur[j]=n;
    for (dim_ta i=0; i<n; i++) {
      ind[i]=ind[0]+D*i;
      for (nel_ta j=0; j<D; j++) {
        ind[i][j]=i;
      }
    }

    //initialize "displaced" vectors:
    for (nel_ta i=0; i<n; i++) {
      for (dim_ta j=0; j<D; j++) {
        //take advantage of already allocated space:
        bf_param.x0[i][j]=border[i][j];
        bf_param.v[i][j]=border[i][j];
      }
    }
    //start the iteration to find non-zero numerical derivatives:
    for (int iter=0; iter<diter; iter++) {
      //printf("%d %g\n", iter, hrel);
      for (dim_ta j=0; j<D; j++) {
        //displace borders in one dimension:
        for (nel_ta i=0; i<ncur[j]; i++) {
          k=ind[i][j];
          bf_param.x0[k][j]=border[k][j]+hrel*border[k][j];
          bf_param.v[k][j]=border[k][j]-hrel*border[k][j];
          //store change in independent variable (features data):
          gradient[k][j]=bf_param.v[k][j]-bf_param.x0[k][j];
        }

        //"condense" the input data:
        for (nel_ta i=0; i<ncur[j]; i++) bf_param.x[i]=bf_param.x0[ind[i][j]];
        //and call the external probability calculator:
        bf_param.twoclass->batchR(bf_param.x, r1, ncur[j], D);
        //for each of the displacements:
        for (nel_ta i=0; i<ncur[j]; i++) bf_param.x[i]=bf_param.v[ind[i][j]];
        bf_param.twoclass->batchR(bf_param.x, r2, ncur[j], D);

        //simple finite difference:
        for (int i=0; i<ncur[j]; i++) {
          k=ind[i][j];
          //printf("%g ", r2[i]-r1[i]);
          gradient[k][j]=(r2[i]-r1[i])/gradient[k][j];
        }
        //printf("\n");

        //remove the displacement to reduce copying:
        for (nel_ta i=0; i<ncur[j]; i++) {
          k=ind[i][j];
          bf_param.x0[k][j]=border[k][j];
          bf_param.v[k][j]=border[k][j];
        }
      }
      for (dim_ta j=0; j<D; j++) nnew[j]=0;
      //recalculate all the derivatives which weren't signficant:
      for (nel_ta i=0; i<ncur[0]; i++) {
        flag=0;
        for (dim_ta j=0; j<D; j++) {
          //(at least one)/all the elements of the gradient vector have to be non-zero:
          //(only recalculate non-zero elements:)
          k=ind[i][j];
          if (gradient[k][j]==0) {
            //ind[nnew[j]][j]=k;
            //nnew[j]++;
            flag=1;
            break;
          }
        }
        if (flag) {
          for (dim_ta j=0; j<D; j++) {
            k=ind[i][j];
            ind[nnew[j]][j]=k;
            nnew[j]++;
          }
        }
      }
      flag=1;
      for (dim_ta j=0; j<D; j++) {
        if (nnew[j]!=0) flag=0;
        ncur[j]=nnew[j];
      }
      //we need to increase the relative difference for next iteration:
      hrel=hrel*2;

      //for (dim_ta j=0; j<D; j++) {
      //  for (nel_ta i=0; i<n; i++) printf("%12.5g ", gradient[i][j]);
      //  printf("\n");
      //}
      if (flag) break;
    }

    //samples with insignificant gradients must be removed:
    nnew[0]=0;
    for (nel_ta i=0; i<n; i++) {
      //flag=0;
      for (dim_ta j=0; j<D; j++) {
        if (gradient[i][j]!=0) {
          for (dim_ta j=0; j<D; j++) {
            border[nnew[0]][j]=border[i][j];
            gradient[nnew[0]][j]=gradient[i][j];
          }
          nnew[0]++;
          break;
        }
      }
    }

    //clean up:
    //printf("starting batch_borders cleanup\n");
    bf_batch_param_free(&bf_param);
    delete [] ind[0];
    delete [] ind;
    //printf("finished batch_borders cleanup\n");

    return nnew[0];
  }
        
  template nel_ta batch_borders<real_a, cls_ta>(char *, char *,
		int (*fsamp) (void *, real_a *, real_a *),
		void *, nel_ta, dim_ta, real_a, iter_ta, real_a,
		real_a **, real_a **, real_a, int, int, iter_ta);
}

