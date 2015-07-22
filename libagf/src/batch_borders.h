#ifndef BATCH_BORDERS__H
#define BATCH_BORDERS__H

namespace libagf {

  template <class real>
  int batch_rootfind(real *t1, 		//left bracket of independent variable
		real *t2, 		//right brackate of independent variable
		int n, 			//number of variables to minimize
		void (* func) (real *, int *, int, real *, void *), //function to zero 
		real *tol, 		//desired tolerance
		int maxiter,		//maximum iterations
		void *param, 		//parameters to pass to func
		real *x1, 		//dependent variable at left bracket
		real *x2);		//dependent variable at right bracket

  template <class real, class cls_t>
  nel_ta batch_borders(char *command, 		//command to return probabilities
		char *model,			//model for probabilties
		int (*fsamp) (void *, real *, real *),		//command to sample from each side
		void *s_param,			//parameters to pass to fsamp
		nel_ta n, 			//desired number of border samples
		dim_ta D, 			//dimensions
		real tol, 			//tolerance of borders samples
		iter_ta maxit, 			//maximum number of iteration finding borders
		real ht,			//relative difference for gradient calculations
		real **border, 			//returned border samples
		real **gradient,		//returned gradients
		real r0=0,			//value of prob. at border
		int Mflag=0,			//use LIBSVM file format
		int Kflag=0,			//keep temporary files
		iter_ta diter=10);		//number of iterations for numerical derivatives

}

#endif

