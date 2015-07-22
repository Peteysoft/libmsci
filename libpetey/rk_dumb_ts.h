#ifndef RK_DUMB_TS__H
#define RK_DUMB_TS__H

//"thread safe" version
namespace libpetey {

  //derivative functions are compatible with the GSL
  //(return value is ignored)

  template <class ind_type, class dep_type>
  void rk_dumb(ind_type t0, 			//initial value of independent variable
		dep_type *xinit, 		//intiall values of dependent variables
		long nx, 			//number of dependent variables
		ind_type dt, 			//time step
		long nt,			//number of time steps
		dep_type ** xvec, 		//returned integration values
		void *param,			//parameters to pass to derivative function
		int (* derivs) (ind_type, dep_type *, dep_type *, void *) );

  template <class ind_type, class dep_type>
  dep_type **rk_dumb(ind_type t0, dep_type *xinit, long nx, ind_type dt, long nt, void *param,
		int (* derivs) (ind_type, dep_type *, dep_type *, void *) );

  //these functions are mainly for testing:

  //generalized power function:
  template <class real>
  int test_rk_pfunc(real t, real *x, real *f, void *param);
    
  //indefinite integral of generalized power function:
  template <class real>
  real test_rk_ipfunc(real t, void *param);

  //test R-K dum with random one-dimensional power functions:
  template <class real>
  int test_rk_ts(int nterm,		//number of terms
			real t1,	//integration limits
			real t2,
			real h);	//step size

  //linear mapping:
  template <class real>
  int test_rk_lODE(real t, real *x, real *f, void *param);
    
  //test R-K integrator with system of linear ODEs:
  template <class real>
  int test_rk_ts2(int n,		//size of problem
			real h,		//step size
			int nt);	//number of time steps

} //end namespace libpetey

#endif

