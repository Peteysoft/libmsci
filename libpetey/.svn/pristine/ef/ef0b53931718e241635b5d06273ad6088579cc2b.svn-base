#include <math.h>
#include <gsl/gsl_linalg.h>
#include "solve_lode.h"
#include "randomize.h"
#include "gsl_util.h"

//"thread safe" version
namespace libpetey {

  //derivative functions are compatible with the GSL
  //(return value is ignored)

  template <class ind_type, class dep_type>
  void rk_dumb(ind_type t0, dep_type *xinit, long nx, ind_type dt, long nt,
		dep_type ** xvec, void *param,
		 int (* derivs) (ind_type, dep_type *, dep_type *, void *) ) {
    ind_type t, tph;
    dep_type dxall[nx], dx[nx];
    dep_type xint[nx];
    long i, j;

    for (j=0;j<nx;j++) xvec[0][j]=xinit[j];

    t=t0;
    for (i=1;i<=nt;i++) {
      //take the first set of derivatives:
      (* derivs) (t, xvec[i-1], dx, param);
      for (j=0;j<nx;j++) {
        xint[j]=xvec[i-1][j]+dx[j]*dt/2;	//calculate intermediate dep. values
        dxall[j]=dx[j]/6;				//start calculating the total deriv
      }

      //take the second set of derivatives:
      tph=dt/2.;
      tph=t+tph;
      (* derivs) (tph, xint, dx, param);
      for (j=0;j<nx;j++) {
        xint[j]=xvec[i-1][j]+dx[j]*dt/2;
        dxall[j]=dxall[j]+dx[j]/3;
      }

      //take the third set of derivatives:
      (* derivs) (tph, xint, dx, param);
      for (j=0;j<nx;j++) {
        xint[j]=xvec[i-1][j]+dx[j]*dt;
        dxall[j]=dxall[j]+dx[j]/3;
      }

      //take the fourth set of derivatives:
      t=t+dt;
      (* derivs) (t, xint, dx, param);
      for (j=0;j<nx;j++) {
        xvec[i][j]=xvec[i-1][j]+(dxall[j]+dx[j]/6)*dt;	//calculate next vector
      }
    }
  }

  template <class ind_type, class dep_type>
  dep_type **rk_dumb(ind_type t0, dep_type *xinit, long nx, ind_type dt, long nt, void *param,
		 int (* derivs) (ind_type, dep_type *, dep_type *, void *) ) {
    dep_type **xvec;
    xvec=new dep_type * [nt+1];
    xvec[0]=new dep_type [nx*(nt+1)];
    for (int i=0; i<=nt; i++) xvec[i]=xvec[0]+i*nx;
    rk_dumb(t0, xinit, nx, dt, nt, xvec, param, derivs);
    return xvec;
  }

  //generalized power function:
  template <class real>
  int test_rk_pfunc(real t, real *x, real *f, void *param) {
    int err=0;
    void **p2=(void **) param;
    int nterm=*(int *) p2[0];		//number of terms
    real *c=(real *) p2[1];		//coefficients
    real t_n;				//t^n
    real fact;				//n!

    t_n=1;
    fact=1;
    *f=c[0];
    for (int i=1; i<nterm; i++) {
      t_n*=t;
      fact*=i;
      *f+=c[i]*t_n/fact;
    }
    return err;
  }
    
  //indefinite integral of generalized power function:
  template <class real>
  real test_rk_ipfunc(real t, void *param) {
    void **p2=(void **) param;
    int nterm=*(int *) p2[0];		//number of terms
    real *c=(real *) p2[1];		//coefficients
    real t_n;				//t^n
    real fact;				//n!
    real f;

    t_n=1;
    fact=1;
    f=0;
    for (int i=1; i<=nterm; i++) {
      t_n*=t;
      fact*=i;
      f+=c[i-1]*t_n/fact;
    }
    return f;
  }
    

  template <class real>
  int test_rk_ts(int nterm,		//number of terms
			real t1,	//integration limits
			real t2,
			real h) {	//step size
    void *param[2];
    real c[nterm];
    real x0=0;
    long nt=(t2-t1)/h;
    real **x;
    real x1, x2;
    real xanal;
    real k;
    real errp, erra;		//predicted and actual error
    int exit_code;

    x=new real *[nt+2];
    x[0]=new real[nt+2];
    for (int i=0; i<nt+2; i++) x[i]=x[0]+i;

    //random coefficients:
    for (int i=0; i<nterm; i++) c[i]=ranu();

    //initialize parameters:
    param[0]=&nterm;
    param[1]=c;

    //given step size:
    rk_dumb(t1, &x0, 1, h, nt, x, param, &test_rk_pfunc<real>);
    //initial condition gets copied to itself (***):
    rk_dumb(h*nt, x[nt], 1, t2-h*nt, 1, x+nt, param, &test_rk_pfunc<real>);
    x1=x[nt+1][0];

    //double step size:
    nt=(t2-t1)/h/2;
    rk_dumb(t1, &x0, 1, 2*h, nt, x, param, &test_rk_pfunc<real>);
    //initial condition gets copied to itself (***):
    rk_dumb(2*nt*h, x[nt], 1, t2-2*h*nt, 1, x+nt, param, &test_rk_pfunc<real>);
    x2=x[nt+1][0];

    //analytic result:
    xanal=test_rk_ipfunc(t2, param)-test_rk_ipfunc(t1, param);

    //predicted error:
    k=(x1-x2)/pow(h, 5)/15;

    errp=(x1-x2)/15;
    erra=xanal-x1;

    printf("R-K integrator:\n");
    printf("numerical   = %g  +   %g\n", x1, (x1-x2)/15);
    printf("analytic    = %g (dif=%g)\n", xanal, xanal-x1);

    //if error is the right order-of magnitude, then we say it passed:
    if (fabs(2*(errp-erra))/(errp+erra)<1) exit_code=0; else exit_code=1;

    delete [] x[0];
    delete [] x;

    return exit_code;
  }

  //linear mapping:
  template <class real>
  int test_rk_lODE(real t, real *x, real *f, void *param) {
    int err=0;
    gsl_matrix *A=(gsl_matrix *) param;

    //for (int i=0; i<A->size1; i++) printf("%g ", x[i]);
    //printf("\n");

    for (int i=0; i<A->size1; i++) {
      f[i]=0;
      for (int j=0; j<A->size2; j++) {
        f[i]+=gsl_matrix_get(A, i, j)*x[j];
        //printf("%g ", gsl_matrix_get(A, i, j));
      }
      //printf("\n");
      //printf("%g ", f[i]);
    }
    //printf("\n");
    return err;
  }
    
  //maybe overkill, but another test routine:
  template <class real>
  int test_rk_ts2(int n,		//size of problem
			real h,		//step size
			int nt) {	//number of time steps
    gsl_matrix *A;			//linear mapping
    real x0[n];
    real **x;				//integration results
    real *x1, *x2;
    real xanal[n];
    real errp[n], erra[n];		//predicted and actual error
    int bind;
    int exit_code=0;

    x=new real *[nt+1];
    x[0]=new real[(nt+1)*n];
    for (int i=1; i<nt+1; i++) x[i]=x[0]+i*n;

    //random coefficients:
    A=random_gsl_matrix(n, n);

    //start x as a basis vector:
    bind=(int) (n*ranu());
    printf("bind=%d\n", bind);
    for (int i=0; i<n; i++) x0[i]=0;
    x0[bind]=1;

    //given step size:
    rk_dumb((real) 0., x0, (long) n, h, (long) nt, x, (void *) A, &test_rk_lODE<real>);
    x1=x[nt];

    //double step size:
    rk_dumb((real) 0., (real *) x0, (long) n, 2*h, (long) nt/2, x, (void *) A, &test_rk_lODE<real>);

    if (nt % 2 != 0 ) {
      rk_dumb(2*(nt/2)*h, x[nt/2], (long) n, h, (long) 1, x+nt/2, (void *) A, &test_rk_lODE<real>);
    }
    x2=x[(nt+1)/2];	//works if n>1

    //analytic result:
    gsl_matrix_complex *v=gsl_matrix_complex_alloc(n, n);
    gsl_matrix_complex *vinv=gsl_matrix_complex_alloc(n, n);
    gsl_vector_complex *eval=gsl_vector_complex_alloc(n);
    gsl_vector *gx0=gsl_vector_alloc(n);
    gsl_vector *gx=gsl_vector_alloc(n);

    gsl_vector_set_zero(gx0);
    gsl_vector_set(gx0, bind, 1);

    //perform an eigenvalue decomposition:
    diagonalize_matrix(A, v, eval, vinv);
    solve_lode(gx0, v, eval, vinv, h*nt, gx);

    gsl_matrix_complex_free(v);
    gsl_matrix_complex_free(vinv);
    gsl_vector_complex_free(eval);

    real err[n];		//truncation error
    real res[n];		//residual
    
    //calculate truncation error and residual:
    for (int i=0; i<n; i++) err[i]=(x1[i]-x2[i])/15;
    for (int i=0; i<n; i++) res[i]=(gsl_vector_get(gx, i)-x1[i]);

    printf("R-K integrator:\n");
    printf("num  ");
    for (int i=0; i<n; i++) printf(" %g", x1[i]);
    printf("\n");
    printf("anal ");
    for (int i=0; i<n; i++) printf(" %g", gsl_vector_get(gx, i));
    printf("\n");
    printf("err  ");
    for (int i=0; i<n; i++) printf(" %g", err[i]);
    printf("\n");
    printf("res  ");
    for (int i=0; i<n; i++) printf(" %g", res[i]);
    printf("\n");

    //if residual is the right order-of magnitude, then we say it passed:
    for (int i=0; i<n; i++) {
      if (fabs(2*(res[i]-err[i]))/(res[i]+err[i])>1) exit_code=1;
    }

    gsl_vector_free(gx);
    gsl_vector_free(gx0);

    gsl_matrix_free(A);

    delete [] x[0];
    delete [] x;

    return exit_code;
  }

  template int test_rk_ts<float>(int, float, float, float);
  template int test_rk_ts<double>(int, double, double, double);

  template int test_rk_ts2<float>(int, float, int);
  template int test_rk_ts2<double>(int, double, int);

} //end namespace libpetey

