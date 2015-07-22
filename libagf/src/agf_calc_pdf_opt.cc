#include <stdio.h>
#include <math.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_poly.h>
#include <gsl/gsl_linalg.h>

#include "randomize.h"
#include "peteys_tmpl_lib.h"
#include "full_util.h"
#include "roots_mins.h"

#include "agf_lib.h"

using namespace std;
using namespace libpetey;

//calculates this function:
//f(x, y) = a*ln x + b*ln^2 x + c*y^2 + d
double aofunc(double coeff[4], double x, double y, double grad[2], double hess[3]) {
  double f;
  double ex=exp(x);

  f=ex*(coeff[0] + coeff[1]*ex) + coeff[2]*y*y + coeff[3];

  //calculate the gradient:
  grad[0]=ex*(coeff[0]*ex + 2*coeff[1]*ex);
  grad[1]=2*coeff[2]*y;

  //calculates the Hessian:
  hess[0]=ex*(coeff[0]+4*coeff[1]*ex);
  hess[1]=2*coeff[2];
  hess[2]=0;

  return f;
}

//cost function:
double aocost(double **coeff, int n, double x, double y, double grad[2], double hess[3]) {
  double c=0;
  grad[0]=0;
  grad[1]=0;
  hess[0]=0;
  hess[1]=0;
  hess[2]=0;
  for (int i=0; i<n; i++) {
    double f;
    double g2[2];
    double h2[3];

    //cost function:
    f=aofunc(coeff[i], x, y, g2, h2);
    c+=f*f;

    //gradient of cost function:
    grad[0]+=f*g2[0];
    grad[1]+=f*g2[1];

    //hessian of cost function:
    hess[0]+=g2[0]*g2[0]+f*h2[0];
    hess[1]+=g2[1]*g2[1]+f*h2[1];
    hess[2]+=g2[0]*g2[1]+f*h2[2];
  }

  grad[0]*=2;
  grad[1]*=2;
  hess[0]*=2;
  hess[1]*=2;
  hess[2]*=2;

  return c;
}

//simple Newton's method:
int solve_agf_opt(double **coeff, int n, double &p, double &pp) {
  double c, cold;
  double grad[2];
  double hess[3];
  double pn, ppn;
  double de;
  double alpha=0.1;
  double tol=1e-13;

  p=log(p);
  c=aocost(coeff, n, p, pp, grad, hess);
  do {
    cold=c;
    de=hess[0]*hess[1]-hess[2]*hess[2];
    pn=p-(hess[1]*grad[0]-hess[2]*grad[1])/de;
    ppn=pp-(hess[0]*grad[1]-hess[2]*grad[0])/de;
    c=aocost(coeff, n, pn, ppn, grad, hess);
    printf("c=%g; p=%g; ddp=%g\n", c, exp(p), pp);
    //printf("grad=[%g, %g]; hess=[%g, %g, %g]\n", grad[0], grad[1], 
//		hess[0], hess[1], hess[2]);

    alpha=0.1;
    while (c>cold) {
      pn=p-alpha*grad[0];
      ppn=pp-alpha*grad[1];
      c=aocost(coeff, n, pn, ppn, grad, hess);
      alpha=alpha/10;
      printf("alpha=%g; c=%g; p=%g; ddp=%g\n", alpha, c, exp(pn), ppn);
    }

    p=pn;
    pp=ppn;
  } while (fabs(grad[0])>tol && fabs(grad[1])>tol);

  p=exp(p);

  return 0;
}


namespace libagf {

//option summary:
//opts=??0	parm[1]=initial filter variance, half it each subsequent trial
//
//opts=??1	parm[0]=lower filter variance
//		parm[1]=upper filter variance
//
//opts=??2	parm[0]=lower filter variance bracket when solving with fixed W
//		parm[1]=upper filter variance bracket when solving with fixed W
//		parm[2]=minimum value for W
//		parm[3]=maximum value for W
//		geometric progression for filter variance
//opts=?1?	parm same as above; geometric progression, random filter variance
//opts=1??		"	; geometric progression, randomize sample order

//Guassian filter with optimal error rate:
template <class real>
real agf_calc_pdf_opt(real **mat,		//matrix of samples
		dim_ta D,		//dimensionality of problem
		nel_ta n,           //number of samples
		real *test,		//test point
		real *parm,			//initial filter variance/range of W
		nel_ta ntest,			//number of test cases
		real &err,		//returned error estimate
		int opts)			//options
{

  //since the matrix has a very large condition number,
  //we do all calculations in double precision:
  double *r2;		//distances squared
  real *r2f;
  double *w;		//weights
  real *wdum;		//dummy weights (aren't used...)
	
  double *pdf;		//array of estimated pdfs
  double *var;		//array of filter variances
  real var0, varf;
  double dvar;		//for calculating variances
	
  double *tw;                     //current value of the total weight
  double n1;			//for calculating the norm
  double norm;			//normalization coefficient

  //estimate of the pdf:
  double p;
  //Laplacian:
  double dp;
  double var_opt[3], var_opt1;
  double eb, ev;		//bias term, variance term
  double et;			//estimated error
  double ptest;			//for checking the fit
  double tw_test;

  long *ind;

  int varvar=opts%10;
  int randomflag=(opts/10)%10;
  int randomizeflag=(opts/100)%10;

  r2=new double[n];
  w=new double[n];

  tw=new double[ntest];
  pdf=new double[ntest];
  var=new double[ntest];

  //calculate distances: (do all calculations double precision...)

  //calculate the weights:
  printf("%d %d %d\n", opts, randomflag, randomizeflag);
  printf("Calculating weights and distances\n");
  //printf("w=[");
  for (nel_ta i=0; i<n; i++) {
    r2[i]=metric2(mat[i], test, D);
  }

  if (varvar==0) {
    var[0]=parm[0];
    for (nel_ta j=1; j<ntest; j++) var[j]=var[j-1]/2;
  } else {
    if (varvar==1) {
      var0=parm[0];
      varf=parm[1];
    } else {
      //find range of filter variance from range of W:
      r2f=new real[n];
      wdum=new real[n];
      for (nel_ta j=0; j<n; j++) r2f[j]=r2[j];
      AGF_CALC_W_FUNC(r2f, n, parm[2], parm, wdum, var0);
      AGF_CALC_W_FUNC(r2f, n, parm[3], parm, wdum, varf);
      delete [] r2f;
      delete [] wdum;
      //printf("var0=%g; varf=%g\n", var0, varf);
    }
    dvar=pow((double) varf/(double) var0, 1./(ntest-1.));

    //use geometric progression, randomize samples, use random filter variance...
    if (randomflag==0) {
      var[0]=var0;
      for (nel_ta j=1; j<ntest; j++) var[j]=var[j-1]*dvar;

      //randomize samples:
      if (randomizeflag) {
        ind=randomize(ntest);
	for (nel_ta j=0; j<ntest; j++) printf("%ld ", ind[j]);
	printf("\n");
	map_vector_inplace(var, ind, ntest);
	delete [] ind;
      }

    //printf("W=%g\n", tw);
    } else {
      for (nel_ta j=0; j<ntest; j++) {
        var[j]=var0*pow(dvar, ntest*ranu());
        //var[j]=gsl_rng_uniform(rann)*(var0-varf)+varf;
        if (randomizeflag) {
          ind=heapsort(var, ntest);
          map_vector_inplace(var, ind, ntest);
          delete [] ind;
        }
      }
    }
  }

  tw[0]=0;
  for (nel_ta i=0; i<n; i++) {  
    w[i]=exp(-r2[i]/var[0]/2);
    //printf("%g ", w[i]);
    tw[0]+=w[i];
  }
  printf("\n\n");

  //repeatedly square the weights and record both the filter width and the result:
  //printf("W=%g\n", tw);

  printf("Calculating pdf estimates\n");
  //calculate pdf:
  n1=sqrt(var[0]*M_PI*2);
  norm=1;
  for (dim_ta j=0; j<D; j++) norm*=n1;
  pdf[0]=tw[0]/norm;

  for (nel_ta i=1; i<ntest; i++) {
    tw[i]=0;
    tw_test=0;
    //printf("w=[");
    for (nel_ta j=0; j<n; j++) {
      //lovely hack:
      if (varvar==0) w[j]*=w[j]; else w[j]=exp(-r2[j]/var[i]/2);
      tw[i]+=w[j];
      //printf("%g ", w[i]);
    }
    n1=sqrt(var[i]*M_PI*2);
    norm=1;
    for (dim_ta j=0; j<D; j++) norm*=n1;
    pdf[i]=tw[i]/norm;
  }

  //estimate the true pdf using a least-squares:
  //once we collect 
  //coefficients for a set of equations of the form:
  //a*P + b*P^2 + c*(grad^2)^2 + d = 0
  //we have everything we need to know about the problem:
  //pass the coefficients off another routine to solve them...
  double **coeff;

  coeff=new double*[ntest];
  coeff[0]=new double[ntest*4];
  for (int i=0; i<ntest; i++) {
    //coefficients for the basic equation:
    coeff[i]=coeff[0]+i*4;
    coeff[i][0]=-2*n*pdf[i];
    coeff[i][1]=n*n-1/pow(var[i], D);
    coeff[i][2]=-var[i]*var[i];
    coeff[i][3]=pdf[i]*pdf[i];

    //printf("coefficients: A=%g, B=%g, C=%g, D=%g\n", 
//		coeff[i][0], coeff[i][1], coeff[i][2], coeff[i][3]);
    //printf("n^2 g_%d=%g P + %g P^2 + %g (\\nabla^2)^2 + %g \\\\\n", i+1,
//		coeff[i][0], coeff[i][1], coeff[i][2], coeff[i][3]);


  }

  p=agf_calc_pdf_opt2(mat, D, n, test, var, ntest, err);

  //p=1;
  //dp=1;
  //solve_agf_opt(coeff, ntest, p, dp);
    
  return p;

  if (opts==2) {
    ind=heapsort(var, ntest);
  } else {
    ind=new long[ntest];
    for (long i=0; i<ntest; i++) ind[i]=i;
  }
  /*
  //for (int i=0; i<ngood; i++) {
    dp1=(-p1*(b1+p1*c1)-d1)/a1;
    dp[i]=
    //calculate optimal filter width:
    var_opt[i]=pow(D*p1/p2/2, 1./(D+2.));
  //use this to calculate the pdf:
  tw_test=0;
  for (nel_ta i=0; i<n; i++) {
    w[i]=exp(-r2[i]/var_opt1/2);
    tw_test+=w[i];
  }
  n1=sqrt(var_opt1*M_PI*2);
  norm=1;
  for (dim_ta j=0; j<D; j++) norm*=n1;

  //calculate error:
  ev=p1/pow(var_opt1, D)/n/n;
  eb=p2*var_opt1*var_opt1/n/n;

  et=ev+eb;
  err=sqrt(et);

  printf("Errors: bias=%g; variance=%g; total=%g\n", eb, ev, et);

  p3=tw_test/norm/n;

  printf("2 pdf estimates; p^2: %g %g %g\n", sqrt(p1), p3, sqrt(p2));
  printf("optimal band width: %g\n", var_opt1);


  //for testing:
  if (opts==2) {
    ind=heapsort(var, ntest);
    for (nel_ta i=0; i<ntest; i++) {
      ptest=sqrt(p1/n/n/pow(var[ind[i]], D)+p2*var[ind[i]]/n/n);
      printf("%g %g %g %g\n", var[ind[i]], pdf[ind[i]]/n, tw[ind[i]], ptest);
    }
    delete [] ind;
  } else {
    for (nel_ta i=0; i<ntest; i++) {
      ptest=sqrt(p1/n/n/pow(var[i], D)+p2*var[i]/n/n);
      printf("%g %g %g %g\n", var[i], pdf[i]/n, tw[i], ptest);
    }
  }

  //clean up:
  printf("Cleaning up...\n");
  */
  
  delete [] r2;
  delete [] w;
  delete [] var;
  delete [] pdf;
  delete [] tw;

  /*if (opts&1) {
    delete [] d2f;
    delete [] wdum;
  }*/

  //return p3;

}

//calculate the coefficients:
//c[0]=-2f(h)
//c[1]=1/(n^2*h^(2*D)
//c[2]=-h^2/n^2
//c[3]=f^2(h)
template <class real>
int agf_opt_calc_coeff(real *d2, 	//distances squared
		nel_ta n, 		//number of samples
		dim_ta D, 		//number of dimensions
		real h2, 		//squared filter width
		real *coeff) {		//returned coefficients
  real f;		//estimated pdf at bandwidth sqrt(h2)
  real w;		//individual weight
  real n1;
  real norm;		//normalization coefficient
  real h2D;		//h^(2D)

  f=0;
  for (nel_ta i=0; i<n; i++) {
    w=exp(-d2[i]/h2/2);
    f+=w;
  }

  n1=sqrt(h2*M_PI*2);
  norm=n1;
  h2D=h2;
  for (dim_ta j=1; j<D; j++) {
    norm*=n1;
    h2D*=h2;
  }
  //f=f/norm/n;
  f=f/norm;

  coeff[0]=-2*n*f;
  coeff[1]=-1/h2D;
  coeff[2]=-h2*h2;
  coeff[3]=f*f;

}

template <class real>
void ** agf_opt_scale_coefficient_matrix(real **c0, 	//initial coef. matrix
			int ntrial,			//number of rows
			dim_ta D,			//number of dimensions
			real sf2, 			//scale factor squared
			real **c, 			//transformed matrix
			nel_ta n=1) {			//number of samples
  real sf_2D;			//scale factor^(2D)
  real sf_D;			//scale factor^D
  real sf4=sf2*sf2;		//scale factor^4

  sf_2D=pow(sf2, D);
  sf_D=sqrt(sf_2D);

  //fill new coefficients, transformed by scale factor:
  for (int i=0; i<ntrial; i++) {
    c[i][0]=c0[i][0]/sf_D;
    c[i][1]=n*n+c0[i][1]/sf_2D;		//add the one
    c[i][2]=c0[i][2]*sf4;
    c[i][3]=c0[i][3]/sf_2D;
    //printf("%g %g %g %g\n", gsl_matrix_get(c, i, 0), c0[i][1]/sf_2D, gsl_matrix_get(c, i, 2), gsl_matrix_get(c, i, 3));
  }
}

template <class real>
void agf_opt_cond_num_from_scale_factor(real sf2, 	//scale factor squared
		void * param, 			//original coefs., #
		real *cond) {			//return condition number
  void **p2=(void **) param;
  real **c0;		//original coefficients (scale factor 1)
  real **c;		//transformed coefficients
  int ntrial;		//number of sets of coefficients
  nel_ta n;		//number of samples
  dim_ta D;		//number of dimensions
  gsl_matrix *u;	//transformed coefficients
  gsl_matrix *vt;	//singular vectors
  gsl_vector *s;	//singular values
  gsl_vector *work;

  c0=(real **) p2[0];		//unscaled coefficient matrix
  ntrial=*(int *) p2[1];	//number of rows
  D=*(dim_ta *) p2[2];		//dimensionality of problem
  n=*(nel_ta *) p2[3];		//number of samples

  //fill new coefficients, transformed by scale factor:
  c=allocate_matrix<real, int>(ntrial, 4);
  agf_opt_scale_coefficient_matrix(c0, ntrial, D, sf2, c, n);
  //transfer to gsl-compatible matrix:
  u=gsl_matrix_alloc(ntrial, 4);
  for (int i=0; i<ntrial; i++) {
    for (int j=0; j<4; j++) gsl_matrix_set(u, i, j, c[i][j]);
  }

  //find the singular values:
  vt=gsl_matrix_alloc(4, 4);
  s=gsl_vector_alloc(4);
  work=gsl_vector_alloc(4);
  gsl_linalg_SV_decomp(u, vt, s, work);

  //calculate the condition number from the singular values:
  *cond=gsl_vector_get(s, 0)/gsl_vector_get(s, 3);
  //printf("condition number=%g\n\n", *cond);

  gsl_matrix_free(u);
  gsl_matrix_free(vt);
  gsl_vector_free(s);
  gsl_vector_free(work);

  delete_matrix(c);
}

//Guassian filter with optimal error rate:
template <class real>
real agf_calc_pdf_opt2(real **mat,		//matrix of samples
		dim_ta D,		//dimensionality of problem
		nel_ta n,           //number of samples
		real *test,		//test point
		double *var,		//list of variance (bandwidth squared)
		int ntrial,		//number of trials to use
		real &err)		//returned error estimate
{
  //while the inputs might be single precision, we do all the calculations
  //double:
  double r2[n];		//distance squared
  double r22[n];	//distance squared and scaled
  double **c0;		//initial coefficients
  double **cscaled;	//coefficients from scaled distances
  double **c;		//coefficients corrected by scale factor
  double sf2;		//scale factor squared
  double sf2a, sf2b;	//brackets
  double cnum;		//condition number
  double cnum_a, cnum_b;	//condition number at brackets
  long errstatus;	//error status from minimization routine
  void *param[4];	//parameters
  double p, dp;		//the stuff we want

  for (nel_ta i=0; i<n; i++) r2[i]=metric2(mat[i], test, D);

  c0=allocate_matrix<double, int>(ntrial, 4);

  for (int i=0; i<ntrial; i++) agf_opt_calc_coeff(r2, n, D, var[i], c0[i]);

  printf("coefficient matrix:\n");
  print_matrix(stdout, c0, ntrial, 4);
  printf("\n");

  //fill parameters for optimization routines:
  param[0]=c0;
  param[1]=&ntrial;
  param[2]=&D;
  param[3]=&n;
  sf2=1;
  agf_opt_cond_num_from_scale_factor(sf2, param, &cnum);
  sf2a=sf2;
  sf2b=sf2;
  do {
    sf2a/=10;
    agf_opt_cond_num_from_scale_factor(sf2a, param, &cnum_a);
    printf("%g %g\n", sf2a, cnum_a);
    if (cnum_a < cnum) {
      sf2b=sf2;			//we have all three
      cnum_b=cnum;
      sf2=sf2a;
      cnum=cnum_a;
    }
  } while (cnum_a<=cnum);
  if (sf2==sf2b) {
    do {
      sf2b*=10;
      agf_opt_cond_num_from_scale_factor(sf2b, param, &cnum_b);
      printf("%g %g\n", sf2b, cnum_b);
      if (cnum_b < cnum) {
        sf2a=sf2;			//we have all three
        cnum_a=cnum;
        sf2=sf2b;
        cnum=cnum_b;
      }
    } while (cnum_b <= cnum);
  }

  double sf3=1e-10;
  for (int i=0; i<21; i++) {
    double cn;
    agf_opt_cond_num_from_scale_factor(sf3, param, &cn);
    printf("(2) %g %g\n", sf3, cn);
    sf3*=10;
  }

  printf("minimization brackets: %g %g %g\n", sf2a, sf2, sf2b);
  printf("minimization brackets: %g %g %g\n", cnum_a, cnum, cnum_b);

  //minimize the condition number of the coefficient matrix wrt scale factor:
  sf2=min_golden<double>(&agf_opt_cond_num_from_scale_factor<double>,
		  param, sf2a, sf2, sf2b, 0.001, 1000, errstatus, cnum);

  printf("optimal scale-factor=%g\n", sf2);
  agf_opt_cond_num_from_scale_factor(sf2, param, &cnum);
  printf("optimal condition number=%g\n", cnum);

  //transform matrix based on this factor:
  c=allocate_matrix<double, int>(ntrial, 4);
  //sf2=1;
  agf_opt_scale_coefficient_matrix(c0, ntrial, D, sf2, c, n);

  printf("optimized coefficient matrix:\n");
  print_matrix(stdout, c, ntrial, 4);
  printf("\n");

  //solve for the unknowns:
  p=1;
  dp=1;
  solve_agf_opt(c, ntrial, p, dp);

  p=exp(p)/pow(sf2, D)/n;

  delete_matrix(c0);
  delete_matrix(c);

  return p;

  print_matrix(stdout, c0, ntrial, 4);
  printf("\n");

  //try different scale factors to make sure coefficients are being
  //scaled correctly:
  cscaled=allocate_matrix<double, int>(ntrial, 4);
  sf2=0.1;
  agf_opt_cond_num_from_scale_factor(sf2, param, &cnum);
  for (nel_ta i=0; i<n; i++) r22[i]=sf2*r2[i];
  for (int i=0; i<ntrial; i++) agf_opt_calc_coeff(r22, n, D, sf2*var[i], cscaled[i]);
  print_matrix(stdout, cscaled, ntrial, 4);
  printf("\n");
  sf2=10;
  agf_opt_cond_num_from_scale_factor(sf2, param, &cnum);
  for (nel_ta i=0; i<n; i++) r22[i]=sf2*r2[i];
  for (int i=0; i<ntrial; i++) agf_opt_calc_coeff(r22, n, D, sf2*var[i], cscaled[i]);
  print_matrix(stdout, cscaled, ntrial, 4);
  printf("\n");

  delete_matrix(cscaled);

  return 0;

}


template float agf_calc_pdf_opt<float>(float **mat, dim_ta D, nel_ta n, 
		float *test, float *parm, nel_ta ntest, float &err, int opts);
template double agf_calc_pdf_opt<double>(double **mat, dim_ta D, nel_ta n, 
		double *test, double *parm, nel_ta ntest, double &err, int opts);

template float agf_calc_pdf_opt2<float>(float **mat, dim_ta D, nel_ta n, 
		float *test, double *parm, int ntest, float &err);
template double agf_calc_pdf_opt2<double>(double **mat, dim_ta D, nel_ta n, 
		double *test, double *parm, int ntest, double &err);

}

