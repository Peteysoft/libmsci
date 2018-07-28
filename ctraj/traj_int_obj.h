#ifndef TRAJ_INT_INCLUDED
#define TRAJ_INT_INCLUDED 1

#include <stdio.h>

#include "time_class.h"
#include "simple_temp.h"
#include "dependent_swap.h"
#include "composite_dataset.h"

class traj_int_obj {
  private:
    dependent_swap<float> *u;       //Velocity in the x direction.
    dependent_swap<float> *v;       //Velocity in the y direction.

    simple<float> *xgrid;          //Gridding in the x direction
    simple<float> *ygrid;          //Gridding in the y direction

    dataset *param[4];

    FILE *vfield_swap;
    composite_dataset vdata;
  protected:
    simple<time_class> *tgrid;        //The time gridding

    float **result;
    long res_size;
    long nres;

    interpol_index t0;
    interpol_index dt;

    time_class *tt;

    //sets the size of the "results" array:
    void setresize(long nt);

  public:

    //traj_int_obj();
    traj_int_obj(char *filename);
    virtual ~traj_int_obj();

    void set_page_size(sub_1d_type page_size);

    ind_type nt();

    virtual void integrate(float p0[2], float pf[2], double tstart, double tstep, long nt);
    virtual void integrate(float p0[2], float **result1, double tstart, double tstep, long nt);

    //integrates along pre-determined time values:
    float ** int_tarray(float p0[2], time_class t0, time_class *tarr, long nt);

    //uses a Taylor expansion to perform the same as above
    //for integrating over short distances...
    //float ** int_tarray_approx(float p0[2], time_class t0, time_class *tarr, long nt, double toffs=0.25);

    interpol_index get_tind(time_class t);
    time_class get_time(interpol_index tind);

    int get_result(time_class date, float p[2]);

    //get the Jacobi matrix of the winds (in the transformed system):
    void get_Jmatrix(float x, float y, interpol_index tind, float *j);

    //uses previously integrated trajectory to integrate
    //the "H" matrix:
    virtual void integrate_Hmatrix(float *hmatrix);
    virtual void integrate_Hcomp_matrix(float *hmatrix);

    virtual void integrate_Hmatrix(double *hmatrix);
    virtual void integrate_Hcomp_matrix(double *hmatrix);

    long print_result(FILE *fs);
    long print_result(FILE *fs, short hemi);

};

inline interpol_index traj_int_obj::get_tind(time_class t) {
  return tgrid->interp(t);
}

#endif
