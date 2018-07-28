#include <stdio.h>

#include "simple_temp.h"
#include "dependent_swap.h"
#include "time_class.h"

#include "composite_dataset.h"

#include "traj_int_obj.h"

struct boundary_element {
  float x;
  float y;
  short hemi;
  boundary_element *next;
};

class boundary_swap {
  private:
    traj_int_obj *tintn;
    traj_int_obj *tints;

  protected:
    boundary_element *list;		//contour
    long n;		//number of nodes

    double tind;	//current time index
    double toffs;	//difference btw. N. and S.
    double tstep;	//time step
    long nrk;		//number of runge-kutta steps for each step

    float min_spac;	//minimum spacing between nodes
    float max_spac;	//maximum spacing between nodes
    int wrap_flag;	//is it a closed contour??

    FILE *swap;		//swap file for coordinates
    char *sbase;	//base name of swap file
    long find;		//index of current file
    long bsizen;	//desired block size
    long bstart;	//start of current block
    long bsize;		//size of current block
    long nextb;		//pointer to next block

    time_class advance_current();
    void fix_current();
    void store_current();
    void get_next();
  public:
    boundary_swap();
    boundary_swap(char *sfile, char *nfile, char *bfile, long bs);
    virtual ~boundary_swap();

    void set_parm(double tcoarse, long nfine, float minspc, float maxspc);
    ind_type nt(interpol_index &offset);
    interpol_index set(time_class t0);
    time_class set(interpol_index tind1);

    interpol_index get_tind();
    time_class get_t();

    long init_circle(float x0, float y0, float r);

    time_class advance();

    long store(char *fname);		//creates a copy of the swap file
    long print(FILE *fs);
    void print_current();

};
