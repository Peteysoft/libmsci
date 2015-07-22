#include <stdio.h>

#include "time_class.h"
#include "simple_temp.h"
#include "dependent_swap.h"
#include "time_class.h"

class surgery_obj {
  private:
    composite_dataset ndata;
    composite_dataset sdata;

    simple<float> *xgrid_N;
    simple<float> *ygrid_N;

    simple<float> *xgrid_S;
    simple<float> *ygrid_S;

    simple<float> *zgrid;
    simple<time_class> *tgrid;

    dependent_swap<float> *u_N;		//Northern hemisphere v field
    dependent_swap<float> *v_N;
    dependent_swap<float> *u_S;		//Southern hemisphere v field
    dependent_swap<float> *v_S;

    FILE *N_fs;				//swap file stream for N. hemisphere
    FILE *S_fs;				//swap file stream for S. hemisphere

    float **result;
  protected:
    double t0N;				//initial time index in N. hemi. data
    double t0S;				//   "      "    "    " S.  "     "
    float dt;				//course time step
    long nrk;				//number of intermediate Runge-Kutta steps
    surgery_element *list;		//first one in sorted list

  public:
    surgery_obj();
    surgery_obj(char *nfile, char *sfile);
    ~surgery_obj();

    void set(char *t0, float tcoarse, long nfine, float dmin);
    long init_circle(float x0, float y0, float r);

    time_class advance();		//advances each point by long time step
    long add_new();			//adds new points according to min. dist.
    long surgery();			//performs surgery 

    virtual long read(FILE *fs, int wrap_flag);
    virtual long write(FILE *fs);
    virtual long print(FILE *fs);

};

