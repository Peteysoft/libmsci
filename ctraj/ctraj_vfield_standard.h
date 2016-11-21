#ifndef CTRAJ_VFIELD_STANDARD__H
#define CTRAJ_VFIELD_STANDARD__H 1

#include <stdio.h>

#include "time_class.h"

#include "simple_temp.h"
#include "dependent_swap.h"
#include "composite_dataset.h"

#include "ctraj_defaults.h"
#include "ctraj_vfield_base.h"

#include "az_eq_t.h"

namespace ctraj {
  using namespace libpetey;
  using namespace datasets;

  template <class real>
  class ctraj_vfield_standard:public ctraj_vfield_base<real> {

    protected:
      int32_t refd;

      simple<real> *z;

      composite_dataset *all[2];
      simple<time_class> *t[2];
      simple<real> *x[2];
      simple<real> *y[2];
      dependent<real> *U[2];
      dependent<real> *V[2];
      FILE *fs[2];

      double tdiff;
      double tmult;

      az_eq_t<real> *metric;

      //for inputting data:
      simple<real> *lon;           //lon. grid
      simple<real> *lat;           //lat. grid

      dependent<interpol_index> *c1[2];        //interpolation coefficients (longitude)
      dependent<interpol_index> *c2[2];        //interpolation coefficients (latitude

      void init_vvar();

      //redundant, but needed for set-up:
      int64_t page_size;

    public:
      ctraj_vfield_standard();
      virtual ~ctraj_vfield_standard();

      //initializes the object, returns 0 for success:
      int init(char **fname,              //S. hemi. file, N. hemi. file
                int64_t ps, 		     //pages size for swapping datasets
                const char *mode,       //"r" read only access, 
                                        //"r+" to add to an existing dataset
                                        //"w" to start a new one
                int32_t ref=1);         //reference domain (1 for North, 0 for South)

      //for building up datasets from scratch:
      //all methods return 0 for success

      //sets the horizontal grid (on an azimuthal-equidistant coord. system):
      int set_outgrid(real zlev,          //vertical level
                int32_t ngrid,          //number of grid-points per side
                real sl2=SIDELENGTH);   //side-length/2

      //sets the time grid:
      int set_tgrid(const simple<time_class> &t1);

      //sets the horizontal input grid:
      int set_ingrid(const simple<real> &lonin,   //longitude grid
                const simple<real> &latin);     //latitude grid

      //starts the write process:
      int write(int keep_input_grids=0);

      //returns an empty field on the lon-lat input grid:
      dependent<real> * empty_field();

      //adds a field for a single time grid:
      int add_field(int32_t tind,                         //time index
                        dependent<real> *uu,            //zonal wind field
                        dependent<real> *vv);           //meridional wind field


      virtual int setup(int argc, char **argv);

      virtual void help(FILE *fs);

      //convert to and from "absolute" coords:
      virtual int32_t absolute(int32_t domain, real *x);

      //convert to and from "reference" coords:
      virtual int32_t reference(int32_t domain, real *x);

      virtual int v(int32_t domain, double tind, real *x1, real *v);
      virtual int32_t fix(int32_t domain, double tind, real *x1);

      virtual int32_t ndim();
      virtual double maxt();

      virtual double get_tind(char *date);
      virtual int get_t(double tind, char *date);

      virtual int jmat(int32_t domain, double tind, real *loc, real *jmat);

      double get_tind(time_class date);
      time_class get_t(double tind);
      real get_lev();

      az_eq_t<real> *get_metric();

  };

} //end namespace ctraj

#endif

