#ifndef CTRAJ_VFIELD_BASE_H
#define CTRAJ_VFIELD_BASE_H 1

#include <stdio.h>
#include <stdint.h>

namespace ctraj {

  template <class real>
  class ctraj_vfield_base {
    public:
      virtual ~ctraj_vfield_base();

      //read arguments from command line and uses these to setup the object
      //returns number of time steps/grids and starting time grid
      virtual int setup(int argc, char **argv)=0;
 
      virtual void help(FILE *fs)=0;

      //convert to and from "absolute" coords: (i.e. lon-lat)
      virtual int32_t absolute(int32_t domain, 	//-1: to transformed coords
		      				//>=0: to absolute coords;
						//number determines domain
		      real *x)=0;		//coords in and out

      //convert to and from "reference" coords: (i.e. for interpolation within advected contour)
      //returns domain
      virtual int32_t reference(int32_t domain, //-1 return to "correct" domain
		      				//>=0 move to "reference" domain
		      real *x)=0;		//coords in and out

      //returns the velocity field at a given location:
      virtual int v(int32_t domain, 	//domain of point
		double tind, 			//relative time index
		real *x, 			//relative coordinates
		real *v)=0; 			//returned velocity field

      //determines whether a given location is within or without its domain
      //adjusts relative coordinates if necessary, returns new domain
      //if location is outside simulation region, returns -1 
      //or other negative error code (to be determined)
      virtual int32_t 		//returns new domain or error
		fix(int32_t domain, 	//domain of point
		double tind, 		//relative time index
		real *x)=0; 		//relative coordinates

      virtual int32_t ndim()=0;		//number of dimensions
      virtual double maxt()=0;		//maximum time index

      virtual double get_tind(char *date)=0;		//convert from date string to t-index
      virtual int get_t(double tind, char *date)=0;	//convert from t-index to date string

      //Jacobi matrix: 
      virtual int jmat(int32_t domain,
		double tind,
		real *x,
		real *jmat);

  };

  //this is the base class for velocity fields contained in a single domain,
  //i.e., defined over only one coordinate system--domain parameter is ignored
  // really all we need is the following function prototype:
  // int v(double t, float *x, float *v);
  //but it's good to encapsulate it in an object so that the parameters are
  //also encapsulated...
  //this class defines all the methods that don't really do anything...
  template <class real>
  class ctraj_vfield_single_domain:public ctraj_vfield_base<real> {
    public:
      ctraj_vfield_single_domain();

      virtual int32_t absolute(int32_t domain, real *x);
 
      virtual int32_t reference(int32_t domain, real *x);

      virtual int32_t 		//returns new domain or error
		fix(int32_t domain, 	//domain of point
		double tind, 		//relative time index
		real *x); 		//relative coordinates

      virtual double get_tind(char *date);
      virtual int get_t(double tind, char *date);

  };

} //end namespace ctraj

#endif
