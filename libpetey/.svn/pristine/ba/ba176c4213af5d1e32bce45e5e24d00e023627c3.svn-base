#ifndef COMPOSITE_DEFINED

#include <stdio.h>

#include "string_petey.h"
#include "dataset.h"

#define DELETE_VARS 0
#define NO_DELETE_VARS 1

namespace libpetey {
namespace datasets {

class composite_dataset:public dataset {
  protected:
    long n_var;			//number of variables
    string **symbol;		//list of symbols
    dataset **variable;		//array of variables
    long *index;		//index sorting the symbols
    long table_size;		//actual size of the table (as opposed to that used)
    long *ndep;		//number of dependents
    long **dependencies;	//index of dependencies
    int update_flag;
    short nodeleteflag;
  public:
    composite_dataset();
    composite_dataset(short nodelete);
     ~composite_dataset();

    void update();		//create sorted index

    //read or write to a binary file:
    virtual long read(FILE *fileptr);
    virtual long write(FILE *fileptr);

    //add a new variable:
    long add_var(string name);

    //change a variable at a given location:
    long cvar(long loc, dataset *new_var);

    //delete a variable:
    long del_var(string name);

    //search for a variable:
    long search_var(string name, long &loc);

    //get a variable:
    dataset * get_var(long loc);

    //print to std. out:
    void print();

};

} //end namespace datasets
} //end namespace libpetey

#define COMPOSITE_DEFINED
#endif
