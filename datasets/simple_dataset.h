//
// This software is released under the following terms:
//
// 1. No commercial use.
// 2. Copies and derivative works are free to use and modify.
// 3. Attribution must be given to all contributors of both original and derivative works.
//
// Authors:
//
// 2017-07-16 Peter Mills: added license information 
//

//
// Defines the base virtual class of the family of "simple dataset" classes.
//

#ifndef SIMPLE_DATASET_INCLUDED
#define SIMPLE_DATASET_INCLUDED

#include "dataset.h"
namespace libpetey {
  namespace datasets {
    class simple_dataset;
  }
}
#include "dependent_dataset.h"

namespace libpetey {
namespace datasets {

//class dependent_dataset_abstract;
//class simple_dataset<dependent_dataset_abstract *>;

class simple_dataset: public dataset {

  protected:
    dependent_dataset **dependents;		//the dependents
    rank_type *rank;				//the ranks of each of the dependents
    long ndep;				//the number of dependents
  public:
    simple_dataset();
    virtual ~simple_dataset();

    virtual long read(FILE *fileptr)=0;
    virtual long write(FILE *fileptr)=0;

    virtual void print()=0;

    long add_dependent(dependent_dataset *dep, rank_type r);	//adds a dependent
    long del_dependent(dependent_dataset *dep);			//removes a dependent
};

} //end namespace datasets
} //end namespace libpetey

#endif

