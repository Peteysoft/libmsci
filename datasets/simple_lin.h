#include "simple_dataset.h"

#pragma interface

//dirt simple version of the "simple dataset" class:
//any new data elements are simply inserted into their
//positions in the list after re-allocating the array...

template <class primitive>
class simple : public simple_dataset {
//  friend main;
  protected:
    primitive *data;
  public:
    simple();
    //copy constructor:
    simple(const simple<primitive> &other);
    simple operator =(const simple &other);
//    long size_of();
    long add_el(primitive value);
    long del(primitive value);
    long search(primitive value);
    long get(primitive &value, long ind);
    ~simple();
};

/*typedef simple<float> simple_float;

typedef simple<double> simple_double;

typedef simple<long> simple_long;
*/
