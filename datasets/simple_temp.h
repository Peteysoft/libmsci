#ifndef SIMPLE_TEMP_INCLUDED
#define SIMPLE_TEMP_INCLUDED

#include "peteys_tmpl_lib.h"
#include "simple_dataset.h"

#pragma interface

namespace libpetey {
namespace datasets {

//dirt simple version of the "simple dataset" class:
//any new data elements are simply inserted into their
//positions in the list after re-allocating the array...

template <class primitive>
class simple : public simple_dataset {
//  friend class simple_dataset;
//  friend class simple<float>;
//  friend class simple<double>;
//  friend class simple<long>;
//  friend class simple<time>;
//  friend class simple<string>;

  private:
    void set_type();
  protected:
    primitive *data;
    ind_type last_search;	//last index searched
  public:
    simple();
    simple(primitive *new_data, ind_type n, char sort_flag=1);
    simple(const primitive &x0, const primitive &xf, ind_type n);
    virtual ~simple();

    //copy constructor:
    simple(const simple<primitive> &other);

    template <class p2>
    simple operator =(const simple<p2> &other);

    virtual long read(FILE *fileptr);
    virtual long write(FILE *fileptr);
    virtual void print();

//    long size_of();
    ind_type add_el(primitive value);
    ind_type del(primitive value);
    ind_type search(primitive value);
    ind_type get(primitive &value, ind_type ind);
    ind_type get(primitive &value, interpol_index ind);

    interpol_index interp(primitive value);

    //comparison operator:
    template <class p2>
    int operator == (const simple<p2> &other);

};

/*typedef simple<float> simple_float;

typedef simple<double> simple_double;

typedef simple<long> simple_long;
*/

} //end namespace datasets
} //end namespace libpetey

#endif
