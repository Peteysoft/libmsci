#include <stdio.h>
#include <stdint.h>

#include "linked.h"
#include "gbin_base.h"

template <class gbin_index, class coord_t, class binsub_t>
class gbin_s:public gbin_base<gbin_index, coord_t, binsub_t> {

  public:
    gbin_s();
    gbin_s(void *params);
    gbin_s(coord_t *v, int32_t n, int32_t nsamp=50);
    ~gbin_s();

    void init(coord_t *v, int32_t n);
    int32_t add_el(global_bin<gbin_index, coord_t> *newel); 

};

