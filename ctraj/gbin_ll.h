#include <stdio.h>
#include <stdint.h>

#include "linked.h"
#include "gbin_base.h"

template <class gbin_index, class coord_t, class binsub_t>
class gbin_ll:public gbin_base<gbin_index, coord_t, binsub_t> {
  protected:
    linked_list<global_bin<gbin_index, coord_t> *> **bininit;
    int update_flag;

  public:
    gbin_ll();
    gbin_ll(void *params);
    ~gbin_ll();

    int update();
    int32_t add_el(global_bin<gbin_index, coord_t> *newel); 

};

