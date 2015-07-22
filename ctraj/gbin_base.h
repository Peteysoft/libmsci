#include <stdio.h>
#include <stdint.h>

#include "linked.h"

template <class gbin_index, class coord_t>
struct global_bin {
  coord_t coords;
  gbin_index ind;
};

template <class gbin_index, class coord_t, class binsub_t>
class gbin_base {
  protected:
    int32_t nbin;
    binsub_t *binsub;
    global_bin<gbin_index, coord_t> ***bins;
    int32_t *nel;

    float *calc_distances(coord_t &coords, int32_t sub);
  public:
    gbin_base();
    ~gbin_base();

    void default_init(void *params);
    void clear();
    void empty();

    float nn1(coord_t &coords, global_bin<gbin_index, coord_t> *&nn);
    float nn(coord_t &coords, global_bin<gbin_index, coord_t> *&nn);
    float tn(coord_t &coords, float lat, global_bin<gbin_index, coord_t> *&tn);

    size_t save(FILE *fs);
    size_t load(FILE *fs);

};

