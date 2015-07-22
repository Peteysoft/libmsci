#include <stdint.h>

struct lonlat_coord {
  float lon;
  float lat;
};

class azeq_binsub {
  protected:
    int32_t nx, ny;
    float sidelength;
    float xoffset, yoffset;

    void getsub(float x, float y, int32_t &xind, int32_t &yind);
    int32_t convertsub(int32_t xind, int32_t yind, short hemi);

    float xsub2x(float xind);
    float ysub2y(float yind);

  public:
    //initialize with default bins:
    azeq_binsub();
    //file contains bins:
    azeq_binsub(void *params);
    ~azeq_binsub();

    int32_t nbin();

    float metric(lonlat_coord &c1, lonlat_coord &c2);

    //get the subscript:
    int32_t getsub(lonlat_coord &coord);
    void convertsub(int32_t sub, int32_t &xind, int32_t &yind, short &hemi);
    int32_t edge_distances(lonlat_coord &coord, int32_t *sub, float *d);
    float min_edge_distance(lonlat_coord &coord);

};

