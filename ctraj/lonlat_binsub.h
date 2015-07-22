#include <stdint.h>

struct lonlat_coord {
  float lon;
  float lat;
};

class lonlat_binsub {
  protected:
    float offset;		//size of polar bin
    int32_t nlat;
    int32_t *nlon;
    int32_t *nlon_cum;

    void cumulate_nlon();
    void getsub(float lon, float lat, int32_t &lonind, int32_t &latind);
    int32_t convert_sub(int32_t lonind, int32_t latind);

    void lonlat_from_sub(float lonind, float latind, float &lon, float &lat);
    float lat_from_sub(float latind);
  public:
    //initialize with default bins:
    lonlat_binsub();
    //file contains bins:
    lonlat_binsub(void *params);
    ~lonlat_binsub();

    int32_t nbin();

    float metric(const lonlat_coord &c1, const lonlat_coord &c2);

    //get the subscript:
    int32_t getsub(const lonlat_coord &coord);
    int32_t edge_distances(const lonlat_coord &coord, int32_t *sub, float *d);

};

inline void lonlat_binsub::lonlat_from_sub(float lonind, float latind, 
			float &lon, float &lat) {
  lat=offset+(180-2*offset)*(latind-1)/nlat-90;
  lon=360*lonind/nlon[(int) latind];
}

inline float lonlat_binsub::lat_from_sub(float latind) { 
  return offset+(180-2*offset)*(latind-1)/nlat-90;
}
