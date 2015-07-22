#include <math.h>

#include "ctraj_defaults.h"

namespace ctraj {

//approximate squared distance for lon-lat coordinates:
float sdist(float lon1, float lat1, float lon2, float lat2) {
  float latm;
  float diff1, diff2;
  int div;

  latm=(lat1+lat2)/2;

  if (fabs(latm) > LAT_THRESH) {
    float r1, r2;
    float x1, x2, y1, y2;
    //use a locally cartesian system close to the poles:
    r1=REARTH*DEG2RAD*(90-fabs(lat1));
    r2=REARTH*DEG2RAD*(90-fabs(lat2));
    x1=r1*cos(DEG2RAD*lon1);
    y1=r1*sin(DEG2RAD*lon1);
    x2=r2*cos(DEG2RAD*lon2);
    y2=r2*sin(DEG2RAD*lon2);
    diff1=x1-x2;
    diff2=y1-y2;
  } else {
    float dlon;

    dlon=fabs(lon1-lon2);
    div=dlon/360;
    dlon=dlon - 360*div;
    if (dlon > 180) dlon=360-dlon;

    diff1=dlon*KMPERDEG*cos(latm*DEG2RAD);
    diff2=(lat1-lat2)*KMPERDEG;
  }

  return diff1*diff1+diff2*diff2;

}

//interpolate between two longitude and latitude coordinates:
void interpolate_lon_lat(float lon1, float lat1, float lon2, float lat2,
		double frac, float &lon, float &lat) {
	
  if (fabs(lat1) > LAT_THRESH || fabs(lat2) > LAT_THRESH) {
    float r1, r2;
    float x1, x2, y1, y2;
    float x, y;
    //use a locally cartesian system close to the poles:
    r1=REARTH*DEG2RAD*(90-fabs(lat1));
    r2=REARTH*DEG2RAD*(90-fabs(lat2));
    x1=r1*cos(DEG2RAD*lon1);
    y1=r1*sin(DEG2RAD*lon1);
    x2=r2*cos(DEG2RAD*lon2);
    y2=r2*sin(DEG2RAD*lon2);

    //do the interpolation:
    x=(1-frac)*x1+frac*x2;
    y=(1-frac)*y1+frac*y2;

    //convert back:
    lon=atan2f(y, x)/DEG2RAD;
    lat=90-sqrt(x*x+y*y)/REARTH/DEG2RAD;
    if (lat1 < 0) lat=-lat;
  } else {
    float dlon;

    dlon=lon2-lon1;
    if (fabs(dlon) > 180) dlon=dlon-fabs(dlon)*360/dlon;

    //interplation is simple linear away from poles:
    lon=lon1+frac*dlon;
    lat=(1-frac)*lat1+frac*lat2;
  }
}

} //end namespace ctraj

