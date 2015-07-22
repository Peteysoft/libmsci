#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv) {

  float c0;

  if (argc != 2) {
    printf("Usage: c0 lat\n");
    printf("--creates contour at latitude lat\n");
    exit(-1);
  }

  sscanf(argv[1], "%f", &c0);

  for (float lon=0.; lon < 360.; lon++) {
    printf("%f %f\n", lon, c0);
  }

}

