#include <stdio.h>
#include <string.h>
#include <Magick++.h>

using namespace std;
using namespace Magick;

#define NCHAN 2
#define GRAIN 1000

int main(int argc, char **argv) {
  FILE *fs;
  char *fname;
  Image chan;
  PixelPacket *row;
  Geometry size;
  int blen;
  unsigned char **data;
  int n;		//total number of data elements
  int count;		//cumulative count of data elements read in

  blen=strlen(argv[1]);

  for (int i=0; i<NCHAN; i++) {
    fname=new char[blen+8];
    sprintf(fname, "%s_B%1.1d.TIF", argv[1], i+1);

    try {
      chan.read(fname);
    }
    catch(Magick::WarningCoder &warning) {
      fprintf(stderr, "ImageMagick threw the following warning: %s\n", warning.what());
      if (i!=0) fprintf(stderr, "%5.1f%%", 100.*count/n);
    }

    //do this the easy way:
    if (i==0) {
      size=chan.size();
      data=new unsigned char * [size.height()*size.width()];
      n=size.height()*size.width()*NCHAN;
      count=0;
      data[0]=new unsigned char[n];
      for (int j=1; j<size.height()*size.width(); j++) data[j]=data[0]+j*NCHAN;
      fprintf(stderr, "%5.1f%%", 0.);
    }

    //do this the easy way:
    for (int j=0; j<size.height(); j++) {
      row=chan.getPixels(0, j, size.width(), 1);
      for (int k=0; k<size.width(); k++) {
        data[size.width()*j+k][i]=row[k].blue;
	//float val=data[size.width()*j+k][i];
	//printf("%hhd %f\n", row[k].blue, val);
	count++;
	if (count % GRAIN == 0) {
          fprintf(stderr, "\b\b\b\b\b\b%5.1f%%", 100.*count/n);
	}
      }
    }
  }

  if (argc < 3) {
    for (int i=0; i<size.height()*size.width(); i++) {
      for (int j=0; j<NCHAN; j++) {
        printf("%5hhd", data[i][j]);
      }
      printf("\n");
    }
  } else {
    float rec[NCHAN];
    int32_t nchan=NCHAN;
    fs=fopen(argv[2], "w");
    fwrite(&nchan, sizeof(nchan), 1, fs);
    printf("%dx%d\n", size.width(), size.height());
    for (int i=0; i<size.height()*size.width(); i++) {
      for (int j=0; j<nchan; j++) rec[j]=data[i][j];
      printf("%hhd %f\n", data[i][0], rec[0]);
      fwrite(rec, sizeof(float), nchan, fs);
    }
    fclose(fs);
  }

  delete fname;
  delete data[0];
  delete data;

}
