#include "../simple_temp.h"
#include <stdio.h>
#include <iostream.h>

#pragma interface

int main() {
  simple<float> data;
  float value;
  char let;
  char line[100];
  long size, i;

  do {
//    let=getc(stdin);
    cin >> let;
    switch (let) {
      case 'n':
		printf("%d\n", data.nel());
		break;
      case 'a':
//      		fgets(line, size, stdin);
//		sscanf("%f", line, &value);
		cin >> value;
		printf("%d\n", data.add_el(value));
		break;
      case 'p':
      		for (i=0; i<data.nel(); i++) {
		  data.get(value, i);
		  printf("%f ", value);
		}
		printf("\n");
		break;
      case 'd':
      		//fgets(line, size, stdin);
		//sscanf("%f", line, &value);
		cin >> value;
		printf("%d\n", data.del(value));
		break;
      case 's':
      		//fgets(line, size, stdin);
		//sscanf("%f", line, &value);
		cin >> value;
		printf("%d\n", data.search(value));
		break;
      case 'g':
      		//fgets(line, size, stdin);
		//sscanf("%d", line, &i);
		cin >> i;
		data.get(value, i);
		printf("%f\n", value);
		break;
      case 'i':
      		cin >> value;
		printf("%f\n", data.interp(value));
		break;
    }
  } while (let != 'q');

}

