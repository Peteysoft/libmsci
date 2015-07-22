#include <iostream.h>
#include "simple_temp.h"

//test the "simple_dataset" class:
int main() {

  simple<float> list;
  long i, index, n;
  float value;

  list.add_el(4.0);
//  cout << list[0] << "\n";
  cout << list.nel() << "\n";
  cout << list.search(4.0) << "\n";
  cout << list.search(3.0) << "\n";
  cout << list.search(5.0) << "\n";
  list.add_el(2.0);
  list.add_el(3.0);
  list.add_el(1.0);
  list.add_el(6.0);

  for (i=0;i<list.nel();i++) {
	  list.get(value, i);
	  cout << value << " ";
  }
  cout << "\n";
  cout << "\n";

  //  index=list.get_ind(5);

//  cout << index << "\n";

  list.add_el(5.0);
  index=list.search(5.0);
  cout << index << "\n";
  for (i=0;i<list.nel();i++) {
	  list.get(value, i);
	  cout << value << " ";
  }

  list.del(4.0);
  for (i=0;i<list.nel();i++) {
	  list.get(value, i);
	  cout << value << " ";
  }
  list.del(2.0);
  list.del(0.0);
  list.del(6.0);
  list.del(3.0);
  list.del(4.0);
  list.del(5.0);

  cout << list.nel() << "\n";

  list.del(1.0);
  cout << list.nel() << "\n";

  return 1;
}

