#include <iostream.h>
#include "simple2.cpp"

//test the "simple_dataset" class:
int main() {

  simple_dataset<float> list;
  long i, index, n;
  long * mapping;

  list.add(4.0);
  list.add(2.0);
  list.add(3.0);
  list.add(1.0);
  list.add(6.0);

  for (i=0;i<list.size_of();i++) cout << list[i] << " ";
  cout << "\n";
  
  //  index=list.get_ind(5);

//  cout << index << "\n";

  list.add(5.0);
  index=list.search(5.0);
  cout << index << "\n";
  for (i=0;i<list.size_of();i++) cout << list[i] << " ";
  return 1;
}
