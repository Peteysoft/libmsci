template <class type>
void naive_sort(type *array, int n) {
  type swp;

  for (int i=0; i<n; i++) {
    for (int j=0; j<i; j++) {
      if (array[i]<array[j]) {
        swp=array[j];
        array[j]=array[i];
        array[i]=swp;
        //we're deliberately writing an inefficient sorting routine,
        //why not throw in a goto while we're at it?
        goto skip;
      }
    }
    for (int j=i+1; j<n; j++) {
      if (array[j]<array[i]) {
        swp=array[j];
        array[j]=array[i];
        array[i]=swp;
        //we're deliberately writing an inefficient sorting routine,
        //why not throw in a goto while we're at it?
        goto skip;
      }
    }
    skip:
  }
}
        
