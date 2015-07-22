
template <class dtype>
long dependent<dtype>::interpol(dtype &value, double *indices) {
  long lind[rank], hind[rank];
  long l, h;
  long sub, mult;
  long n=1 << rank;
  long order;
  double dist, norm, frac;
  double weights[n];

  if (n_data == 0) return -1;

  //determine the upper and lower values for the indices:
  for (long i=0; i<rank; i++) {
    l=(long) indices[i];
    if (dim[i]==1) {
      l=0;
      h=0;
    } else {
      if (l>=dim[i]-1) l=dim[i]-2; else if (l<0) l=0;
      h=l+1;
    }
    lind[i]=l;
    hind[i]=h;
  }

  //determine the weighting coefficients for the array elements at each vertex:
  for (long i=0; i<n; i++) weights[i]=1;

  for (long i=0; i<n; i++) {
    sub=0;
    mult=1;
    dist=0;
    for (long j=0; j<rank; j++) {
      if (dim[j] == 1) {
        dist=dist+0.01;
	continue;
      }
      order= 1 << j;
      if ((i & order) == 0) {
        frac=indices[j]-(float) lind[j];
        dist=dist+frac*frac;
      } else {
        frac=indices[j]-(float) hind[j];
        dist=dist+frac*frac;
      }
      mult=mult*dim[j];
    }
    for (long j=0; j<i; j++) weights[j]*=dist;
    for (long j=i+1; j<n; j++) weights[j]*=dist;
  }

  //apply the weights and calculate the normalization coefficient:
  norm=0;
  value=0;
  for (long i=0; i<n; i++) {
    sub=0;
    mult=1;
    for (long j=0; j<rank; j++) {
      order= 1 << j;
      if ((i & order) == 0) {
        sub=sub+mult*lind[j];
      } else {
        sub=sub+mult*hind[j];
      }
      mult=mult*dim[j];
   }
    value+=(dtype) weights[i]*data[sub];
    norm+=weights[i];
  }
  value=value/(dtype) norm;

  return 0;
}
