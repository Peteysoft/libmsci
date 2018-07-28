#include <stdint.h>

namespace libpetey {

template <typename real>
real nd_interpol(real *data, int32_t rank, int64_t *sub, double *coeff) {
  int32_t n=1<<rank;
  real val=0;

  for (int32_t i=0; i<n; i++) {
    val+=data[sub[i]]*coeff[i];
  }
  return val;
}

int interpol_coeff(int32_t rank, int32_t *dim, double *indices,
		int64_t *subscripts, double *coeffs) {

  int32_t lind[rank], hind[rank];
  int32_t l, h;
  int64_t sub, mult;
  long n=1 << rank;
  long order;
  double frac;
  long norm;
  double weight;

//  printf("(%f, %f)\n", indices[0], indices[1]);

  //determine the upper and lower values for the indices:
  norm=1;
  for (int32_t i=0; i<rank; i++) {
    if (dim[i]==1) {
//      l=0;
//      h=0;
      norm*=2;		//also determine the normalization coeff.
    } else {
      l=(int32_t) indices[i];
      if (l>=dim[i]-1) l=dim[i]-2; else if (l<0) l=0;
      lind[i]=l;
      hind[i]=l+1;
    }
  }

  //calculate the weights and apply them:
  for (long i=0; i<n; i++) {
    weight=1;
    sub=0;
    mult=1;
    for (int32_t j=0; j<rank; j++) {
      if (dim[j] == 1) continue;
      order= 1 << j;
      if ((i & order) == 0) {
        frac=(float) hind[j]-indices[j];
        sub=sub+mult*lind[j];
      } else {
        frac=indices[j]-(double) lind[j];
        sub=sub+mult*hind[j];
      }
      weight*=frac;
      mult=mult*dim[j];
   }
   subscripts[i]=sub;
   coeffs[i]=weight/(double) norm;
//   printf("Weight %d= %f\n", i, weight/(double) norm);
  }

//  printf("Normalization coef: %f\n", norm);

  return 0;
}

} //end namespace libpetey

