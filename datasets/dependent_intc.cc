#include "dependent_intc.h"

namespace libpetey {
namespace datasets {

dependent_intc::dependent_intc(stype **s, rank_type ndep) {
  rank_type i;

  rank=ndep;
  dim=new ind_type[rank];
  dependencies=new stype *[rank];

  n_data=1;
  for (i=0; i<rank; i++) {
    dependencies[i]=s[i];
    s[i]->add_dependent(this, i);
    dim[i]=dependencies[i]->nel();
    n_data*=dim[i];
  }

}

dependent_intc::dependent_intc(stype *s1) {

	rank=1;
	dim=new ind_type[rank];
	dependencies=new stype * [rank];

	n_data=s1->nel();
	dim[0]=n_data;
	if (n_data==0) return;
	dependencies[0]=s1;

}

dependent_intc::dependent_intc(stype *s1, stype *s2) {

	rank=2;
	dim=new ind_type[rank];
	dependencies=new stype * [rank];

	dim[0]=s1->nel();
	dim[1]=s2->nel();
	dependencies[0]=s1;
	dependencies[1]=s2;

	n_data=dim[0]*dim[1];
}
dependent_intc::dependent_intc(stype *s1, stype *s2, stype *s3) {

	rank=3;
	dim=new ind_type[rank];
	dependencies=new stype * [rank];

	dim[0]=s1->nel();
	dim[1]=s2->nel();
	dim[2]=s2->nel();
	dependencies[0]=s1;
	dependencies[1]=s2;
	dependencies[2]=s3;

	n_data=dim[0]*dim[1]*dim[2];
}

dependent_intc::dependent_intc(stype *s1, stype *s2, stype *s3, stype *s4) {

	rank=4;
	dim=new ind_type[rank];
	dependencies=new stype * [rank];

	dim[0]=s1->nel();
	dim[1]=s2->nel();
	dim[2]=s3->nel();
	dim[3]=s4->nel();
	dependencies[0]=s1;
	dependencies[1]=s2;
	dependencies[2]=s3;
	dependencies[3]=s4;

	n_data=dim[0]*dim[1]*dim[2]*dim[3];
}

errtype dependent_intc::interpol_coeff(interpol_index *indices,
		sub_1d_type *subscripts, double *coeffs) {

  ind_type lind[rank], hind[rank];
  ind_type l, h;
  sub_1d_type sub, mult;
  long n=1 << rank;
  long order;
  interpol_index frac;
  long norm;
  double weight;

  if (n_data == 0) return NO_DATA;

//  printf("(%f, %f)\n", indices[0], indices[1]);

  //determine the upper and lower values for the indices:
  norm=1;
  for (rank_type i=0; i<rank; i++) {
    if (dim[i]==1) {
//      l=0;
//      h=0;
      norm*=2;		//also determine the normalization coeff.
    } else {
      l=(ind_type) indices[i];
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
    for (rank_type j=0; j<rank; j++) {
      if (dim[j] == 1) continue;
      order= 1 << j;
      if ((i & order) == 0) {
        frac=(interpol_index) hind[j]-indices[j];
        sub=sub+mult*lind[j];
      } else {
        frac=indices[j]-(interpol_index) lind[j];
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

  return NO_PROBLEM;
}

void dependent_intc::print_meta() {
  printf("Type: %d\n", type);
  for (rank_type i=0; i<rank-1; i++) printf("%d X ", dim[i]);
  printf("%d = %d\n", dim[rank-1], n_data);
  for (rank_type i=0; i<rank; i++) {
    dependencies[i]->print();
  }
}

} //end namespace datasets
} //end namespace libpetey

