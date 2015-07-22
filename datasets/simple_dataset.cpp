#include <stdlib.h>
#include <iostream>
#include "simple_dataset.h"
#include "dependent_dataset.h"

namespace libpetey {
namespace datasets {

// not sure why we need this since the class is virtual so an instance
// cannot be created:
simple_dataset::simple_dataset() {
//  printf("Simple base constructor called\n");

  this->type=SIMPLE_BASE;
  ndep=0;
  dependents=NULL;
  rank=NULL;
}

simple_dataset::~simple_dataset() {
//  printf("Simple base destructor called\n");
  if (ndep != 0) {
    delete [] dependents;
    delete [] rank;
  }
}

// adds a dependent
long simple_dataset::add_dependent(dependent_dataset *dep, rank_type r) {
  dependent_dataset **new_dep;
  rank_type *new_rank;
  long i;

  new_dep=new dependent_dataset *[ndep+1];
  new_rank=new rank_type[ndep+1];

  for (i=0; i<ndep;i++) {
    new_dep[i]=dependents[i];
    new_rank[i]=rank[i];
  }

  new_rank[ndep]=r;
  new_dep[ndep]=dep;

  delete [] rank;
  delete [] dependents;

  rank=new_rank;
  dependents=new_dep;

  ndep++;
  
  //cout << "New dependent added: " << ndep << "\n";
  fprintf(stderr, "New dependent added: %ld\n", (int64_t) ndep);

  return ndep;

}

long simple_dataset::del_dependent(dependent_dataset *dep) {
  dependent_dataset **new_dep;
  rank_type *new_rank;
  long i, j;

  //search for the dependent:
  if (ndep==0) return -1;
  for (i=0; i<ndep;i++) {
    if (dep == dependents[i]) break;
  }

  if (i==ndep) return -1;

  if (ndep == 1) {
    delete [] rank;
    delete [] dependents;
    ndep=0L;
    return ndep;
  }

  new_dep=new dependent_dataset*[ndep-1];
  new_rank=new rank_type[ndep-1];

  for (j=0; j<i;j++) {
    new_dep[j]=dependents[j];
    new_rank[j]=rank[j];
  }

  for (j=i+1; j<ndep ; j++) {
    new_dep[j-1]=dependents[j];
    new_rank[j-1]=rank[j];
  }

  delete [] rank;
  delete [] dependents;

  rank=new_rank;
  dependents=new_dep;

  ndep--;

  return ndep;

}

} //end namespace datasets
} //end namespace libpetey
