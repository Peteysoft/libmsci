#include "datasets_lang.h"
#include "simple_temp.h"
#include "time_class.h"
#include "string_petey.h"
#include "dependent_temp.h"

#define MAXNARG 100
#define MAXNVAR 1000
#define MAXNHEAP 1000

//datasets:
//this is the "base" "directory":
composite_dataset *base_ds;

//tells us whether we're on the LHS or RHS:
int qualflag;
//when "full" names are encountered, we store the paths here:
//for defining dependent ds's, full names on RHS:
//(scope of datasets outside sub-expressions)
composite_dataset *first_qualifier;
//full names on LHS (scope of datasets within subexpession):
composite_dataset *second_qualifier;

//crude, but works:
//here we store datasets within a sub-expression that need
//to be remembered:
dataset *arglist[MAXNARG];
//here we store any known datasets within a sub-expression that need
//to be remembered:
simple_dataset *argkey[MAXNARG];
int narg;

//list of primitive types:
symbol_table primitive_types;

//the heap--stores temporaries within expressions:
dataset *heap[MAXNHEAP];
int nheap;

//start-up:
void ds_initialize() {
  primitive_types.add("float");
  primitive_types.add("double");
  primitive_types.add("int");
  primitive_types.add("long");
  primitive_types.add("time");
  primitive_types.add("string");
  stack_ptr=0;
  base_ds=new composite_dataset();

}
//call at the end of a statement:
void ds_end_statement();

//allocate datasets of a particular type:
dataset *allocate_simple(char *type) {
  int type;
  dataset *ds;

  type=primitive_types.lookup(type);
  switch (type) {
    case (0): 
      ds=new simple<float>();
      break;
    case (1): 
      ds=new simple<double>();
      break;
    case (2): 
      ds=new simple<int32_t>();
      break;
    case (3): 
      ds=new simple<int64_t>();
      break;
    case (4): 
      ds=new simple<time_class>();
      break;
    case (5): 
      ds=new simple<string_petey>();
      break;
    default:
      ds=NULL;
  }

  return ds;

}

dataset *allocate_dependent(char *type) {
  int type;
  dataset *ds;

  type=primitive_types.lookup(type);
  switch (type) {
    case (0): 
      ds=new dependent<float>(argkey, narg);
      break;
    case (1): 
      ds=new dependent<double>(argkey, narg);
      break;
    case (2): 
      ds=new dependent<int32_t>(argkey, narg);
      break;
    case (3): 
      ds=new dependent<int64_t>(argkey, narg);
      break;
    case (4): 
      ds=new dependent<time_class>(argkey, narg);
      break;
    case (5): 
      ds=new dependent<string_petey>(argkey, narg);
      break;
    default:
      ds=NULL;
  }
  narg=0;
  return ds;

}

//push an intermediate expression onto the heap:
int push_heap(dataset *ds) {
  if (nheap>=MAXNHEAP) {
    yyerror("datasets: Heap overflow\n");
    return -1;
  }
  heap[nheap]=ds;
  nheap++;
  return 0;
}

//clear the heap for re-use:
void clear_heap() {
  for (int i=0; i<nheap; i++) {
    delete heap[i];
  }
  nheap=0;
}

