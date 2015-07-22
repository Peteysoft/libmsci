#include "dataset.h"
#include "simple_dataset.h"
#include "symbol_table.h"

#define STACK_HEIGHT 100
#define MAXNVAR 1000

//datasets:
//this is the "base" "directory":
extern composite_dataset *base_ds;

//tells us whether we're on the LHS or RHS:
extern int qualflag;
//when "full" names are encountered, we store the paths here:
//for defining dependent ds's, full names on RHS:
//(scope of datasets outside sub-expressions)
extern composite_dataset *first_qualifier;
//full names on LHS (scope of datasets within subexpession):
extern composite_dataset *second_qualifier;

//crude, but works:
//here we store datasets within a sub-expression that need
//to be remembered:
extern dataset *arglist[STACK_HEIGHT];
//here we store any known datasets within a sub-expression that need
//to be remembered:
extern simple_dataset *argkey[STACK_HEIGHT];
extern int narg;

//list of primitive types:
extern symbol_table primitive_types;

//start-up:
void ds_initialize();
//call at the end of a statement:
void ds_end_statement();

//allocate datasets of a particular type:
dataset *allocate_simple(char *type);
dataset *allocate_dependent(char *type);


