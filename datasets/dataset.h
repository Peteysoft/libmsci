//
// Defines the base virtual class of all the dataset classes.
//
// Inheritance structure of the dataset classes goes something like this:
//
// dataset
//    |
//    | --  simple_dataset
//    |           |
//    |           | --  simple<float>
//    |           | --  simple<long>
//    |           | --  simple<time>
//    |           | --  simple<... >
//    |
//    | --  dependent_dataset
//    |           |
//    |           | --  dependent<float>
//    |           | --  dependent<long>
//    |           | --  dependent<time>
//    |           | --  dependent<... >
//    |
//    | --  composite_dataset
//
// Compilation dependencies go something like this:
//
//                          /---->  simple<...>
//                          |
//             /---> simple_dataset -----\
//             |         /| |            |
//             |          | |/           |
// dataset ----|---> dependent_dataset -------->  composite_dataset
//             |           |             |
//             |      dependent<...>     |
//             |                         |
//             \-------------------------/
//
// Planned extensions include simple dataset classes for holding only linearly
// gridded data and dependent dataset classes that swap from a file, but where these
// will fit in the hierarchy is as yet uncertain.  There is also a dependent class
// which holds no data and includes only one method for returning the subscripts and
// coefficients of an interpolated point, but this is meant as a one-off for the
// the tracer simulation.
//
// Other extensions might include making insertions more efficient for both the
// simple and dependent classes, although at the expense of storage and possibly
// retrieval time.

#ifndef DATASET_INCLUDED
#define DATASET_INCLUDED

#include <stdio.h>
#include <stdint.h>

#define UNDEF 0

#define SIMPLE_BASE 100

#define SIMPLE_FLOAT 101
#define SIMPLE_DOUBLE 102
#define SIMPLE_INT32 103
#define SIMPLE_INT64 104
#define SIMPLE_TIME 105
#define SIMPLE_STRING 106
#define SIMPLE_LONG 110

#define DEPENDENT 200

#define DEPENDENT_FLOAT 201
#define DEPENDENT_DOUBLE 202
#define DEPENDENT_INT32 203
#define DEPENDENT_INT64 206
#define DEPENDENT_TIME 204
#define DEPENDENT_STRING 205

#define DEPENDENT_FLOAT_S 211
#define DEPENDENT_DOUBLE_S 212
#define DEPENDENT_LONG_S 213

#define COMPOSITE 300

#define FLOAT_MISSING 0.0
#define DOUBLE_MISSING 0.0
#define INT32_MISSING 0
#define INT64_MISSING 0

//error codes:
#define NO_PROBLEM 0
#define DATA_ELEMENT_NOT_FOUND -1
#define RANK_ERROR -2
#define SUBSCRIPT_OUT_OF_RANGE -3
#define MISC_INDEXING_ERROR -4
#define NO_DATA -5
#define INDICES_HAVE_WRAPPED -200
#define INDEX1_OUT_OF_RANGE -101
#define INDEX2_OUT_OF_RANGE -102
#define INDEX3_OUT_OF_RANGE -103
#define INDEX4_OUT_OF_RANGE -104
#define INDEX5_OUT_OF_RANGE -105
#define INDEX6_OUT_OF_RANGE -106
#define INDEX7_OUT_OF_RANGE -107
#define INDEX8_OUT_OF_RANGE -108
#define INDEX9_OUT_OF_RANGE -109
#define INDEX10_OUT_OF_RANGE -110
#define INDEX11_OUT_OF_RANGE -111
#define INDEX12_OUT_OF_RANGE -112
#define INDEX13_OUT_OF_RANGE -113
#define INDEX14_OUT_OF_RANGE -114
#define INDEX15_OUT_OF_RANGE -115
#define INDEX16_OUT_OF_RANGE -116
#define INDEX17_OUT_OF_RANGE -117
#define INDEX18_OUT_OF_RANGE -118
#define INDEX19_OUT_OF_RANGE -119
#define INDEX20_OUT_OF_RANGE -120

/*
typedef int16_t dstype_type;	//for holding dataset types
typedef int32_t sub_1d_type;	//for holding 1-dimensional subscripts
typedef int32_t ind_type;		//for holding multi-dimensional indices
typedef int16_t rank_type;		//for holding the rank of the dataset
*/

namespace libpetey {
namespace datasets {

typedef int dstype_type;	//for holding dataset types
typedef int32_t sub_1d_type;	//for holding 1-dimensional subscripts
typedef int32_t ind_type;		//for holding multi-dimensional indices
typedef int16_t rank_type;		//for holding the rank of the dataset

typedef sub_1d_type errtype;	//for returning errors
typedef double interpol_index;	//for doing interpolation
typedef long counter_type;	//for doing anything else!

extern errtype INDEX_OUT_OF_RANGE[20];

extern FILE *datasets_log;

class dataset {
protected:
  dstype_type type;	//the type of the dataset.  Types are defined above.
  sub_1d_type n_data;	//the number of data elements.

public:
  dataset();
  virtual ~dataset();
  virtual long read(FILE *fileptr)=0;		//read from a file
  virtual long write(FILE *fileptr)=0;		//write to a file
/*
  virtual dataset *add (const dataset &other)=0;
  virtual dataset *subtract (const dataset &other)=0;
  virtual dataset *multiply (const dataset &other)=0;
  virtual dataset *divide (const dataset &other)=0;
*/
//	dataset()=0;
//	~dataset()=0;

  sub_1d_type nel();		//returns the number of elements
  dstype_type type_of();	//returns the type
};

} //end namespace datasets
} //end namespace libpetey

#endif

