//Copyright (C) 2007 Peter Mills.  All rights reserved.

#ifndef _LIBAGF__AGF_FCONV_H__INCLUDED
#define _LIBAGF__AGF_FCONV_H__INCLUDED 1

#include <stdio.h>

#include "agf_defs.h"

namespace libagf {

  //scans a single line for features data in LIBSVM format:
  template <class real>
  dim_ta scan_svm_features(char *line, dim_ta *&ind2, real *&raw);

  //gets format code for different data types:
  template <class type>
  void get_format_code(char *code);

  //read classification data containing in an ASCII file (LVQPAK format):
  //returns number of elements or negative error code=line where error occurred
  nel_ta read_lvq(const char *fname, 		//file name
		real_a **&train, 		//returned training data
		cls_ta *&cls, 			//returned class data
		dim_ta &nvar, 			//number of dimensions
		int flags=0);			//1st bit: no header
						//2nd bit: no class data
						//3rd bit omit class data

  //same as above, read from a stream:
  nel_ta read_lvq(FILE *fs, real_a **&train, cls_ta *&cls, dim_ta &nvar, int flags=0);

  //read classification data containing in an ASCII file (LIBSVM format):
  //returns number of elements or negative error code=line where error occurred
  template <class real, class cls_t>
  nel_ta read_svm(const char *fname, 		//file name
		real **&train, 			//returned training data
		cls_t *&cls, 			//returned class data
		dim_ta &nvar, 			//number of dimensions
		real missing,			//value for missing data
		int Uflag);			//re-label classes to go from 0..nc-1

  //same as above, read from a stream:
  template <class real, class cls_t>
  nel_ta read_svm(FILE *fs, real **&train, cls_t *&cls, dim_ta &nvar, real missing, int Uflag);

  //SVM format is actually always floating point for the ordinates:
  template <class real, class cls_t>
  nel_ta read_svm(FILE *fs, real **&train, cls_t *&ord, dim_ta &nvar, real missing=0);

  //print to an ascii file (LVQPAK or LIBSVM format):
  template <class cls_t, class real>
  void print_lvq_svm(FILE *fs, 		//stream
		real **vec,		//features data
		cls_t *cls,		//class data
					//(lvq only: no class data printed if null)
		nel_ta n,		//number of samples
		dim_ta nvar,		//number of features
		int svmflag=0,		//libsvm format
		int nhflag=0);		//do not print header (lvq only)

  //read the output from svm-predict:
  template <class cls_t, class real>
  nel_ta read_svmout(FILE *fs,
		cls_t *&cls,
		real **&p,
		cls_t &ncls,
		nel_ta n=-1);		//if n is greater than 0, cls & p are not 
  					//allocated and exactly n
					//entries are read (if they exist)

}

#endif
