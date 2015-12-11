//Copyright (C) 2007 Peter Mills.  All rights reserved.

#ifndef _LIBAGF__AGF_IO_H__INCLUDED
#define _LIBAGF__AGF_IO_H__INCLUDED 1

#include <stdio.h>

#include "agf_defs.h"

namespace libagf {

  //-returns null pointer on failure
  //-sample number is -1 if there is bad data

  //read vector and return (coordinate) data from a file:
  real_a ** read_vecfile(const char *filename, 	//datafile name
		dim_ta &m, 			//number of samples
		nel_ta &n);			//number of dimensions

  //reads and returns class data from a file
  cls_ta * read_clsfile(const char *filename,	//name of data file 
		nel_ta &n);			//number of samples

  //read and return scalar (floating point) data:
  real_a * read_datfile(const char *filename,	//name of data file
		 nel_ta &n);			//number of samples in file


  //read and return normalization/transformation matrix:
  real_a ** read_stats2(const char *filename, 	//name of file
		real_a *&ave, 			//constant term
		dim_ta &m, 			//first dimension of matrix
		dim_ta &n);			//second      "

  //prints out averages and standard deviations:
  int print_stats(FILE *fs, real_a *ave, real_a *std, dim_ta ndim);
  //reads file containing data written by print_stats:
  int read_stats(const char *filename, real_a *ave, real_a *std, dim_ta ndim);
  //read training data (features plus classes):
  int agf_read_train(const char *fbase, real_a **&train, cls_ta *&cls, nel_ta &n, dim_ta &nvar);
  //read borders data (border samples plus gradient vectors)
  int agf_read_borders(const char *fbase, real_a **&brd, real_a **&grd, nel_ta &n, dim_ta &nvar);

  //gets features training data plus any transformations:
  real_a ** agf_get_features(const char *fbase, agf_command_opts *opts, nel_ta &n, dim_ta &nvar, flag_a sufflag=0);
  //get statistical classification training data:
  int agf_get_train(const char *fbase, agf_command_opts *opts, real_a **&x,
		  cls_ta *&cls, nel_ta &n, dim_ta &nvar);

  //compiles the command for preprocessing the feature data:
  char *compile_precommand(const char *fname, agf_command_opts *opts);

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
