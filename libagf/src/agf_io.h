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

}

#endif
