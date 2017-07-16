//
// This software is released under the following terms:
//
// 1. No commercial use.
// 2. Copies and derivative works are free to use and modify.
// 3. Attribution must be given to all contributors of both original and derivative works.
//
// Authors:
//
// 2017-07-16 Peter Mills: added license information 
//

//
// I/O routines for libAGF binary formats.
//

#ifndef _LIBAGF__AGF_IO_H__INCLUDED
#define _LIBAGF__AGF_IO_H__INCLUDED 1

#include <stdio.h>

#include "agf_defs.h"

namespace libagf {

  //-returns null pointer on failure
  //-sample number is -1 if there is bad data

  //read vector and return (coordinate) data from a file:
  template <typename real>
  real ** read_vecfile(const char *filename, 	//datafile name
		dim_ta &m, 			//number of samples
		nel_ta &n);			//number of dimensions

  //reads and returns class data from a file
  template <typename cls_t>
  cls_t * read_clsfile(const char *filename,	//name of data file 
		nel_ta &n);			//number of samples

  //read and return scalar (floating point) data:
  template <typename real>
  real * read_datfile(const char *filename,	//name of data file
		 nel_ta &n);			//number of samples in file


  //read and return normalization/transformation matrix:
  template <typename real>
  real ** read_stats2(const char *filename, 	//name of file
		real *&ave, 			//constant term
		dim_ta &m, 			//first dimension of matrix
		dim_ta &n);			//second      "

  //prints out averages and standard deviations:
  int print_stats(FILE *fs, float *ave, float *std, dim_ta ndim);
  int print_stats(FILE *fs, double *ave, double *std, dim_ta ndim);
  //reads file containing data written by print_stats:
  int read_stats(const char *filename, real_a *ave, real_a *std, dim_ta ndim);
  //read training data (features plus classes):
  template <typename real, typename cls_t>
  int agf_read_train(const char *fbase, real **&train, cls_t *&cls, nel_ta &n, dim_ta &nvar);
  //read borders data (border samples plus gradient vectors)
  template <typename real>
  int agf_read_borders(const char *fbase, real **&brd, real **&grd, nel_ta &n, dim_ta &nvar);

  //gets features training data plus any transformations:
  template <typename real>
  real ** agf_get_features(const char *fbase, agf_command_opts *opts, nel_ta &n, dim_ta &nvar, flag_a sufflag=0);
  //get statistical classification training data:
  int agf_get_train(const char *fbase, agf_command_opts *opts, real_a **&x,
		  cls_ta *&cls, nel_ta &n, dim_ta &nvar);

  //compiles the command for preprocessing the feature data:
  char *compile_precommand(const char *fname, agf_command_opts *opts);

}

#endif
