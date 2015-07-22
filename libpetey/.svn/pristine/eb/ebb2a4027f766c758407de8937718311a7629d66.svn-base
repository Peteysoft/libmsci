#ifndef READ_ASCII_ALL_H_INCLUDED
#define READ_ASCII_ALL_H_INCLUDED

#include <stdio.h>

namespace libpetey {

  //reads an entire file in and returns it in an array of strings:
  //uses a linked list...
  char ** read_ascii_all(FILE *fs, long *n, int flag=0);
  char ** read_ascii_all(const char *filename, long *n, int flag=0);

  //better version of fgets:
  char * fget_line(FILE *fs, int flag=0);

  //splits a string based on whitespace
  //returns a set of integer indexing into the beginning of each "word"
  //(seems kind of useless; what about the lengths and endings of the
  //substrings?)
  int * split_string(const char *string, 	//string
		  int &n);			//number of substrings

  //lets make a better one:
  char ** split_string_destructive(char *string, 	//string
		  int &n);			//number of substrings

}

#endif

