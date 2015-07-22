#ifndef STRINGS_INCLUDED
#define STRINGS_INCLUDED

//defines a string type and a set of operators
#include <stdio.h>
#include <iostream>

namespace libpetey {

  class string_petey {
    //friend std::ostream &operator << (ostream &outfile, string_petey &s);
    protected:
      char *stuff;
      long size;
    public:
      //constructors and destructors:
      string_petey();
      string_petey(const char *astring);		//construction from char arr
      string_petey(const string_petey &astring);		//construction from another string
      string_petey &operator = (const char *astring);	//assignment from char array
      string_petey &operator = (const string_petey &astring);	//assignment
      ~string_petey();

      //"arithmetic":
      string_petey &operator + (const string_petey &string2);	//concatination

      //comparators:
      int cmp (const string_petey &string2);		//general...
      int operator > (const string_petey &string2);
      int operator < (const string_petey &string2);
      int operator >= (const string_petey &string2);
      int operator <= (const string_petey &string2);
      int operator == (const string_petey &string2);
      int operator != (const string_petey &string2);

      //substrings:
      char el(long index);
      char operator [] (long index);		//get an element of the string
      string_petey mid(long start, long length);	//get a substring
      //separate into substrings base on a delimiter:
      string_petey * sep(string_petey &del);
      string_petey * sep(char del);

      //type conversion:
      operator char* ();
      operator double ();

      //miscellaneous:
      void print();
      long length();

      //io:
      long read(FILE *fileptr);
      long write(FILE *fileptr);
  };

  //std::ostream &operator << (std::ostream &outfile, string_petey &s);
  //std::istream &operator >> (std::istream &outfile, string_petey &s);

  #ifndef __STRING__
    typedef string_petey string;
  #endif

} //end namespace libpetey

#endif
