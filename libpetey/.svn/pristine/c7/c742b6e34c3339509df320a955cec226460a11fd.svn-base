//defines a string type and a set of operators
#include <stdlib.h>
#include <iostream>
#include "string_petey.h"
#include "linked.h"

namespace libpetey {

string_petey::string_petey() {
  stuff=new char[1];
  stuff[0]='\x0';
  size=0;
}

string_petey::string_petey(const char *astring) {
  long n,i;
  for (n=0;astring[n] != '\x0';n++) {
//    cout << n;
  }
  stuff=new char[n+1];
  size=n;
  for (i=0;i<=n;i++) stuff[i]=astring[i];
  //cout << "Copy constructor (char *) copied: " << stuff << "; size: " << size << "\n";
}

//copy constructor:
string_petey::string_petey(const string_petey &astring) {
  long i;
//  if (stuff != NULL) delete stuff;
  size=astring.size;
  stuff=new char[size+1];
  for (i=0;i<=size;i++) stuff[i]=astring.stuff[i];
  //cout << "Copy constructor copied: " << stuff << "; size: " << size << "\n";
}

string_petey &string_petey::operator = (const string_petey &astring) {
  long n,i;
  delete [] stuff;
  size=astring.size;
  stuff=new char[size+1];
  for (i=0;i<=size;i++) stuff[i]=astring.stuff[i];
  //cout << "= operator copied: " << astring.stuff << "; size: " << size << "\n";
//  if (stuff[size]!='\0') cout << "Null char. not copied\n";
  return *this;
}

string_petey::~string_petey() {
  //cout << "Deleting: " << stuff << "\n";
  delete [] stuff;
}

string_petey &string_petey::operator = (const char *astring) {
  long n,i;
  delete [] stuff;
  for (n=0;astring[n] != '\x0';n++) {
//    cout << n;
  }
  size=n;
  stuff=new char[n+1];
  for (i=0;i<=n;i++) stuff[i]=astring[i];
  //cout << "= (char *) operator copied: " << stuff << "; size: " << size << "\n";
  //printf("= (char *) operator copied: %s; size: %d\n", stuff, size);
  return *this;
}

string_petey &string_petey::operator + (const string_petey &string2) {
  string_petey *newstring;
  long i;
  newstring=new string_petey;

  newstring->size=size+string2.size;
//  cout << newstring->size << "\n";
  newstring->stuff=new char[newstring->size];

  for (i=0;i<size;i++) {
    newstring->stuff[i]=stuff[i];
  }
  for (i=0;i<string2.size;i++) {
    newstring->stuff[i+size]=string2.stuff[i];
  }
  newstring->stuff[newstring->size]='\x0';
  //cout << newstring->stuff << "\n";
  //printf("%s\n", newstring->stuff);
  return *newstring;
}

/*int strcmp(string string1, string string2) {
  long i, min_size;
  int exit;

  if (string1.size > string2.size) {
    exit=1;
    min_size=string2.size;
  } else if (string1.size == string2.size) {
    exit=0;
    min_size=string2.size;
  } else {
    exit=-1;
    min_size=string1.size;
  }
  for (i=0;i<min_size;i++) {
    if (string1.stuff[i] > string2.stuff[i]) {
      exit=1;
      break;
    } else if (string1.stuff[i] < string2.stuff[i]) {
      exit=-1;
      break;
    }
  }
  return exit;
}*/

int strcmp(string_petey string1, string_petey string2) {
  return string1.cmp(string2);
}

int string_petey::cmp(const string_petey &string2) {
  long i, min_size;
  int exit;

  if (size > string2.size) {
    exit=1;
    min_size=string2.size;
  } else if (size == string2.size) {
    exit=0;
    min_size=string2.size;
  } else {
    exit=-1;
    min_size=size;
  }
  for (i=0;i<min_size;i++) {
    if (stuff[i] > string2.stuff[i]) {
      exit=1;
      break;
    } else if (stuff[i] < string2.stuff[i]) {
      exit=-1;
      break;
    }
  }
  return exit;
}

int string_petey::operator > (const string_petey &string2) {
  if (cmp(string2)>0) return 1;
  return 0;
}

int string_petey::operator < (const string_petey &string2) {
  if (cmp(string2)<0) return 1;
  return 0;
}

int string_petey::operator >= (const string_petey &string2) {
  if (cmp(string2)>=0) return 1;
  return 0;
}

int string_petey::operator <= (const string_petey &string2) {
  if (cmp(string2)<=0) return 1;
  return 0;
}

int string_petey::operator == (const string_petey &string2) {
  if (cmp(string2)==0) return 1;
  return 0;
}

int string_petey::operator != (const string_petey &string2) {
  if (cmp(string2)!=0) return 1;
  return 0;
}

char string_petey::el(long index) {
  if (index > size || index < 0) return '\0';
  return stuff[index];
}

char string_petey::operator [] (long index) {
  if (index > size || index < 0) return '\0';
  return stuff[index];
}


//get a substring:
string_petey string_petey::mid(long start, long length) {
  int index;
  string_petey newstring;
//  newstring=new string_petey;
  if (start >= size) return newstring;
  if (start+length >= size) length=size-start;
  newstring.size=length+1;
  delete [] newstring.stuff;
  newstring.stuff=new char[length+1];
  for (index=0; index < length; index++) {
    newstring.stuff[index]=stuff[index+start];
  }
  newstring.stuff[length]='\0';
  return newstring;
}

/*
//separate into substrings based on a delimiter:
string_petey * string_petey::sep(string_petey &del, long &n) {
  long dsize;
  long i, j;
  linked_list<string_petey> sarr;

  dsize=del.length();
  i=0; j=0;
  for (;;) {
    for (j=i; j<size; j++) {
      if (mid(j, dsize)==del) break;
    }
    sarr.add(mid(i, j-i));
    if (j>=size) break;
    i=j+dsize;
  }

  return sarr.make_array(n);

}

string_petey * string_petey::sep(char del, long &n) {
  long dsize;
  long i, j;
  linked_list<string_petey> sarr;

  dsize=del.length();
  i=0; j=0;
  for (;;) {
    for (j=i; j<size; j++) {
      if (el(j)==del) break;
    }
    sarr.add(mid(i, j-i));
    if (j>=size) break;
    for (i=j; i<size; i++) {
      if (el(i)!=del) break;
    }
    if (i>=size) break;
  }

  return sarr.make_array(n);

}
*/

string_petey::operator char* () {
  char * result;
  
  result=new char[size+1];
  
  for (long i=0;i<=size;i++) {
    result[i]=stuff[i];
  }

  return result;

}

string_petey::operator double() {
  double result;

  sscanf(stuff, "%lg", &result);

  return result;

}

void string_petey::print() {
  //std::cout << stuff << "\n";
  printf("%s\n", stuff);
}

long string_petey::length() {
  return size;
}

long string_petey::read(FILE *fileptr) {
  long nread;

  if (stuff != NULL) delete [] stuff;

  nread=fread(&size, 1, sizeof(size), fileptr);
  stuff=new char[size+1];
  nread+=fread(stuff, 1, sizeof(char)*size, fileptr);

  stuff[size]='\0';

  return nread;
}

long string_petey::write(FILE *fileptr) {
  long nwritten;

  nwritten=fwrite(&size, 1, sizeof(size), fileptr);
  nwritten+=fwrite(stuff, 1, sizeof(char)*size, fileptr);

  return nwritten;
}

/*
std::ostream &operator << (ostream &outfile, string_petey &s) {
  outfile << s.stuff;
  return outfile;
}
*/

} //end namespace libpetey

