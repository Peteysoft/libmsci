//
// Copywrite 2004 Peter Mills.  All rights reserved.
//
// A class definition for working with dates and times.
//

#ifndef TIME_INCLUDED
#define TIME_INCLUDED

//#include <iostream>

//port of "time_lib.pro" to C++

#define SEPARATORS "//-::"

// Returns the number of leap years between y1 and y2:

namespace libpetey {

const int MAX_TSTR_WIDTH=23;

inline int nleap(int y1, int y2) {
  int y1m1, y2m1;
  y1m1=y1-1;
  y2m1=y2-1;
  return y2m1/4 - y1m1/4 + y1m1/100 - y2m1/100 + y2m1/400 - y1m1/400;
}

class time_class {

  protected:
    int yr;
    double dy;
  public:

    //constructors:
    time_class();
    time_class(int y, int m, int d, int h, int min, float s);
    time_class(int y, double d);
    time_class(const char *t);
    time_class(const double t);
    time_class(const time_class &t);		//copy constructor

    //todays date:
    void now();

    //assignment operators:
    time_class & operator = (const time_class &t);
    time_class & operator = (const char *t);
    time_class & operator = (const double t);

    //initialization functions:
    int init(int y, int mon, int d, int h, int min, float s);
    int read_string(const char *t, const char *sep=SEPARATORS);
    int read_string(const char *t, long &i, const char *sep=SEPARATORS);

    double diff(const time_class &t2) const;	//difference between two time_classs in days
    time_class & add(double d);		//adds d days to time

    //comparators:
    int operator == (const time_class &t2);
    int operator != (const time_class &t2);
    int operator < (const time_class &t2);
    int operator > (const time_class &t2);
    int operator >= (const time_class &t2);
    int operator <= (const time_class &t2);

    //arithmetic operators:
    time_class operator + (time_class &t2);
    time_class operator - (const time_class &t2);
    time_class operator * (double n);
    time_class operator / (double n);
    double operator / (time_class &t2);

    //type conversion:
    char *write_string(char *time_string, const char *sep=SEPARATORS);
    operator double();
    operator char *();

    //extract fields:
    int year();
    int month();
    int day();
    double doy();	//day if year
    int dow();		//day of week
    int hour();
    int minute();
    float second();

    template <class integer, class real>
    void get_fields(integer &y, integer &mon, integer &d,
    		integer &h, integer &min, real &s);

    //modify fields:
    time_class & year(int y);
    time_class & month(int m);
    time_class & day(int d);
    time_class & hour(int h);
    time_class & minute(int m);
    time_class & second(float s);

};

extern time_class ref_date;

// Initialize as 0:
inline time_class::time_class() {
  yr=0;
  dy=0;
}

// Initialize from a string:
inline time_class::time_class(const char *t) {
  read_string(t);
}

// Initialize the fields directly:
inline time_class::time_class(int y, double d) {
  yr=y;
  dy=d;
}

// Specify all the date fields in turn:
inline time_class::time_class(	int y,		//year
			int m,		//month
			int d,		//day of the month
			int h,		//hour
			int min,	//minutes
			float s) {	//seconds
  init(y, m, d, h, min, s);
}

//initialize from a floating point value giving the number of days elapsed
//from the global variable ref_date:
inline time_class::time_class(const double t) {
  yr=ref_date.yr;
  dy=ref_date.dy;

  add(t);
}

//copy constructor:
inline time_class::time_class(const time_class &t) {
  yr=t.yr;
  dy=t.dy;
}

inline double time_class::doy() {
  return dy;
}

//returns the year field:
inline int time_class::year() {
  return yr;
}

//assignment operator:
inline time_class & time_class::operator = (const time_class &t) {
  yr=t.yr;
  dy=t.dy;
  return *this;
}

//assign from string:
inline time_class & time_class::operator = (const char *t) {
  read_string(t);
  return *this;
}

//assign from floating point.  This gives the number of days elapsed since ref_date.
inline time_class & time_class::operator = (const double t) {
  yr=ref_date.yr;
  dy=ref_date.dy;

  add(t);

  return *this;
}

//equality comparator:
inline int time_class::operator == (const time_class &t2) {
  if (yr == t2.yr && dy == t2.dy) return 1;
  return 0;
}

//non-equality comparator:
inline int time_class::operator != (const time_class &t2) {
  if (yr != t2.yr || dy != t2.dy) return 1;
  return 0;
}

//greater-than operator:
inline int time_class::operator > (const time_class &t2) {
  if (yr > t2.yr) return 1;
  if (yr == t2.yr && dy > t2.dy) return 1;
  return 0;
}

//less-than operator:
inline int time_class::operator < (const time_class &t2) {
  if (yr < t2.yr) return 1;
  if (yr == t2.yr && dy < t2.dy) return 1;
  return 0;
}

//greater-than or equal-to operator:
inline int time_class::operator >= (const time_class &t2) {
  if (yr > t2.yr) return 1;
  if (yr == t2.yr) {
    if (dy > t2.dy || dy == t2.dy) return 1;
  }

  return 0;
}

//less-than or equal-to operator:
inline int time_class::operator <= (const time_class &t2) {
  if (yr < t2.yr) return 1;
  if (yr == t2.yr) {
    if (dy < t2.dy || dy == t2.dy) return 1;
  }

  return 0;
}

//returns the number of days elapsed between this and the argument
//as a floating point value:
inline double time_class::diff(const time_class &t2) const {

  return 365*(yr-t2.yr)+dy-t2.dy+nleap(t2.yr, yr);
}

//converts to a floating point.  Gives the number of days elapsed since ref_date.
inline time_class::operator double() {

  if (yr != 0) return diff(ref_date);
  	else return dy;

//  return 365*(yr-REF_YEAR)+self.dy-REF_DAY-nleap(yr, REF_YEAR);

}

//converts to a string.  Remember that it should be deleted after use!
inline time_class::operator char * () {
  char *time_string;
  time_string=new char[MAX_TSTR_WIDTH];
  return write_string(time_string);
}

inline int time_class::read_string(const char *t, const char *sep) {
  long i=0;
  return read_string(t, i, sep);
}

//io operators:
//std::ostream &operator << (std::ostream &outfile, time_class &t);
//std::istream &operator << (std::istream &infile, const time_class &t);

}

#ifndef _TIME_H
//prevents conflict with C library "time" function, while maintaining
//backwards compatibility with older programs:
//typedef time_class time;
#endif

#endif
