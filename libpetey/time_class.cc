//
// Copywrite 2004 Peter Mills.  All rights reserved.
// 
//port of "time_lib.pro" to C++
//

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <time.h>

#include "peteys_tmpl_lib.h"
#include "time_class.h"

#define SPERYEAR 315360000
#define SPERDAY 86400
#define REF_YEAR 1900
#define REF_DAY 0

namespace libpetey {

//global array giving the number of days in each month:
int32_t days_in_month[13]={0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};

//global array giving the total number of days from the beginning
//of the year at the start of each month:
int32_t month_days[13]={0, 0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334};

//returns a 1 if year is leap, 0 if not:
inline int isleap(int year) {
  if (year % 400 == 0) return 1;
  if (year % 100 == 0) return 0;
  if (year % 4 == 0) return 1;
  return 0;
  //return (year/4 == year/4.0);
}

//reference date for converting dates into floating point:
time_class ref_date(1900, 1, 1, 0, 0, 0);

#define MAXFIELDSIZE 100

//current time:
void time_class::now() {
  time_t timer;
  struct tm *t2;

  time(&timer);

  t2=localtime(&timer);

  init(t2->tm_year+1900, t2->tm_mon+1, t2->tm_mday, t2->tm_hour, t2->tm_min, t2->tm_sec);

  free(t2);
}

//Converts a date from a string.  Must be in the following format:
//	year/month/day-hour:minute:second
//Whitespace is allowed.

int time_class::read_string(const char *t, long &i, const char *sep) {
  //char *sep=SEPARATORS;			//list of separators
  int32_t field[MAXFIELDSIZE];		//string containing the current field
  int i1, j, isep;			//counters
  int f[5]={0, 0, 0, 0, 0};		//year, month, day, hour, minute fields
  float s=0.0;				//seconds field
//  int days_in_month[12]={31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
  int oldyear;
  double dd;
  int y, m, d;

  i=0;

  for (isep=0; sep[isep]!=0; isep++) {
    for (j=i; t[j]!=0 && t[j]!=sep[isep];j++) {
      //printf("%c%5i\n", t[j], j);
      if ((t[j] < 48 || t[j] > 57) && t[j]!=32 && t[j]!=9) {
//        printf("Bad character in date string\n");
        if (j > i) sscanf(&t[i], "%d%n", &f[isep], &i1);
        init(f[0], f[1], f[2], f[3], f[4], s);
        i=i+i1;
        return 1;
//      } else {
//        field[j-i]=t[j];
      }
/*      if (j >= MAXFIELDSIZE) {
        printf("Field width exceeds maximum\n");
	return 2;
      }*/
    }
    sscanf(&t[i], "%d", &f[isep]);
//    field[j-i]=0;
//    printf("%d\n", f[isep]);
//    sscanf(field, "%d", &f[isep]);

    if (t[j] == 0) break;
    i=j+1;
  }

  //printf("%s %f\n", t[j], s);
  if (t[j] == sep[4]) sscanf(&t[i], "%f%n", &s, &i1);
  i=i+i1;

  //done parsing the string, here we initialize the fields:

  return init(f[0], f[1], f[2], f[3], f[4], s);
}


//Converts a date from a string.  Must be in the following format:
//	year/month/day-hour:minute:second
//Whitespace is allowed.
/*
int time_class::read_string(const char *t, long &i) {
  char *sep=SEPARATORS;			//list of separators
  char field[MAXFIELDSIZE];		//string containing the current field
  int i1, j, isep;			//counters
  int f[5]={0, 0, 0, 0, 0};		//year, month, day, hour, minute fields
  float s=0.0;				//seconds field
//  int days_in_month[12]={31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
  int oldyear;
  double dd;
  int y, m, d;

  i=0;

  for (isep=0; sep[isep]!=0; isep++) {
    //skip leading whitespace:
    for (j=1; t[j] == 32 || t[j] == 9; j++) {
      if (t[j]==0) {
        i=j;
        init(f[0], f[1], f[2], f[3], f[4], s);
      }
    }
    //get the "meat":
    for (j=i; t[j]!=0 && t[j]!=sep[isep] && t[j] != 32 && t[j] !=9;j++) {
      //printf("%c%5i\n", t[j], j);
      if (t[j] < 48 || t[j] > 57) {
//        printf("Bad character in date string\n");
        if (j > i) sscanf(&t[i], "%d%n", &f[isep], &i1);
        init(f[0], f[1], f[2], f[3], f[4], s);
        i=i+i1;
        return 1;
//      } else {
//        field[j-i]=t[j];
      }
      if (j >= MAXFIELDSIZE) {
        printf("Field width exceeds maximum\n");
	return 2;
      }
    }
    //scan trailing whitespace:
    if (t[j]=32 || t[j] == 9) {
      for (j=i; t[j]!=0 && t[j]!=sep[isep];j++) {
        if (t[j] != 32 && t[j] != 9) {
          if (j > i) sscanf(&t[i], "%d%n", &f[isep], &i1);
          init(f[0], f[1], f[2], f[3], f[4], s);
          i=i+i1;
          return 1;
        }
      }
    }
    sscanf(&t[i], "%d", &f[isep]);
//    field[j-i]=0;
//    printf("%d\n", f[isep]);
//    sscanf(field, "%d", &f[isep]);

    if (t[j] == 0) break;
    i=j+1;
  }

  //printf("%s %f\n", t[j], s);
  if (t[j] == sep[4]) sscanf(&t[i], "%f%n", &s, &i1);
  i=i+i1;

  //done parsing the string, here we initialize the fields:

  return init(f[0], f[1], f[2], f[3], f[4], s);
}
*/

//Given year, month, day, hour, minute, second, initialized an
//instance of the time class:

int time_class::init(	int y, 		//the year--no non-y2k compliant, please
		int mon, 	//the month
		int d, 		//the day of the month
		int h, 		//the hour
		int min, 	//minutes
		float s) {	//seconds

//  int month_days[12]={0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334};
  double dd;

  if (y==0) {
    yr=0;
    dy=mon*30+d+((s/60+min)/60.+h)/24;
    return 0;
  }

  if (d != 0) d=d-1;

  yr=y+(mon-1) / 12;
  mon=(mon-1) % 12 + 1;

  d=d+month_days[mon];
/*
  for (long i=0; i<mon; i++) {
    d=d+days_in_month[i];
  }
*/
//  printf("Year, month, day, hour, minute, second: %d %d %d %d %d %f\n",
//  		yr, m, d, f[3], f[4], s);
  dd=d+((s/60.+min)/60.+h)/24.;
//  printf("Day of year=%f\n", dd);

  dy=0;
  add(dd);
//  printf("Day of year=%f\n", dy);

  if (isleap(yr) && mon>2) dy++;


  return 0;

}

// Writes an object instance of type time to a string:

char * time_class::write_string(char * time_string, const char *sep) {
  //char *sep=SEPARATORS;
//  int month_days[12]={0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334};
  short int y, mon, d, h, min;
  short int dy_int;
  float s, dayfrac;

  get_fields(y, mon, d, h, min, s);

  if (y == 0 && mon == 0 && d == 0) {
    if (s == 0) {
      sprintf(time_string, "%2.2d%1c%2.2d", h, sep[3], min);
    } else {
      sprintf(time_string, "%2.2d%1c%2.2d%1c%5.2f", h, sep[3], min, sep[4], s);
    }
  } else {
    if (s == 0) {
      if (min == 0) {
        if (h == 0) {
	  if (d == 0) {
	    sprintf(time_string, "%4.4d%1c%2.2d", y, sep[0], mon);
	  } else {
	    sprintf(time_string, "%4.4d%1c%2.2d%1c%2.2d", y, sep[0], mon, sep[1], d);
	  }
	} else {
	  sprintf(time_string, "%4.4d%1c%2.2d%1c%2.2d%1c%2.2d", y, sep[0], mon, sep[1], d, sep[2], h);
	}
      } else {
	sprintf(time_string, "%4.4d%1c%2.2d%1c%2.2d%1c%2.2d%1c%2.2d", y, sep[0], mon, sep[1], d, sep[2], h, sep[3], min);
      }
    } else {
      sprintf(time_string, "%4.4d%1c%2.2d%1c%2.2d%1c%2.2d%1c%2.2d%1c%5.2f",
      		y, sep[0], mon, sep[1], d, sep[2], h, sep[3], min, sep[4], s);
    }
  }

  return time_string;
}

// Adds d days to the object instance:

time_class & time_class::add(double d) {
  int oldyear;

  dy=dy+d;
  if (yr != 0 && dy >= 365 + isleap(yr)) {
    oldyear=yr;
    yr=yr+(int) dy / 365;
    dy=dy-365.*(int) (dy / 365.);
    dy=dy-nleap(oldyear, yr);
  }
  while (dy < 0) {
    yr--;
    dy=365+dy+isleap(yr);
  }

  return *this;

}

// Addition operator.  May not be well defined if the second argument has
// a year field.

time_class time_class::operator + (time_class &t2) {
  time_class newt;
  char *this_str, *t2_str;

  if (dy !=0 && yr !=0 && t2.yr != 0 && t2.dy != 0) {
    this_str=(char *) *this;
    t2_str=(char *) t2;
    printf("Addition of two times (%s, %s) is ambiguous\n", this_str, t2_str);
    delete this_str;
    delete t2_str;
  }

  newt.yr=yr+t2.yr;
  newt.dy=dy;

  newt.add(t2.dy);

  return newt;
}

// Subtraction operator.  The returned time will not have a year field.

time_class time_class::operator - (const time_class &t2) {
  time_class newt;

  if (t2.yr != 0) {
    newt.yr=0;
    newt.dy=diff(t2);
    return newt;
  } else {
    newt.yr=yr;
    newt.dy=dy;
    newt.add(-t2.dy);
  }

  return newt;
}

// Multiplication operator.  May not be well defined if the first argument
// has a year field.

time_class time_class:: operator * (double n) {
  double s;
  time_class newt;
  char *this_str;

  if (yr != 0 && dy !=0) {
    this_str=(char *) *this;
    printf("Multiplication of time and double (%s * %g) is ambiguous\n", (char *) this_str, n);
    newt.dy=((double) *this)*n;
    delete this_str;
  } else if (yr==0) {
    newt.dy=this->dy*n;
  } else {
    newt.dy=this->yr*n;
  }

  return newt;
}  

// Division operator.  May not be well defined if the first argument has a
// year field.

time_class time_class::operator / (double n) {
  double s;
  time_class newt;
  char *this_str;

  if (yr != 0 && dy !=0) {
    this_str=(char *) *this;
    printf("Division of time and double (%s / %g) is ambiguous\n", this_str, n);
    newt.dy=((double) *this)/n;
    delete this_str;
  } else if (yr==0) {
    newt.dy=this->dy/n;
  } else {
    newt.dy=this->yr/n;
  }
  
  return newt;
}

// Division operator.  May not be well defined if the second field has a
// year field.

double time_class::operator / (time_class &t2) {
  double n;
  char *this_str, *t2_str;

  if ((yr != 0 || t2.yr !=0) && (dy != 0 || t2.dy != 0)) {
    this_str=(char *) *this;
    t2_str=(char *) t2;
    printf("Division of time and time (%s / %s) is ambiguous\n", this_str, t2_str);
    n=((double) *this)/((double) t2);
    delete this_str;
    delete t2_str;
  } else if (yr==0) {
    n=this->dy/t2.dy;
  } else {
    n=this->yr/t2.yr;
  }

  return n;
}

// Returns all the fields of the date.  Two versions--one for longword and one for
// short.
template <class integer, class real>
void time_class::get_fields(integer &y, 	//year
			integer &mon, 	//month
			integer &d,	//day of the month
			integer &h, 	//hour
			integer &min, 	//minutes past the hour
			real &s) {	//seconds

//  int month_days[12]={0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334};
  real dayfrac;
  integer dy_int;

  y=yr;

  dayfrac=dy-(integer) dy;
  dayfrac=dayfrac*24.;
  h=(integer) dayfrac;
  dayfrac=(dayfrac-h)*60;
  min=(integer) dayfrac;
  s=(dayfrac-min)*60;

  if (y == 0) {
    d=(integer) dy;
    //printf("Day: %d %f\n", d, dy);
    mon=0;
  } else {
    //do a binary search on the month:

    dy_int=(integer) dy;

    mon=(integer) (bin_search(month_days+1, 12L, (int32_t) dy_int, -1L))+1;

    d=1;
    if (isleap(yr)) {
      if (mon > 2) {
        if (dy_int == month_days[mon]) {
          mon--;
	}
	if (dy_int != 59) d=0;
      }
/*      if (dy_int == 59) {
        mon--;
	d=0;
      }
*/
    }
    d+=dy_int-month_days[mon];
  }
}

template void time_class::get_fields<int16_t, float>(int16_t &, int16_t &, int16_t &,
		int16_t &, int16_t &, float &);

template void time_class::get_fields<int32_t, float>(int32_t &, int32_t &, int32_t &,
		int32_t &, int32_t &, float &);
template void time_class::get_fields<int32_t, double>(int32_t &, int32_t &, int32_t &,
		int32_t &, int32_t &, double &);

template void time_class::get_fields<int64_t, float>(int64_t &, int64_t &, int64_t &,
		int64_t &, int64_t &, float &);
template void time_class::get_fields<int64_t, double>(int64_t &, int64_t &, int64_t &,
		int64_t &, int64_t &, double &);

// Returns the month field:

int time_class::month() {
  int mon, dy_int;
//  int month_days[12]={0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334};
  float dayfrac;

  dy_int=(int) dy;

  if (yr == 0) {
    mon=0;
  } else {
    //do a binary search on the month:
    mon=(int) bin_search(month_days+1, 12L, (int32_t) dy, -1L)+1;

    if (isleap(yr)) {
      if (mon > 2 && dy_int == month_days[mon]) {
        mon--;
      }
//      if (dy_int==60) mon--;
    }
  }

  return mon;
}

// Returns the day-of-the-month field:

int time_class::day() {
  int mon, d, dy_int;
//  int month_days[12]={0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334};

  if (yr == 0) {
    d=(int) dy;
  } else {
    //do a binary search on the month:
    dy_int=(int) dy;

    mon=(int) (bin_search(month_days+1, 12L, (int32_t) dy_int, -1L))+1;

    d=1;
    if (isleap(yr)) {
      if (mon > 2) {
        if (dy_int == month_days[mon]) {
          mon--;
	}
	if (dy_int != 59) d=0;
      }
      /*
      if (dy_int == 60) {
        mon--;
	d=0;
      }
      */
    }
    d+=dy_int-month_days[mon];

  }

  return d;
}

//returns the day of the week:
int time_class::dow() {
  double dd;

  dd=(double) *this;

  if (dd >= 0) {
    //add one because 1900 was not a leap year
    //(and this program thinks it was...)
    //(**--fixed)
    return ((long(dd))%7)+1;
  } else {
    return 7+long(floor(dd))%7;
  }

}

// Returns the hour field:

int time_class::hour() {
  int h;
  float dayfrac;

  dayfrac=dy-(int) dy;
  dayfrac=dayfrac*24.;
  h=(int) dayfrac;

  return h;
}

// Returns the minutes field:

int time_class::minute() {
  int h, min;
  float dayfrac;

  dayfrac=dy-(int) dy;
  dayfrac=dayfrac*24.;
  h=(int) dayfrac;
  dayfrac=(dayfrac-h)*60;
  min=(int) dayfrac;

  return min;
}

// Returns the seconds field:

float time_class::second() {
  float min;
  float dayfrac, s;

  min=(dy-(int) dy)*24*60;
  s=(min-(int) min)*60;

  return s;

}

time_class & time_class::second(float s) {
  float sold;
  sold=second();

  add((s-sold)/SPERDAY);

  return *this;
}


/*
//Defines the output operator:
ostream & operator << (ostream &outfile, time_class &t) {
  char tstring[32];
  t.write_string(tstring);
  outfile << tstring;
  return outfile;
}
*/

}  //end namespace libpetey

