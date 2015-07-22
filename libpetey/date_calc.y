%{
#include <stdio.h>
#include <string.h>
#include <time.h>

#include "time_class.h"

#define DATESTR_LEN 100

using namespace std;
using namespace libpetey;

extern "C" {
  int yylex(void);
  int yyerror(char *);
}

char datestr[DATESTR_LEN];

char *dateformat;
char *dcinput;
int dciptr;

void date_calc_helpscreen(FILE *docfs) {
  fprintf(docfs, "\n");
  fprintf(docfs, "SYNOPSIS\n");
  fprintf(docfs, "\n");
  fprintf(docfs, "date_calc [-f [\"]format[\"]] [-h] [[\"]commands[\"]]\n");
  fprintf(docfs, "\n");
  fprintf(docfs, "DESCRIPTION\n");
  fprintf(docfs, "\n");
  fprintf(docfs, "The date calculator, date_calc, is an interactive utility for performing\n");
  fprintf(docfs, "calculations with dates and times just as you would with floating point\n");
  fprintf(docfs, "numbers.  date_calc also handles ordinary floating point calculatins\n");
  fprintf(docfs, "\n");
  fprintf(docfs, "OPTIONS\n");
  fprintf(docfs, "\n");
  fprintf(docfs, "-f [\"]format[\"]\n");
  fprintf(docfs, "   Specifies the output format for calculations.  Uses the same format\n");
  fprintf(docfs, "   codes as the date command\n");
  fprintf(docfs,  "\n");
  fprintf(docfs, "-h\n");
  fprintf(docfs, "  Print this help screen\n");
  fprintf(docfs, "\n");
  fprintf(docfs, "INPUT FORMAT\n");
  fprintf(docfs, "\n");
  fprintf(docfs, "Dates must be in the following format, as specified using the date\n");
  fprintf(docfs, "style format codes:\n");
  fprintf(docfs, "\n");
  fprintf(docfs, "%%Y/%%m[/%%d[-%%H[:%%M[:%%S]]]]\n");
  fprintf(docfs, "\n");
  fprintf(docfs, "where the seconds field (%%S) can be expressed as a decimal.\n");
  fprintf(docfs, "Alternatively:\n");
  fprintf(docfs, "\n");
  fprintf(docfs, "%%H:%%M[:%%S]\n");
  fprintf(docfs, "\n");
  fprintf(docfs, "OPERATORS\n");
  fprintf(docfs, "\n");
  fprintf(docfs, "The date calculator accepts the four arithmetic operators.  Addition\n");
  fprintf(docfs, "and multiplication are given by the usual symbols: + and *,\n");
  fprintf(docfs, "respectively, however, because the symbols normally reserved for\n");
  fprintf(docfs, "subtraction and division are already used to specify the dates\n");
  fprintf(docfs, "we use instead: _ and |, respectively.  Note that multiplication\n");
  fprintf(docfs, "must combine either a date and a number or two numbers.  It follows\n");
  fprintf(docfs, "then that if you divide a date by a date, you get back a number\n");
  fprintf(docfs, "while dividing a date with a number returns a date.\n");
  fprintf(docfs, "\n");
  fprintf(docfs, "EXAMPLE\n");
  fprintf(docfs, "\n");
  fprintf(docfs, "I was born on December 22, 1973 at eight o'clock in the morning.\n");
  fprintf(docfs, "How many hours have I been alive?\n");
  fprintf(docfs, "\n");
  fprintf(docfs, "$ date_calc \"(2013/8/1_1973/12/22-8)|1:0\"\n");
  fprintf(docfs, "\n");
  fprintf(docfs, "which returns:\n");
  fprintf(docfs, "\n");
  fprintf(docfs, "347200\n");
  fprintf(docfs, "\n");
  fprintf(docfs, "I have some climate data starting on January 1, 1948, with one field\n");
  fprintf(docfs, "for every six hours.  I need the date of the 20165th field (zero-based)\n");
  fprintf(docfs, "and I need to print it out in a format convenient for turning into\n");
  fprintf(docfs, "a file name:\n");
  fprintf(docfs, "\n");
  fprintf(docfs, "$ date_calc -f \"%%Y.%%m.%%d.%%H\" \"1948/1/1+6:0*20165\"\n");
  fprintf(docfs, "\n");
  fprintf(docfs, "1961.10.20.06\n");
  fprintf(docfs, "\n");
  fprintf(docfs, "SEE ALSO: date\n");
  fprintf(docfs, "\n");
  fprintf(docfs, "AUTHOR: Peter Mills (petey@peteysoft.org)\n");
  fprintf(docfs, "\n");
}

%}

%union {
  time_class *date;
  double number;
}

%token <date> DATE
%token <number> NUMBER
%token HELP

%left '_' '+'
%left '*' '|'

%type <date> date_exp
%type <number> number_exp

%%
statement_list: statement
	|     statement_list statement ;

statement: date_exp '\n' {
    struct tm cptm;
    float sec;

    if (dateformat==NULL) {
	 $1->write_string(datestr);
	printf("%s\n", datestr);
	delete $1;
    } else {
      $1->get_fields(cptm.tm_year, cptm.tm_mon, cptm.tm_mday, cptm.tm_hour,
		cptm.tm_min, sec);
      cptm.tm_year-=1900;
      cptm.tm_mon--;
      cptm.tm_wday=$1->dow();
      cptm.tm_yday=$1->doy();
      cptm.tm_sec=sec;
      cptm.tm_isdst=0;
      strftime (datestr, DATESTR_LEN, dateformat, &cptm);
      printf("%s\n", datestr);
      delete $1;
    }
  }
| number_exp '\n' {printf("%lg\n", $1);}
| HELP '\n' {date_calc_helpscreen(stdout);}
| error '\n' {
    yyclearin;
    yyerrok;
  }
| '\n';

date_exp:  
 date_exp '+' date_exp { $$ = new time_class(*$1 + *$3);
				delete $1; delete $3; }
	| date_exp '_' date_exp { $$ = new time_class(*$1 - *$3); 
				delete $1; delete $3;} 
	| date_exp '*' number_exp { $$ = new time_class(*$1 * $3);
				delete $1;} 
	| number_exp '*' date_exp { $$ = new time_class(*$3 * $1);
				delete $3;}
	| date_exp '|' number_exp { $$ = new time_class(*$1/$3);
				delete $1;} 
	| '(' date_exp ')' { $$ = $2; }
 	| DATE {$$=$1;};

number_exp: 
	 number_exp '+' number_exp { $$ = $1 + $3; }
	| number_exp '_' number_exp { $$ = $1 - $3; }
	| number_exp '*' number_exp { $$ = $1 * $3;}
	| date_exp '|' date_exp { $$ = *$1 / *$3;
				delete $1; delete $3;} 
	| number_exp '|' number_exp { $$ = $1 / $3;}
	| '(' number_exp ')' { $$ = $2; }
	| NUMBER {$$=$1;};
%%

#include "parse_command_opts.h"

int main(int argc, char **argv) {
  void * optarg[20];
  int optflag[20];
  int dcilen;

  dateformat=NULL;
  argc=parse_command_opts(argc, argv, "fh", "%s%", optarg, optflag, 1);
  if (argc<0) {
    fprintf(stderr, "date_calc: error parsing command line\n");
    fprintf(stderr, "           use -h for help\n");
    return -1;
  }
  if (optflag[1]) {
    date_calc_helpscreen(stdout);
    return 0;
  }
    
  if (optflag[0]) dateformat=(char *) optarg[0];

  if (argc>1) {
    dcilen=0;
    for (int i=1; i<argc; i++) dcilen+=strlen(argv[i]);
    dcinput=new char [dcilen+argc];
    dcinput[0]='\0';
    for (int i=1; i<argc; i++) {
      strcat(dcinput, argv[i]);
      strcat(dcinput, "\n");
    }
    dciptr=0;
  } else {
    dcinput=NULL;
  }

  yyparse();
}

