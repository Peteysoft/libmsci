%{
#include "time_class.h"
#include "linked.h"

using namespace std;
using namespace libpetey;

linked_list<time_class> date_stack;
linked_list<double> number_stack;
linked_list<long> type_stack;

enum token_type {DATE, NUMBER, OPERATOR, ERROR};

int push_date (char *date) {
  token_type ttype=DATE;
  time_class time;
  char test[30];

  time.read_string(date);
  date_stack.push(time);
  type_stack.push(ttype);

  time.write_string(test);
  printf("%s\n", test);
  return ttype; 
}

int push_number (char *number) {
  token_type ttype=NUMBER;
  double value;

  sscanf(number, "%lg", &value);
  number_stack.push(value);
  type_stack.push(ttype);
  return ttype; 
}

int calculate(char *op) {
  long t1;
  long t2;
  long n;
  char rstring[30];

  n=type_stack.pop(t2);
  if (n <= 0) {
    fprintf(stderr, "Only %d elements on stack\n", n+1);
    type_stack.push(t2);
    return ERROR;
  }
  type_stack.pop(t1);

  if ((t1 == DATE && t2 == NUMBER) || (t1 == NUMBER && t2 == DATE)) {
    time_class date;
    double number;
    time_class result;
    date_stack.pop(date);
    number_stack.pop(number);

    switch (op[0]) {
      case '+':
        fprintf(stderr, "Cannot add a date and a number\n");
        return ERROR;
      case '-':
        fprintf(stderr, "Cannot subtract a number from a date (or vice versa)\n");
        return ERROR;
      case '*':
        result=date*number;
        result.write_string(rstring);
        printf("%s\n", rstring);
        date_stack.push(result);
        type_stack.push(DATE);
        return OPERATOR;
      case '/':
        result=date/number;
        result.write_string(rstring);
        printf("%s\n", rstring);
        date_stack.push(result);
        type_stack.push(DATE);
        if (t1==NUMBER && t2==DATE) {
          fprintf(stderr, "Warning: swapped operands (date with number)\n");
          return ERROR;
        }
        return OPERATOR;
    }
  } else if (t1 == DATE && t2 == DATE) {
    time_class date1, date2;
    time_class result1;
    double result2;
    date_stack.pop(date2);
    date_stack.pop(date1);

    switch (op[0]) {
      case '+':
        result1=date1+date2;
        result1.write_string(rstring);
        printf("%s\n", rstring);
        date_stack.push(result1);
        type_stack.push(DATE);
        return OPERATOR;
      case '-':
        result1=date1-date2;
        result1.write_string(rstring);
        printf("%s\n", rstring);
        date_stack.push(result1);
        type_stack.push(DATE);
        return OPERATOR;
      case '*':
        fprintf(stderr, "Cannot multiply two dates\n");
        return ERROR;
      case '/':
        result2=date1/date2;
        printf("%lg\n", result2);
        number_stack.push(result2);
        type_stack.push(NUMBER);
        return OPERATOR;
    }
  } else {
    double value1;
    double value2;
    double result;
    number_stack.pop(value2);
    number_stack.pop(value1);
    switch (op[0]) {
      case '+':
        result=value1+value2;
        printf("%lg\n", result);
        return OPERATOR;
      case '-':
        result=value1-value2;
        printf("%lg\n", result);
        return OPERATOR;
      case '*':
        result=value1*value2;
        printf("%lg\n", result);
        return OPERATOR;
      case '/':
        result=value1/value2;
        printf("%lg\n", result);
        return OPERATOR;
    }
  }
  return ERROR;
}

int push_date2(char *date) {
  char *date2;
  int err;

  date2=new char[strlen(date)+7];
  strcpy(date2, "0/0/0-");
  strcat(date2, date);
  err=push_date(date2);
  delete [] date2;

  return err;
}

int push_date3(char *date) {
  char *date2;
  int err;

  date2=new char[strlen(date)+7];
  strcpy(date2, "0/0/");
  strcat(date2, date);
  err=push_date(date2);
  delete [] date2;

  return err;
}

%}

YEAR	[0-9]{1,4}
DMHM	[0-9]{1,2}
TIME	({DMHM}\:{DMHM})|({DMHM}\:{DMHM}\:{DMHM})|({DMHM}\:{DMHM}\:{DMHM}\.[0-9]*)

%%

[\t ]+ ;

[0-9]+ |
[0-9]+\.[0-9]* push_number(yytext);

{YEAR}\/{DMHM} |
{YEAR}\/{DMHM}\/[0-9]{1,2} |
{YEAR}\/{DMHM}\/{DMHM}\-{DMHM} |
{YEAR}\/{DMHM}\/{DMHM}\-{TIME} push_date(yytext);

{TIME}	push_date2(yytext);

[0-9]+\-{TIME} push_date3(yytext);

[+-/*] calculate(yytext);
%%

int main() {
  printf("RPN style date calculator\n");
  yylex();
}

