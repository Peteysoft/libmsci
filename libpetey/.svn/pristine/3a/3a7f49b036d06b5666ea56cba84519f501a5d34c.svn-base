%{
#include <readline/readline.h>
#include <readline/history.h>

#include "time_class.h"

using namespace std;
using namespace libpetey;

extern "C" {
  int yylex(void);
}

#include "date_calc.tab.h"

extern char *dcinput;
extern int dciptr;

char date2[100];

int lineno=1;

#undef YY_INPUT
#define YY_INPUT(b, r, ms) (r=my_input(b, ms))

int my_input(char *b, int ms) {
  char *line;
  int len=0;

  if (dcinput!=NULL) {
    len=strlen(dcinput+dciptr);
    if (len > ms) len=ms;
    strncpy(b, dcinput+dciptr, len);
    dciptr+=len;
    return len;
  } else {
    rl_instream=yyin;
    line=readline("\%d\%c>");
    if (feof(yyin)) return YY_NULL;
    if (line == NULL) return YY_NULL;
    len=strlen(line)+1;
    if (len > ms) len=ms;
    strncpy(b, line, len);
    b[len-1]='\n';
    //b[len+1]='\0';
    if (len>1) add_history(line);
    free(line);

    return len;
  }
}

extern "C" {
  int yyerror(const char *s);
}

%}

YEAR	[0-9]{1,4}
DMHM	[0-9]{1,2}
TIME	({DMHM}\:{DMHM})|({DMHM}\:{DMHM}\:{DMHM})|({DMHM}\:{DMHM}\:{DMHM}\.[0-9]*)

%%

[\t ]+ ;

help {return HELP;}

{YEAR}\/{DMHM} |
{YEAR}\/{DMHM}\/{DMHM} |
{YEAR}\/{DMHM}\/{DMHM}\-{DMHM} |
{YEAR}\/{DMHM}\/{DMHM}\-{TIME} {
  yylval.date=new time_class(yytext);
  yylval.date->write_string(date2);
  return DATE;
}

{TIME}	{
  strcpy(date2, "0/0/0-");
  strcat(date2, yytext);
  yylval.date=new time_class(date2);
  yylval.date->write_string(date2);
  return DATE;
}

[0-9]+\-{TIME} {
  strcpy(date2, "0/0/");
  strcat(date2, yytext);
  yylval.date=new time_class(date2);
  yylval.date->write_string(date2);
  return DATE;
}

[0-9]+ |
[0-9]+\.[0-9]* {
  yylval.number=atof(yytext);
  return NUMBER;
}

\n {
  lineno++;
  return yytext[0];
}

. return yytext[0];

%%

int yyerror(const char *s) {
  fprintf(stderr, "date_calc, line %d: %s at %s\n", lineno, s, yytext);
}

