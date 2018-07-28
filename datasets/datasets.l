%{
#include <stdio.h>
#include <string.h>

extern "C" {int yylex(void);}

#include "datasets.tab.h"

int lineno=1;

extern "C" {
  int yyerror(const char *s);
}

%}

CD_DEL \.
STATEMENT_DELIM ;
STRING	[a-zA-Z][a-zA-Z0-9_]*

%%

{STRING} {
  yylval.symbol=new char[strlen(yytext)+1];
  strcpy(yylval.symbol, yytext);
  return SYMBOL;
};

[\t ]+ ;

\-?[0-9]+ | 
\-?[0-9]+\.[0-9]* |
[0-9]+[eE]\-?[0-9]+ |
[0-9]+\.[0-9]*[eE]-?[0-9]+ {
  yylval.scalar=atof(yytext);
  return SCALAR;
};

\n {
  lineno++; 
}

{STATEMENT_DELIM} {
  lineno++; 
  return STMT_DELIM;
}

{CD_DEL} {
  return CD_DELIM;
}

[+-/*] {
  yylval.op=new char[strlen(yytext)+1];
  strcpy(yylval.op, yytext);
  return INFIX_OP;
}

\=\(\)\, {
  return yytext[0];
};
 
%%

int yyerror(const char *s) {
  fprintf(stderr, "%d: %s at %s\n", lineno, s, yytext);
}

