%{
#include <stdio.h>
#include <string.h>
#include <readline/readline.h>
#include <readline/history.h>

extern "C" {int yylex(void);}

#include "sparse_calc_defs.h"

using namespace sparse_calc;
#include "sc_parse.tab.h"

#undef YY_INPUT
#define YY_INPUT(b, r, ms) (r=my_input(b, ms))


int my_input(char *b, int ms) {
  char *line;
  int len=0;

  if (scalc_istream!=NULL) {
    if (fgets(b, ms, scalc_istream)==NULL) {
      fclose(scalc_istream);
      len=0;
    } else {
      len=strlen(b);
      if (sc_stringbuf!=NULL) delete [] sc_stringbuf;
      strncpy(line, b, len);
      sc_stringbuf=line;
      sc_buflen=len;
    }
  } else {

    //interactive interpretter:
    rl_instream=yyin;
    line=readline("\%s\%c>");
    if (feof(yyin)) return YY_NULL;
    if (line == NULL) return YY_NULL;
    len=strlen(line)+1;
    if (len > ms) len=ms;
    strncpy(b, line, len);
    b[len-1]='\n';
    //b[len+1]='\0';
    if (len>1) add_history(line);
    free(line);

  }

  return len;
}

int lineno=1;

extern "C" {
  int yyerror(const char *s);
}

%}

STRING	[a-zA-Z][a-zA-Z0-9\._]*

%%

[\t ]+ ;

gt {
  yylval.id=-SC_N_RESERVED-funtab.lookup(yytext);
  return GT;
};

ge {
  yylval.id=-SC_N_RESERVED-funtab.lookup(yytext);
  return GE;
};

lt {
  yylval.id=-SC_N_RESERVED-funtab.lookup(yytext);
  return LT;
};

le {
  yylval.id=-SC_N_RESERVED-funtab.lookup(yytext);
  return LE;
};

eq {
  yylval.id=-SC_N_RESERVED-funtab.lookup(yytext);
  return EQ;
};

ne {
  yylval.id=-SC_N_RESERVED-funtab.lookup(yytext);
  return NE;
};

{STRING} {
    sc_type_base * var=sc_var_lookup(yytext);
    if (var==NULL) {
      yylval.var=new sc_type_literal(yytext);
      return VARIABLE;
    } else {
      yylval.var=var;
      if (sc_type_of(var) == SUBROUTINE) {
        return SUBROUTINE;
      } else {
        return VARIABLE;
      }
    }
  }
          
};

\"[^\n^\"]+\" {
  yylval.var=new sc_type_literal(yytext);
  return VARIABLE;
};

\( {
  return LEFT_BRAC;
};

\) {
  return RIGHT_BRAC;
};

\, {
  return COMMA;
};

\^ {
  return POW;
}

\[ {
  return LEFT_SUB;
};

\] {
  return RIGHT_SUB;
};

\+ {
  return PLUS;
};

\- {
  return MINUS;
};

\* {
  return PROD;
};

\# {
  return CPROD;
};

\/ {
  return QUOTIENT;
};

\% {
  return MOD;
};

\= {
  return ASSIGN;
};

\| {
  return NORMBRAC;
};

\: {
  return RANGE;
};

[0-9]+ | 
[0-9]+\.[0-9]* |
[0-9]+[eE]\-?[0-9]+ | 
[0-9]+\.[0-9]*[eE]-?[0-9]+ {
  yylval.var=new sc_scalar_t(yytext);
  return VARIABLE;
};

\$_ {
  sc_type_base *var=sc_var_lookup(&sc_state, yytext);
  if (var==NULL) {
    yyerror("Parameter not defined\n");
    return ERROR;
  }
  yylval.var=var;
  return VARIABLE;
};

\n {
  return DELIM;
};
 
%%

int yyerror(const char *s) {
  fprintf(stderr, "sparse_calc, line %d: %s at %s\n", lineno, s, yytext);
}

