%{
#include <stdio.h>
#include <string.h>
#include <readline/readline.h>
#include <readline/history.h>

extern "C" {int yylex(void);}

#include "sparse_calc_defs.h"

using namespace sparse_calc;
#include "sparse_calc.tab.h"

#undef YY_INPUT
#define YY_INPUT(b, r, ms) (r=my_input(b, ms))


int my_input(char *b, int ms) {
  char *line;
  int len=0;

  if (scalc_istream!=NULL) {
    if (fgets(b, ms, scalc_istream)==NULL) {
      fclose(scalc_istream);
      if (sc_stream_stack.pop(scalc_istream)<0) {
        if (sc_scriptflag) return 0;
        scalc_istream=NULL;
      }
      strcpy(b, "\n");
      //if running script from command line, end parsing:
      return 1;
    }
    len=strlen(b);
    if (sc_firstpass==1) {
      if (sc_stringbuf!=NULL) delete [] sc_stringbuf;
      strncpy(line, b, len);
      sc_stringbuf=line;
      sc_buflen=len;
    } else {
      delete [] line;
    }
    return len;
  }

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

subroutine {return SUBROUTINE_T;}
vector {return VECTOR_T;}
full {
    yylval.id=SC_FULL_T;
    return MATRIX_T;
  };
sparse {
    yylval.id=SC_SPARSE_T;
    return MATRIX_T;
  };
sparse_array {
    yylval.id=SC_SPARSE_ARRAY_T;
    return MATRIX_T;
  };

gt {
  yylval.id=-SC_N_RESERVED-funtab.lookup(yytext);
  return COMPARATOR;
};

ge {
  yylval.id=-SC_N_RESERVED-funtab.lookup(yytext);
  return COMPARATOR;
};

lt {
  yylval.id=-SC_N_RESERVED-funtab.lookup(yytext);
  return COMPARATOR;
};

le {
  yylval.id=-SC_N_RESERVED-funtab.lookup(yytext);
  return COMPARATOR;
};

eq {
  yylval.id=-SC_N_RESERVED-funtab.lookup(yytext);
  return COMPARATOR;
};

ne {
  yylval.id=-SC_N_RESERVED-funtab.lookup(yytext);
  return COMPARATOR;
};

^[ \t]*{STRING} {
  id=funtab.lookup(yytext);
  if (id>=0) {
    yylval.id=-SC_N_RESERVED-id;
    return SUBROUTINE;
  }
  yylval.symbol=new char[strlen(yytext)+1];
  strcpy(yylval.symbol, yytext);
  return SYMBOL;
}

{STRING} {
  FILE *fs;
  char **line;
  int id;
  int err;
  long n;
  char *fname;
  string_petey symbol(yytext);
  time_t moddate;
  int vt;

  id=funtab.lookup(symbol);
  if (id>=0) {
    yylval.id=-SC_N_RESERVED-id;
    return SUBROUTINE;
  }

  id=symtab.lookup(symbol);
  //printf("id of the symbol, %s, is %d\n", yytext, id);

  if (id < 0) {
    yylval.symbol=new char[strlen(yytext)+1];
    strcpy(yylval.symbol, yytext);
    return SYMBOL;
  } else {
    vt=abs(vartyp[id]);
    if (vt==SC_SPARSE_T || 
		vt==SC_SPARSE_ARRAY_T || 
		vt==SC_FULL_T) {
      yylval.id=id;
      delete [] fname;
      return MATRIX;
    } else if (vt==SC_VECTOR_T) {
      yylval.id=id;
      return VECTOR;
    } else if (vt==SC_SCALAR_T) {
      yylval.scalar=scaltab[id];
      return SCALAR;
    } else {
      yylval.symbol=new char[strlen(yytext)+1];
      strcpy(yylval.symbol, yytext);
      return SYMBOL;
    }
  }
          
};

\"[^\n^\"]+\" {
  int len;
  len=strlen(yytext)-1;
  //copy, stripping the quotation marks:
  yylval.literal=new char[len];
  strncpy(yylval.symbol, yytext+1, len-1);
  yylval.literal[len-1]='\0';
  return LITERAL;
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
  yylval.id=-SC_N_RESERVED-funtab.lookup(yytext);
  return POW;
}

\[ {
  return LEFT_SUB;
};

\] {
  return RIGHT_SUB;
};

\+ {
  yylval.id=-SC_N_RESERVED-funtab.lookup(yytext);
  return PLUS;
};

\- {
  yylval.id=-SC_N_RESERVED-funtab.lookup(yytext);
  return MINUS;
};

\* {
  yylval.id=-SC_N_RESERVED-funtab.lookup(yytext);
  return PROD;
};

\# {
  yylval.id=-SC_N_RESERVED-funtab.lookup(yytext);
  return CPROD;
};

\/ {
  yylval.id=-SC_N_RESERVED-funtab.lookup(yytext);
  return QUOTIENT;
};

\% {
  yylval.id=-SC_N_RESERVED-funtab.lookup(yytext);
  return MOD;
};

\= {
  return ASSIGN;
};

\| {
  return NORMBRAC;
};

\: {
  yylval.id=funtab.lookup(yytext);
  return RANGE;
};

[0-9]+ | 
[0-9]+\.[0-9]* |
[0-9]+[eE]\-?[0-9]+ | 
[0-9]+\.[0-9]*[eE]-?[0-9]+ {
  yylval.scalar=atof(yytext);
  return CONSTANT;
};

\$[0-9]+ {
  long id;
  id=symtab.lookup(yytext);
  if (id<0) {
    id=symtab.add(yytext);
    vartyp[id]=SC_PARAMETER_T;
  }
  return PARAMETER;
};

\n {
  return DELIM;
};
 
%%

int yyerror(const char *s) {
  fprintf(stderr, "sparse_calc, line %d: %s at %s\n", lineno, s, yytext);
}

