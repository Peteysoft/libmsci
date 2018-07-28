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

print {return PRINT;};
save {return SAVE;}
path {return PATH;};
help {return HELP;};
sys {return SYSTEM;};
com {return COMMENT;};
run {return RUN;};
size {return SIZE;};

delete {return DELETE;};

scalar {return SCALAR_T;}
vector {return VECTOR_T;}
full {
    yylval.type_id=FULL_T;
    return MATRIX_T;
  };
sparse {
    yylval.type_id=SPARSE_T;
    return MATRIX_T;
  };
sparse_array {
    yylval.type_id=SPARSE_ARRAY_T;
    return MATRIX_T;
  };

^[ \t]*{STRING} {
  yylval.symbol=new char[strlen(yytext)+1];
  strcpy(yylval.symbol, yytext);
  return SYMBOL;
}

all {return ALL;}

transpose {return TRANSPOSE;};

identity {return IDENTITY;};

evd {return EVD;};

svd {return SVD;};

solve {return SOLVE;};

invert {return INVERT;};

empty {return EMPTY;};

read {return READ;};

sqrt {return SQRT;}
sin {return SIN;}
cos {return COS;}
tan {return TAN;}
log {return LOG;}

{STRING} {
  FILE *fs;
  char **line;
  int id;
  int err;
  long n;
  char *fname;
  string_petey symbol(yytext);
  time_t moddate;

  id=symtab.lookup(symbol);
  //printf("id of the symbol, %s, is %d\n", yytext, id);

  if (id < 0) {
    yylval.symbol=new char[strlen(yytext)+1];
    strcpy(yylval.symbol, yytext);
    return SYMBOL;
  } else if (vartyp[id]==SPARSE_T || 
		vartyp[id]==SPARSE_ARRAY_T || 
		vartyp[id]==FULL_T) {
    fname=add_path(yytext);
    yylval.matrix=read_matrix(fname, vartyp[id], BIN);
    if (yylval.matrix==NULL) return ERROR;
    delete [] fname;
    return MATRIX;
  } else if (vartyp[id]==VEC_T) {
    fname=add_path(yytext);
    if (sc_firstpass==1) {
      moddate=fmoddate(fname);
      if (moddate>sc_mod1 || moddate==-1) {
        delete [] fname;
        yylval.vector.data=NULL;
        yylval.vector.n=0;
        return VEC;
      }
    }
    err=read_vector(fname, &yylval.vector, BIN);
    delete [] fname;
    if (err==0) return VEC; else return ERROR;
  } else if (vartyp[id]==SCALAR_DATA) {
    yylval.scalar=scaltab[id];
    return SCALAR;
  } else {
    yylval.symbol=new char[strlen(yytext)+1];
    strcpy(yylval.symbol, yytext);
    return SYMBOL;
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
  yylval.scalar=atof(yytext);
  return SCALAR;
};

\n {
  if (sc_firstpass!=1) lineno++; 
  return DELIM;
};
 
%%

int yyerror(const char *s) {
  fprintf(stderr, "sparse_calc, line %d: %s at %s\n", lineno, s, yytext);
}

