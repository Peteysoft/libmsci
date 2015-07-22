%{


extern "C" {int yylex(void);}

#undef YY_INPUT
#define YY_INPUT(b, r, ms) (r=my_input(b, ms))

int my_input(char *b, int ms) {

}

int lineno=1;

extern "C" {
  int yyerror(const char *s);
}

%}

%%

[\t ]+ ;

\n {
  lineno++;
}

[a-zA-Z][a-zA-Z0-9\._\-]* {
  yylval.fname=new char[strlen(yytext)+1];
  strcpy(yylval.fname, yytext);
  return FNAME;
}

[0-9]+ {
  sscanf(yytext, "%d", &yylval.cls);
  return CLS;
}

\{ {
  return yytext[0];
}

\} {
  return yytext[0];
}

\| {
  return yytext[0];
}

%%

int yyerror(const char *s) {
  fprintf(stderr, "sparse_calc, line %d: %s at %s\n", lineno, s, yytext);
}

