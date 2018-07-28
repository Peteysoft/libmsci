%{
extern "C" {
  int yylex(void);
  int yyerror(const char *);
}

#include "parse_command_opts.h"

#include "datasets_lang.h"

#define YYDEBUG 1


%}

%union {
  double scalar;
  char *symbol;
  char *op;
  dataset *ds;
}

%token STMT_DELIM
%token CD_DELIM

%token <scalar> SCALAR
%token <symbol> SYMBOL
%token <op> INFIX_OP

%type <full_name> symbol;
%type <subexpr> ds;
%type <expression> ds;

%left '-' '+'
%left '*' '/'


%%

statement_list: statement | statement_list statement;

statement: declaration 
  | assignment 
  | error STMT_DELIM {
    yyclearin;
    yyerrok;
  };

full_name: SYMBOL {
    $$=$1;
    if (qualflag==0) second_qualifier=first_qualifier;
    qualflag=1;
  } | full_name CD_DELIM SYMBOL {
    long loc, dum;
    dataset *checkds;
    loc=first_qualifier->search_var($1, dum);
    if (loc < 0) {
      yyerror("Symbol undefined within current scope\n");
      YYERROR;
    }
    //descend down the hierarchy:
    checkds=first_qualifier->get_var(loc);
    //must be a composite dataset:
    if (checkds->typeof() != COMPOSITE) {
      yyerror("Improper type within fully-qualified name\n");
      YYERROR;
    }
    first_qualifier=checkds;
    $$=$2;
    if (qualflag==0) second_qualifier=first_qualifier;
    qualflag=1;
  };

dependent_list: SYMBOL {
    long loc, dum;
    dataset *ds;
    loc=first_qualifier->search_var($1, dum);
    if (loc < 0) {
      yyerror("Symbol undefined within current scope\n");
      YYERROR;
    }
    ds=first_qualifier->get_var(loc);
    //must be a simple dataset:
    if (ds->typeof() >= DEPENDENT || ds->typeof() < SIMPLE_BASE) {
      yyerror("Improper type wthin dependent declaration\n");
      YYERROR;
    }
    argkey[narg]=(simple_dataset *) ds;
    narg++;
  }| dependent_list ',' SYMBOL {
    long loc, dum;
    dataset *ds;
    loc=first_qualifier->search_var($1, dum);
    if (loc < 0) {
      yyerror("Symbol undefined within current scope\n");
      YYERROR;
    }
    ds=first_qualifier->get_var(loc);
    //must be a simple dataset:
    if (ds->typeof() >= DEPENDENT || ds->typeof() < SIMPLE_BASE) {
      yyerror("Improper type wthin dependent declaration\n");
      YYERROR;
    }
    argkey[narg]=(simple_dataset *) ds;
    narg++;
  };

declaration:
  SYMBOL full_name STMT_DELIM {
    long loc;
    dataset *ds;
    loc=first_qualifier->add_var($2);
    ds=allocate_simple($1);
    if (ds==NULL) {
      yyerror("Type name not recognized\n");
      YYERROR;
    }
    first_qualifier->cvar(loc, ds);
  } | SYMBOL full_name '(' dependent_list ')' STMT_DELIM {
    long loc;
    dataset *ds;
    loc=first_qualifier->add_var($2);
    ds=allocate_dependent($1);
    if (ds==NULL) {
      yyerror("Type name not recognized\n");
      YYERROR;
    }
    first_qualifier->cvar(loc, ds);
  };

argument: SYMBOL '=' expression {
    long loc, dum;
    dataset *ds;
    //first we look up the symbol within the current scope:
    loc=second_qualifier->search_var($1, dum);
    if (loc < 0) {
      yyerror("Symbol undefined within current scope\n");
      YYERROR;
    }
    ds=second_qualifier->get_var(loc);
    //must be a simple dataset:
    if (ds->typeof() >= DEPENDENT || ds->typeof() < SIMPLE_BASE) {
      yyerror("Improper type wthin dependent declaration\n");
      YYERROR;
    }
    argkey[narg]=(simple_dataset *) ds;
    arglist[narg]=$2;
    narg++;
  };

argument_list: '(' argument ')' | '(' argument_list ',' argument ')';

subexpr: full_name argument_list {
    long loc, dum;
    dataset *ds;
    loc=first_qualifier->search_var($1, dum);
    if (loc < 0) {
      yyerror("Symbol undefined within current scope\n");
      YYERROR;
    }
    ds=first_qualifier->get_var(loc);
    //must be a dependent dataset:
    if (ds->typeof() >= COMPOSITE || ds->typeof() < DEPENDENT) {
      yyerror("Subexpression has improper type\n");
      YYERROR;
    }
    $$=dataset_interpolate(ds);
  };

expression: SYMBOL {
    long loc, dum;
    dataset *ds;
    loc=second_qualifier->search_var($1, dum);
    if (loc < 0) {
      yyerror("Symbol undefined within current scope\n");
      YYERROR;
    }
    ds=second_qualifier->get_var(loc);
    //must be a dependent or simple dataset:
    if (ds->typeof() >= COMPOSITE || ds->typeof() < SIMPLE_DATASET) {
      yyerror("Type illegal in expression\n");
      YYERROR;
    }
    $$=ds;
  } | subexpr {
    if (push_heap($1)!=0) YYERROR;
    $$=$1;
  } | expression INFIX_OP expression
  | SCALAR;

assignment:
  full_name '=' expression STMT_DELIM
  | full_name '=' subexpr STMT_DELIM
  | subexpr '=' expression STMT_DELIM;

%%

