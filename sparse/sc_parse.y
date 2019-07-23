%{
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <sys/stat.h>
#include <readline/history.h>

extern "C" {
  int yylex(void);
  int yyerror(const char *);
}

#include "parse_command_opts.h"
#include "sparse_calc_defs.h"
#include "error_codes.h"

#define YYDEBUG 1

using namespace sparse_calc;

%}

%union {
  sc_type_base * var;
  sc_fun_t sc_fun;
}

%token ERROR
%token DELIM

%token LEFT_BRAC
%token RIGHT_BRAC

%token LEFT_SUB
%token RIGHT_SUB
%token COMMA
%token NORMBRAC
%token ASSIGN

//infix operators:

%token PLUS
%token PROD
%token CPROD 
%token QUOTIENT
%token MINUS
%token RANGE
%token POW
%token MOD
%token GT
%token LT
%token GE
%token LE
%token EQ

%token <var> SYMBOL
%token <var> LITERAL

%token <var> VARIABLE
%token <sc_fun> SUBROUTINE

%left GT LT GE LE EQ
%left PLUS MINUS
%left PROD CPROD QUOTIENT MOD
%left LEFT_SUB RANGE
%right POW
%left LEFT_BRAC
%nonassoc UMINUS

%type <var> expression
%type <var> list
%type <var> list_expression
%type <var> call

%%

statement_list: statement | statement_list statement;

statement: 
  | assignment
  | call DELIM {
      delete $1;
    }
  | DELIM
  | error DELIM {
    yyclearin;
    yyerrok;
  };

assignment:
  SYMBOL LEFT_SUB expression RIGHT_SUB ASSIGN expression {
      sc_type_base * lval=sc_var_lookup(&sc_state, $1);
      if (lval==NULL) {
        yyerror("Variable undefined\n");
        YYERROR;
        //(looks like a memory leak...)
      }
      lval->subscript_assign($3, $6);
      delete $3;
      delete $6;
      sc_var_assign($1, lval);
      delete lval;
      delete $1;
    }
  | SYMBOL LEFT_SUB expression COMMA expression RIGHT_SUB ASSIGN expression {
      sc_type_base * lval=sc_var_lookup($1);
      if (lval==NULL) {
        yyerror("Variable undefined\n");
        YYERROR;
      }
      lval->subscript_assign($3, $5, $8);
      delete $3;
      delete $5;
      delete $8;
      sc_var_assign($1, lval);
      delete lval;
      delete $1;
    }
  | SYMBOL ASSIGN expression {
      sc_var_assign(&sc_state, $1, $3);
    };

call:
  SYMBOL list_expression {
    $$=sc_call_user_fun(&sc_state, $1, $2);
    delete $1;
    delete $2;
    }
  | SUBROUTINE list_expression {
    $$=(*$1)(&sc_state, $2);
    delete $1;
    delete $2;
  };

list_expression:
  LEFT_BRAC RIGHT_BRAC {
      $$=new sc_type_list();
    }
  | LEFT_BRAC list RIGHT_BRAC {
      $$=$2;
  };

list:
  expression {
      $$=new sc_type_list();
      $$->multiply($1);
    }
  | list COMMA expression {
      $1->multiply($3);
      $$=$1;
    };

expression:
  VARIABLE {
      $$=$1;
    }
  | list_expression {
      $$=$1;
    }
  | LEFT_BRAC expression RIGHT_BRAC {
      $$=$2;
    }
  | MINUS expression %prec UMINUS {
      $$=$2->negate();
      if ($$==NULL) {
        yyerror("Syntax error\n");
        YYERROR;
      }
    }
  | NORMBRAC expression NORMBRAC {
      $$=$2->norm();
      if ($$==NULL) {
        yyerror("Syntax error\n");
        YYERROR;
      }
    }
  | call {
      $$=$1;
    }
  | expression MINUS expression {
      $3->negate();
      $$=$1->add($3);
      delete $1;
      delete $3;
      if ($$==NULL) {
        yyerror("Syntax error\n");
        YYERROR;
      }
    }
  | expression PLUS expression {
      $$=$1->add($3);
      delete $1;
      delete $3;
      if ($$==NULL) {
        yyerror("Syntax error\n");
        YYERROR;
      }
    }
  | expression PROD expression {
      $$=$1->multiply($3);
      delete $1;
      delete $3;
      if ($$==NULL) {
        yyerror("Syntax error\n");
        YYERROR;
      }
    }
  | expression QUOTIENT expression {
      $$=$1->divide($3);
      delete $1;
      delete $3;
      if ($$==NULL) {
        yyerror("Syntax error\n");
        YYERROR;
      }
    }
  | expression POW expression {
      $$=$1->pow($3);
      delete $1;
      delete $3;
      if ($$==NULL) {
        yyerror("Syntax error\n");
        YYERROR;
      }
    }
  | expression GT expression {
      $$=$1->gt($3);
      delete $1;
      delete $3;
      if ($$==NULL) {
        yyerror("Syntax error\n");
        YYERROR;
      }
    }
  | expression LT expression {
      $$=$1->lt($3);
      delete $1;
      delete $3;
      if ($$==NULL) {
        yyerror("Syntax error\n");
        YYERROR;
      }
    }
  | expression GE expression {
      $$=$1->ge($3);
      delete $1;
      delete $3;
      if ($$==NULL) {
        yyerror("Syntax error\n");
        YYERROR;
      }
    }
  | expression LE expression {
      $$=$1->le($3);
      delete $1;
      delete $3;
      if ($$==NULL) {
        yyerror("Syntax error\n");
        YYERROR;
      }
    }
  | expression EQ expression {
      $$=$1->eq($3);
      delete $1;
      delete $3;
      if ($$==NULL) {
        yyerror("Syntax error\n");
        YYERROR;
      }
    }
  | expression NE expression {
      $$=$1->eq($3);
      delete $1;
      delete $3;
      if ($$==NULL) {
        yyerror("Syntax error\n");
        YYERROR;
      }
    }
  | expression MOD expression {
      $$=$1->mod($3);
      delete $1;
      delete $3;
      if ($$==NULL) {
        yyerror("Syntax error\n");
        YYERROR;
      }
    }
  | expression RANGE expression {
      $$=new sc_type_vector($1, $3);
      delete $1;
      delete $3;
      if ($$==NULL) {
        yyerror("Syntax error\n");
        YYERROR;
      }
    }
  | expression CPROD expression {
      $$=$1->multiply_all($3);
      delete $1;
      delete $3;
      if ($$ == NULL) {
        yyerror("Syntax error\n");
        YYERROR;
      }
    }
  | expression LEFT_SUB expression RIGHT_SUB {
      $$=$1->subscript($3);
      delete $3;
      delete $1;
      if ($$==NULL) {
        yyerror("Syntax error\n");
        YYERROR;
      }
    }
  | expression LEFT_SUB expression COMMA expression RIGHT_SUB {
      $$=$1->subscript($3, $5);
      delete $1;
      delete $3;
      delete $5;
      if ($$==NULL) {
        yyerror("Syntax error\n");
        YYERROR;
      }
    }
  | expression LEFT_BRAC expression RIGHT_BRAC {
      $$=$1->get_element($3)
      delete $3;
      delete $1;
      if ($$==NULL) {
        yyerror("Syntax error\n");
        YYERROR;
      }
  };


%%

#include <sys/timeb.h>
#include <sys/types.h>
#include <unistd.h>
#include <getopt.h>

int main (int argc, char **argv) {
  FILE *fs;
  string_petey dum;
  char *line;
  char c;
  long id;

  int hflag;
  int flag[20];
  void *opt_arg[20];
  int loc;

  timeb date;
  pid_t pid;

  argc=parse_command_opts(argc, argv, "dhpb", "%%%s%", opt_arg, flag, 1);

  if (argc<0) {
    fprintf(stderr, "sparse_calc: Error parsing command line.\n");
    argc=-argc;
  }

  //user can make variable declarations at the command line:
  //(repeat options are not useful with parse_command_opts...)
  while ((c = getopt(argc, argv, "s:m:a:v")) != -1) {
    switch (c) {
      case ('s'):
        sc_add_arg(optarg, SPARSE);
        break;
      case ('m'):
        sc_add_arg(optarg, FULL);
        break;
      case ('a'):
        sc_add_arg(optarg, SPARSE_ARRAY);
        break;
      case ('v'):
        sc_add_arg(optarg, VECTOR);
        break;
      case ('?'):
        fprintf(stderr, "sparse_calc: Unknown option -%c -- ignored\n", optopt);
        break;
      default:
        fprintf(stderr, "sparse_calc: Error parsing command line.\n");
        break;
    }
        
  }

  argc-=optind;
  argv+=optind;

  //symtab.print();

  yydebug=flag[0];

  if (flag[2]) {
    path=(char *) opt_arg[2];
    loc=strlen(path);
    if (path[loc-1]!='/') {
      path=new char[loc+2];
      sprintf(path, "%s/", (char *) opt_arg[2]);
      delete [] (char *) opt_arg[2];
    }
  } else {
    path=new char[3];
    strcpy(path, "./");
  }

  if (flag[1]) {
    sc_help(NULL);
    return 0;
  }

  ftime(&date);
  pid=getpid();

  sc_session_id=pid+date.time*100000;		//how big do PID's get??

  if (argc>0) {
    sc_ifstream=fopen(argv[0], "r");
  } else {
    sc_ifstream=stdin;
  }

  yyparse();

  //clean up:
  printf("Deleting unsaved variables:\n");
  for (long i=0; i<sc_vartab.entries(); i++) {
    char *symbol;
    char *command;
    if (sc_vlist[i]!=NULL) delete sc_vlist[i];
    if (sc_vartype[i]!=SCALAR && sc_vartype[i]!=LITERAL) {
      if (delflag[i]) {
        symbol=sc_vartab.get(i);
        command=new char[strlen(symbol)+10];
        //delete variables that haven't been saved
        //(doesn't help you if the process dies or is killed in the middle of a session...)
        sprintf(command, "rm %s", symbol);
        printf("%s\n", command);
        system(command);
        delete [] symbol;
        delete [] command;
      }
    }
  }

}

