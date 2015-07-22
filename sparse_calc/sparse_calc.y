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
  char *symbol;
  real scalar;
  long id;
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

%token <id> PLUS
%token <id> PROD
%token <id> CPROD 
%token <id> QUOTIENT
%token <id> MINUS
%token <id> RANGE
%token <id> POW
%token <id> MOD
%token <id> COMPARATOR

%token VECTOR_T

%token <symbol> SYMBOL
%token <symbol> LITERAL

%token <id> MATRIX_T
%token <id> SCALAR
%token <id> MATRIX
%token <id> VECTOR
%token <id> SUBROUTINE
%token <scalar> CONSTANT

%left COMPARATOR
%left PLUS MINUS
%left PROD CPROD QUOTIENT MOD
%left LEFT_SUB RANGE
%right POW
%left LEFT_BRAC
%nonassoc UMINUS

%type <id> list
%type <id> operator

%%

statement_list: statement | statement_list statement;

statement: 
	| command
	| assignment 
	| DELIM
	| error DELIM {
		  yyclearin;
		  yyerrok;
		};

list: 
  expression {
      $$=1;
    }
  | list COMMA expression {
      $$++;
    };

command:
  SUBROUTINE list DELIM {
      code_push_narg($2);
      code.push($1);
    }
  | VECTOR_T SYMBOL DELIM {
      long id;
      id=symtab.add($2);
      if (id<0) {
        yyerror("Error adding symbol\n");
        YYERROR;
      }
      vartyp[id]=SC_VECTOR_T;
      delete [] $2;
    }
  | MATRIX_T SYMBOL DELIM {
      long id;
      id=symtab.add($2);
      if (id<0) {
        yyerror("Error adding symbol\n");
        YYERROR;
      }
      vartyp[id]=$1;
      delete [] $2;
    };

operator:
  MINUS | PROD | PLUS | MOD | QUOTIENT | CPROD | RANGE | COMPARATOR | POW

expression: 
  SCALAR {
      code.push($1);
    }
  | VECTOR {
      code.push($1);
    }
  | MATRIX {
      code.push($1);
    }
  | CONSTANT {
      code_push_constant($1);
    }
  | PARAMETER {
      code.push($1);
    }
  | LITERAL {
      code_push_literal($1);
      delete [] $1;
    }
  | LEFT_BRAC expression RIGHT_BRAC {
      //do nothing
    }
  | MINUS expression %prec UMINUS {
      code_push_narg(1);
      code.push(sc_neg_code);
    }
  | expression operator expression {
      code_push_narg(2);
      code.push($2);
    }
  | expression LEFT_SUB expression RIGHT_SUB {
      code_push_narg(2);
      code.push(sc_sub_code);
    }
  | expression LEFT_SUB expression COMMA expression RIGHT_SUB {
      code_push_narg(3);
      code.push(sc_sub_code);
    }
  | NORMBRAC expression NORMBRAC {
      code_push_narg(1);
      code.push(sc_norm_code);
    }
  | VECTOR_T LEFT_BRAC SYMBOL RIGHT_BRAC {
      long id;
      id=symtab.add($3);
      if (id<0) {
        yyerror("Error adding symbol\n");
        YYERROR;
      }
      vartyp[id]=SC_VECTOR_T;
      code.push(id);
      delete [] $3;
    }
  | MATRIX_T LEFT_BRAC SYMBOL RIGHT_BRAC {
      long id;
      id=symtab.add($3);
      if (id<0) {
        yyerror("Error adding symbol\n");
        YYERROR;
      }
      vartyp[id]=$1;
      code.push(id);
      delete [] $3;
    }
  | SUBROUTINE LEFT_BRAC list RIGHT_BRAC {
      code_push_narg($3);
      code.push($1);
    };

assignment:
  SYMBOL LEFT_SUB expression RIGHT_SUB ASSIGN expression {
      long id;
      id=symtab.lookup($1);
      if (id<0) {
        yyerror("Variable undefined\n");
        YYERROR;
      }
      //first we subscript:
      code.push(id);
      code.push($3);
      code_push_narg(2);
      code.push(sc_sub_code);

      //element-by-element assignment:
      code.push($6);
      code_push_narg(2);
      code.push(SC_ELBYEL_CODE);
      delete [] $1;
    }
  | SYMBOL LEFT_SUB expression COMMA expression RIGHT_SUB ASSIGN expression {
      id=symtab.lookup($1);
      if (id<0) {
        yyerror("Variable undefined\n");
        YYERROR;
      }
      code.push(id);
      code.push($3);
      code.push($5);
      code_push_narg(3);
      code.push(sc_sub_code);

      code.push($8);
      code_push_narg(2);
      code.push(SC_ELBYEL_CODE);
      delete [] $1;
    }
  | SYMBOL ASSIGN expression {
      //"clobber" assignment:
      id=symtab.lookup($1);
      if (id<0) {
        id=symtab.add($1);
        vartyp[id]=SC_TO_BE_DEFINED;
      }
      code.push(id);
      code.push($3);
      code_push_narg(2);
      code.push(SC_CLOBBER_CODE);
      delete [] $1;
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

  int32_t maxniter=DEF_MAXNITER;
  int32_t nev=DEF_NEV;
  int32_t ncv=DEF_NARNOLDI;
  float tol=DEF_TOL;

  timeb date;
  pid_t pid;

  opt_arg[3]=&nev;
  opt_arg[4]=&ncv;
  opt_arg[5]=&tol;
  opt_arg[6]=&maxniter;

  argc=parse_command_opts(argc, argv, "dhpAvtILe", "%%%s%d%d%d%f%s%s", opt_arg, flag, 1);

  if (argc<0) {
    fprintf(stderr, "sparse_calc: Error parsing command line.\n");
    argc=-argc;
  }

  scalar_assign("MAXNITER", maxniter);
  scalar_assign("NEV", nev);
  scalar_assign("NARNOLDI", ncv);
  scalar_assign("TOL", tol);

  //user can make variable declarations at the command line:
  //(repeat options are not useful with parse_command_opts...)
  while ((c = getopt(argc, argv, "S:M:a:V")) != -1) {
    switch (c) {
      case ('S'):
        id=symtab.add(optarg);
        vartyp[id]=SPARSE_T;
        delflag.off(id);
        break;
      case ('M'):
        id=symtab.add(optarg);
        vartyp[id]=FULL_T;
        delflag.off(id);
        break;
      case ('a'):
        id=symtab.add(optarg);
        vartyp[id]=SPARSE_ARRAY_T;
        delflag.off(id);
        break;
      case ('V'):
        id=symtab.add(optarg);
        vartyp[id]=VEC_T;
        delflag.off(id);
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
    main_help_screen();
    return 0;
  }

  if (flag[7]) {
    read_history((char *) opt_arg[7]);
  }
  if (flag[8]) {
    if (run_script((char *) opt_arg[8])!=0) {
      fprintf(stderr, "sparse_calc: error running script, %s\n", (char *) opt_arg[8]);
      fprintf(stderr, "   ...starting interactive mode.\n");
    }
  }

  ftime(&date);
  pid=getpid();

  scalc_session_id=pid+date.time*100000;	//how big do PID's get??

  if (argc>0) {
    if (run_script(argv[0])!=0) {
      fprintf(stderr, "sparse_calc: error running script, %s\n", argv[0]);
      exit(UNABLE_TO_OPEN_FILE_FOR_READING);
    }
    if (flag[9]) {
      sc_makeflag=1;
    }
    sc_scriptflag=1;
  }

  yyparse();

  //clean up:
  printf("Deleting unsaved variables:\n");
  for (long i=0; i<symtab.entries(); i++) {
    char *symbol;
    char *command;
    if (vartyp[i]!=SCALAR_DATA) {
      if (delflag[i]) {
        symbol=symtab.get(i);
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

