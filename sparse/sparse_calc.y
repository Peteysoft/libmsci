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
  real scalar;
  char *symbol;
  vector_t vector;
  matrix_t *matrix;
  int type_id;
  int funcid;
  char *literal;
}

%token ERROR
%token ENDFIRSTPASS
%token LEFT_BRAC RIGHT_BRAC LEFT_SUB RIGHT_SUB COMMA PLUS PROD CPROD 
		ASSIGN QUOTIENT MINUS DELIM NORMBRAC RANGE POW MOD
%token PRINT PATH HELP SYSTEM COMMENT DELETE RUN
%token SQRT SIN COS TAN LOG
%token SCALAR_T VECTOR_T
%token IDENTITY TRANSPOSE EVD SVD SOLVE INVERT EMPTY READ SIZE SUM SAVE
%token ALL

%token <scalar> SCALAR
%token <symbol> SYMBOL
%token <vector> VEC
%token <matrix> MATRIX
%token <type_id> MATRIX_T
%token <funcid> FUNCTION
%token <literal> LITERAL

%left PLUS MINUS
%left PROD CPROD QUOTIENT MOD
%left LEFT_SUB RANGE
%right POW
%left LEFT_BRAC
%nonassoc UMINUS

%type <vector> vector_exp
%type <matrix> matrix_exp
%type <scalar> scalar_exp

%%

statement_list: statement | statement_list statement;

statement: system 
	| comment
	| help 
	| printout 
	| include
	| addpath 
	| declaration 
        | script
	| assignment 
	| DELIM
        | ENDFIRSTPASS {
		  yyclearin;
                  scalc_istream=fmemopen(sc_stringbuf, sc_buflen, "r");
                  sc_firstpass=2;
                }
	| error DELIM {
		  yyclearin;
		  yyerrok;
		};

comment:
  COMMENT DELIM {}
  | COMMENT LITERAL DELIM {
      delete [] $2;
    };

system:
  SYSTEM LITERAL DELIM {
      int err=system($2);
      if (err!=0) {
        fprintf(stderr, "sh returned %d exit code\n", err);
        yyerror("External error");
      }
      delete [] $2;
    };

help:
  HELP DELIM {
      main_help_screen();
    }
  | HELP HELP DELIM {
      printf("\n$ help <expression|operator|function|type|literal>\n");
      printf("  where:\n");
      printf("  - <expression> is a scalar, vector or matrix\n");
      printf("  - <operator> is one of: ()[]|*/+-#=:\n");
      printf("  - <function> is a function or command name:\n");
      printf("      empty identity read size\n");
      printf("      full sparse sparse_array\n");
      printf("      transpose solve invert evd svd\n");
      printf("      path print help sys com run\n");
      printf("  - <type> is a type name:\n");
      printf("      full sparse sparse_array vector scalar\n");
      printf("  - <literal> is a \"quoted\" string:\n");
      printf("      \"types\" \"vars\" \"funcs\" \"commands\" \"operators\" \"reserved\"\n");
      printf("      \"cg\" \"bcg\" \"evd\" \"svd\"\n");
      printf("      or a variable name\n");
      printf("\n");
    }
  | HELP MATRIX_T DELIM {
      type_help ($2);
    }
  | HELP VECTOR_T DELIM {
      printf("Vector type.  Use in declaration as follows:\n");
      printf("$ vector <variable>\n");
      printf("or as a \'shortcut\' declaration:\n");
      printf("  <vector>=vector(<name>)\n");
      printf("  - which simultaneously declares the variable and returns its contents\n");
      printf("- Vector is read in from a binary file of the same name.\n");
      printf("- Declared variables persist after exit\n");
      break;
    }
  | HELP scalar_exp DELIM {
      printf("Scalar = %g\n", $2);
    }
  | HELP vector_exp DELIM {
      printf("Vector of length, %d\n", $2.n);
      delete [] $2.data;
    }
  | HELP matrix_exp DELIM {
      integer type=matrix_type($2);
      integer m, n, nel;
      $2->dimensions(m, n);
      switch (type) {
        case (FULL_T):
          printf("[%dx%d] full matrix\n", m, n);
          break;
        case (SPARSE_T):
          nel=((sparse_t *) $2)->size();
          printf("[%dx%d] sparse matrix; %d non-zero elements\n", m, n, nel);
          break;
        case (SPARSE_ARRAY_T):
          nel=((sparse_array_t *) $2)->nel();
          printf("[%dx%d] array of sparse matrices; %d members\n", m, n, nel);
          break;
      }
      delete $2;
    }
  | HELP EMPTY DELIM {
      printf("\n<expression>=empty(<m> [,<n>])\n");
      printf("  Returns an empty (zero) vector or matrix of the specified dimensions\n\n");
    }
  | HELP IDENTITY DELIM {
      printf("\n<expression>=identity(<m> [,<n>])\n");
      printf("  Returns an [<m>x<n>] identity matrix\n\n");
    }
  | HELP RANGE DELIM {
      printf("\n<vector>=<i1>:<i2>\n");
      printf("  Returns a consecutive arithmetic sequence\n\n");
    }
  | HELP TRANSPOSE DELIM {
      printf("\n<transpose>=transpose(<matrix>)\n");
      printf("  Returns the transpose of a matrix.\n\n");
    }
  | HELP SOLVE DELIM {
      printf("\n<solution-vector>=solve(<matrix>, <vector>, \"<method\")\n");
      printf("  Solves the matrix equation using <method>\n");
      printf("  Supported solution methods:\n");
      printf("    cg	conjugate gradient\n");
      printf("    bcg	bi-conjugate gradient\n");
      printf("    evd	eigenvalue decomposition\n");
      printf("    svd	singular value decomposition\n");
      printf("    n:<method>	solves the normal equation using <method>\n\n");
    }
  | HELP INVERT DELIM {
      printf("\n<inverse-matrix>=invert(<matrix>, \"<method\")\n");
      printf("  Inverts <matrix> using <method>\n");
      printf("  Supported solution methods:\n");
      printf("    cg	conjugate gradient\n");
      printf("    bcg	bi-conjugate gradient\n");
      printf("    evd	eigenvalue decomposition\n");
      printf("    svd	singular value decomposition\n");
      printf("Also:\n");
      printf("<inverse-vector>=invert(<vector>)\n\n");
    }
  | HELP ASSIGN DELIM {
      printf("\nAssignment operator:\n");
      printf("$ <lvalue>=<expression>\n");
      printf("  - Vectors and matrices are written to a disc in a file of the same name.\n");
      printf("  - Scalars are stored in RAM\n");
      printf("  - Any new assignment clobbers the old variable, regardless of a previous declaration.\n");
      printf("  - For type conversion, see: full sparse sparse_array\n");
      printf("  - L-values may take the following forms:\n");
      printf("    $ <vector>[<scalar>]=<scalar>\n");
      printf("    $ <vector>[<vector>]=<vector>\n");
      printf("    $ <matrix>[<scalar>, <scalar>]=<scalar>\n");
      printf("    $ <matrix>[<vector>, <vector>]=<vector>\n");
      printf("  - Subscripts and right-hand-side must agree in type and size\n\n");
      printf("  - Assigned variables are deleted on exit.\n\n");
    }
  | HELP PROD DELIM {
      printf("\nProduct operator.  Multiplation between all types is supported:\n");
      printf("- vector times a matrix returns a vector\n");
      printf("- full matrix times another matrix or a scalar returns a full matrix\n");
      printf("- sparse matrix times sparse matrix\n");
      printf("  - returns a sparse array if both are square and the dimensions the same,\n");
      printf("  - another sparse matrix otherwise\n");
      printf("- sparse array times sparse matrix or sparse array\n");
      printf("  - returns a sparse array if the dimensions are the same,\n");
      printf("  - full matrix otherwise\n");
      printf("  - vector times vector returns a scalar\n\n");
    }
  | HELP PLUS DELIM {
      printf("\nAddition operator.  Addition between all types is supported:\n");
      printf("- full matrix or sparse matrix plus another matrix returns a full matrix\n");
      printf("- sparse matrix plus sparse matrix returns a sparse matrix\n");
      printf("- addition of a scalar to a matrix or vector adds the scalar to all elements\n");
      printf("- a vector plus a matrix adds the elements of the vector along the diagonal\n\n");
    }
  | HELP MINUS DELIM {
      printf("\n- Minus operator negates and expression or adds the negated expression.\n\n");
    }
  | HELP QUOTIENT DELIM {
      printf("\n- Division currently only supported for scalars.\n");
      printf("- Use solve or invert for vectors and matrices.\n\n");
    }
  | HELP NORMBRAC DELIM {
      printf("\n<scalar>=|<expression>|\n");
      printf("- Returns the norm of <expression>\n\n");
    }
  | HELP CPROD DELIM {
      printf("\nCumulative product operator:\n");
      printf("  <full-matrix>=<sparse-array>#<vector>\n");
      printf("  <full-matrix>=<vector>#<sparse-array>\n\n");
    }
  | HELP LEFT_BRAC DELIM {
      printf("\nBrackets are used:\n");
      printf("  - to over-ride the normal operator precedence.\n");
      printf("  - to enclose function arguments.\n");
      printf("  - to extract elements from a sparse matrix array\n\n");
      printf("Operator precedence is as follows:\n");
      printf("  1. []:\n");
      printf("  2. */#\n");
      printf("  3. +-\n\n");
    }
  | HELP RIGHT_BRAC DELIM {
      printf("\nBrackets are used:\n");
      printf("  - to over-ride the normal operator precedence.\n");
      printf("  - to enclose function arguments.\n");
      printf("  - to extract elements from a sparse matrix array\n\n");
      printf("Operator precedence is as follows:\n");
      printf("  1. []:\n");
      printf("  2. */#\n");
      printf("  3. +-\n\n");
    }
  | HELP LEFT_SUB DELIM {
      printf("\nSubscript operator:\n");
      printf("  <expression>=<matrix>[<expression>]\n");
      printf("  <expression>=<matrix>[<expression>, <expression>]\n");
      printf("  <expression>=<vector>[<expression>]\n");
      printf("  - Both scalars and vectors may be used as subscripts.\n");
      printf("  - Subscript operators may also be used in l-values\n\n");
    }
  | HELP RIGHT_SUB DELIM {
      printf("\nSubscript operator:\n");
      printf("  <expression>=<matrix>[<expression>]\n");
      printf("  <expression>=<matrix>[<expression>, <expression>]\n");
      printf("  <expression>=<vector>[<expression>]\n");
      printf("  - Both scalars and vectors may be used as subscripts.\n");
      printf("  - Subscript operators may also be used in l-values\n\n");
    }
  | HELP EVD DELIM {
     printf("\n<eigenvalue-matrix>=evd(<matrix> [,<nev> [, \"<which>\"]])\n");
     printf("- Calculates the eigenvalues and eigenvectors of a matrix\n");
     printf("  - Returns <nev> eigenvectors in a full matrix.  \n");
     printf("  - Eigenvalues are stored in the [2 x <nev>] matrix EVAL.\n");
     printf("  - Uses the ARPACK routine snaupd\n");
     printf("  - If <nev> is not supplied, takes it from NEV scalar variable, 5 is default\n");
     printf("  - NARNOLDI supplies the number of Arnoldi vectors, 20 is the default\n");
     printf("  - MAXNITER supplies the maximum number of iterations\n");
     printf("  - TOL supplies the tolerance\n");
     printf("  - <which> is which eigenvalues to calculate:\n");
     printf("    - LM (default) largest in magnitude\n");
     printf("    - SM smallest in magnitude\n");
     printf("    - LA largest algebraicly\n");
     printf("    - SA smallest algebraicly\n");
     printf("    - BE takes half from either side of the spectrum\n\n");
    }
  | HELP SVD DELIM {
     printf("\n<eigenvalue-matrix>=svd(<matrix> [,<nev> [,\"<which>\"]])\n");
     printf("- Calculates the eigenvalue and eigenvectors of a symmetric matrix\n");
     printf("  - Returns <nev> eigenvectors in a full matrix.  \n");
     printf("  - Eigenvalues are stored in the vector EVAL.\n");
     printf("  - Uses the ARPACK routine ssaupd\n");
     printf("  - If <nev> is not supplied, takes it from NEV scalar variable, 5 is default\n");
     printf("  - NARNOLDI supplies the number of Arnoldi vectors, 20 is the default\n\n");
     printf("  - MAXNITER supplies the maximum number of iterations\n");
     printf("  - TOL supplies the tolerance\n");
     printf("  - <which> is which eigenvalues to calculate:\n");
     printf("    - LM (default) largest in magnitude\n");
     printf("    - SM smallest in magnitude\n");
     printf("    - LA largest algebraicly\n");
     printf("    - SA smallest algebraicly\n");
     printf("    - BE takes half from either side of the spectrum\n\n");
    }
  | HELP PRINT DELIM {
      printf("\n$ print [\"<filename>\"] <expression>\n");
      printf("  - Prints out the contents of <expression> to standard out\n");
      printf("  - or to <filename>\n");
      printf("\n");
      printf("$ print [\"<filename>\"] <literal>\n");
      printf("  - Where <literal> is \"vars\" or \"history\"\n");
      printf("  - run may be used to re-run a saved command history\n");
      printf("    or re-define listed variables from a previous session\n");
      printf("  - A saved history may be restored using the -H command-line option\n");
      printf("\n");
    }
  | HELP SAVE DELIM {
      printf("\n$ save \"<variable>\"\n");
      printf("  - Save <variable> (don't delete on exit).\n");
    }
  | HELP RUN DELIM {
      printf("\n$ run \"<filename>\"\n");
      printf("  - Runs the commands contained in <filename> just as if they had been\n");
      printf("    typed in the command line\n");
      printf("  - A session may be saved (to be run later) using the print command\n\n");
    }
  | HELP PATH DELIM {
      printf("\n$ path \"<path>\"\n");
      printf("  - Changes the storage directory to <path>\n");
      printf("  - Vector and matrix variable are deleted and stored in a file called:\n");
      printf("    sc.vars.<id>.sc -- where <id> is a unique session id\n");
      printf("  - The current data path is:\n");
      printf("%s\n\n", path);
    }
  | HELP READ DELIM {
      printf("\n<expression>=read (<type>, \"<filename>\"\n");
      printf("  - Reads ASCII data of type <type> from <filename>.\n");
      printf("  - Inverse of print when print is used to print out expressions\n\n");
    }
  | HELP COMMENT DELIM {
      printf("\n$ com [\"<comment>\"]\n");
      printf("  Comment string\n\n");
    }
  | HELP SYSTEM DELIM {
      printf("\n$ sys \"<command>\"\n");
      printf("  Executes a command under the shell (sh)\n\n");
    }
  | HELP DELETE DELIM {
      printf("\n$ delete \"<variable>\"\n");
      printf("  - Deletes <variable> so that it may be redefined.\n");
      printf("  - Data file is deleted on exit unless variable redefined.\n\n");
    }
  | HELP SYMBOL DELIM {
      printf("\n%s is neither a defined variable nor a reserved word\n\n", $2);
      delete [] $2;
    }
  | HELP SIZE DELIM {
      printf("\n<vector>=size(<expression>)\n");
      printf("  Returns a vector containing size and type information:\n");
      printf("  0. total number of elements\n");
      printf("  1. number of rows\n");
      printf("  2. number of columns\n");
      printf("  3. type (%d=sparse, %d=sparse_array, %d=full, %d=vector, %d=scalar\n",
		SPARSE_T, SPARSE_ARRAY_T, FULL_T, VEC_T, SCALAR_DATA);
      printf("  4. for sparse matrices, number of nonzero elements,\n");
      printf("     for sparse arrays, number of matrices\n");
      printf("     number of elements for all others\n\n");
    }
  | HELP LITERAL DELIM {
     literal_help($2);
     delete [] $2;
   };

include:
  RUN LITERAL DELIM {
      if (run_script($2)!=0) {
        yyerror("I/O error");
      }
      delete [] $2;
    };

printout:
  PRINT LITERAL DELIM {
      if (strcmp($2, "vars")==0) {
        print_allvar(stdout);
      } else if (strcmp($2, "saved")==0) {
        print_allvar(stdout, 2);
      } else if (strcmp($2, "history")==0) {
        write_history("/dev/stdout");
      } else {
        fprintf(stderr, "print: \"%s\" not a recognized option\n", $2);
        yyerror("Syntax error");
      }
      delete [] $2;
    }
  | PRINT LITERAL LITERAL DELIM {
      FILE *fs;
      if (strcmp($3, "history")==0) {
        if ($2[0]=='+') {
          if (write_history($2+1)) yyerror("I/O error\n");
        } else {
          fs=fopen($2, "w");
          if (fs==NULL) {
            yyerror("I/O error\n");
          } else {
            fclose(fs);
            if (write_history($2)) yyerror("I/O error\n");
          }
        }
      } else {
        if ($2[0]=='+') fs=fopen($2+1, "a"); else fs=fopen($2, "w");
        if (fs==NULL) {
          yyerror("I/O error\n");
        } else {
          if (strcmp($3, "vars")==0) {
            print_allvar(fs);
          } else if (strcmp($3, "saved")==0) {
            print_allvar(fs, 2);
          } else {
            fprintf(stderr, "print: \"%s\" not a recognized option\n", $3);
            yyerror("Syntax error");
          }
          fclose(fs);
        }
      }
      delete [] $2;
      delete [] $3;
    }
  | PRINT LITERAL matrix_exp DELIM {
      FILE *fs;
      if ($2[0]=='+') fs=fopen($2+1, "a"); else fs=fopen($2, "w");
      if (fs==NULL) {
        yyerror("I/O error\n");
      } else {
        $3->print(fs);
        fclose(fs);
      }
      delete [] $2;
      delete $3;
    }
  | PRINT LITERAL vector_exp DELIM {
      FILE *fs;
      if ($2[0]=='+') fs=fopen($2+1, "a"); else fs=fopen($2, "w");
      if (fs==NULL) {
        yyerror("I/O error\n");
      } else {
        for (integer i=0; i<$3.n; i++) fprintf(fs, "%g\n", $3.data[i]);
        fclose(fs);
      }
      delete [] $2;
      delete [] $3.data;
    }
  | PRINT LITERAL scalar_exp DELIM {
      FILE *fs;
      if ($2[0]=='+') fs=fopen($2+1, "a"); else fs=fopen($2, "w");
      if (fs==NULL) {
        yyerror("I/O error\n");
      } else {
        fprintf(fs, "%g\n", $3);
        fclose(fs);
      }
      delete [] $2;
    }
  | PRINT matrix_exp DELIM {
      $2->print(stdout);
      delete $2;
    }
  | PRINT vector_exp DELIM {
      for (integer i=0; i<$2.n; i++) printf("%g\n", $2.data[i]);
      delete [] $2.data;
    }
  | PRINT scalar_exp DELIM {
      printf("%g\n", $2);
    };

addpath:
  PATH DELIM {
      printf("The current path is:\n%s\n", path);
    }
  | PATH LITERAL DELIM {
      FILE *fs;
      int len;
      struct stat buf;
      int err;
      char savefile0[40];
      char *savefile;

      sprintf(savefile0, "sc.vars.%ld.sc", scalc_session_id);

      err=stat($2, &buf);
      if (err==0 && S_ISDIR(buf.st_mode)) { 
        //save all the variable in the old directory:
        savefile=add_path(savefile0);
        fs=fopen(savefile, "w");
        if (fs!=NULL) {
          print_allvar(fs, 1);
          for (integer i=0; i<symtab.entries(1); i++) {
            char *sym;
            //if (vartyp[i]!=SCALAR_DATA) vartyp[i]=DELETED;
            vartyp[i]=DELETED;
            sym=symtab.get(i);
            //if (delflag[(long) i]) fprintf(fs, "delete %s \"m\"\n", sym);
            if (delflag[i]) fprintf(fs, "delete %s \"m\"\n", sym);
            delete [] sym;
          }
          fclose(fs);
        }

        //set the new directory:
        delete [] path;
        path=$2;
        len=strlen(path);
        if (path[len-1]!='/') {
          path=new char[len+2];
          sprintf(path, "%s/", $2);
          delete [] $2;
        }

        //check for a savefile in the new directory:
        savefile=add_path(savefile0);
        run_script(savefile);
      } else {
        fprintf(stderr, "The path at, %s, is not accessible\n", $2);
        yyerror("Path not accessible");
      }
      delete [] $2;
    };

declaration:
  MATRIX_T SYMBOL DELIM {
      int id;
      //printf("declaration of type %d\n", $1);

      id=symtab.lookup($2);
      if (id<0) {
        id=symtab.add($2);
      }
      vartyp[id]=$1;
      //delflag.off(id);
      delflag[id]=0;
      //printf("current symbol table:\n");
      //symtab.print();
      delete [] $2;
    }
  | VECTOR_T SYMBOL DELIM {
      int id;
      id=symtab.lookup($2);
      if (id<0) {
        id=symtab.add($2);
      }
      vartyp[id]=VEC_T;
      //delflag.off(id);
      delflag[id]=0;
      delete [] $2;
    }
  | DELETE LITERAL DELIM {
      int id=symtab.lookup($2);
      if (id>=0) {
        //if (vartyp[id]!=SCALAR_DATA) delflag.on(id);
        if (vartyp[id]!=SCALAR_DATA) delflag.set(id);
        vartyp[id]=DELETED;
      } else {
        yyerror("Variable not defined");
      }
      delete [] $2;
    }
  | DELETE LITERAL LITERAL DELIM {
      int id=symtab.lookup($2);
      if (id>=0) {
        //if (vartyp[id]!=SCALAR_DATA) delflag.on(id);
        if (vartyp[id]!=SCALAR_DATA) delflag.set(id);
        fprintf(stderr, "Variable, %s, marked for deletion\n", $2);
        if (strcmp($3, "m")!=0) vartyp[id]=DELETED;
      } else {
        yyerror("Variable not defined");
      }
      delete [] $2;
      delete [] $3;
    }
  | DELETE ALL DELIM {
      for (int i=0; i<symtab.entries(); i++) {
        //delflag.on(i);
        delflag.set(i);
        vartyp[i]=DELETED;
      }
    }
  | SAVE LITERAL DELIM {
      int id=symtab.lookup($2);
      if (id<0) {
        yyerror("Variable not defined");
      } else {
        //delflag.off(id);
        delflag[id]=0;
      }
      delete [] $2;
    }
  | SAVE ALL DELIM {
      for (int i=0; i<symtab.entries(); i++) {
        //if (vartyp[(long)i]!=DELETED) delflag.off(i);
        if (vartyp[(long)i]!=DELETED) delflag[i]=0;
      }
    };


vector_exp: VEC {
      if ($$.data==NULL) {
        yyclearin;
        scalc_istream=fmemopen(sc_stringbuf, sc_buflen, "r");
        sc_firstpass=2;
      } else {
        $$=$1;
      }
    }
  | LEFT_BRAC vector_exp RIGHT_BRAC {$$=$2;}
  | MINUS vector_exp %prec UMINUS {
      $$.n=$2.n;
      $$.data=new float[$$.n];
      for (integer i=0; i<$$.n; i++) $$.data[i]=-$2.data[i];
      delete [] $2.data;
    }
  | matrix_exp PROD vector_exp {
      integer m, n;
      $1->dimensions(m, n);
      if (n!=$3.n) {
        fprintf(stderr, "Vector multiply: inner dimensions don't match\n");
        fprintf(stderr, "  matrix [%d x %d] * vector [%d]\n", m, n, $3.n);
        yyerror("Range error");
        delete [] $3.data;
        delete $1;
        YYERROR;
      }
      $$.n=m;
      $$.data = $1->vect_mult($3.data);
      delete [] $3.data;
      delete $1;
    }
  | vector_exp PROD matrix_exp {
      integer m, n;
      $3->dimensions(m, n);
      if (n!=$1.n) {
        fprintf(stderr, "Vector multiply: inner dimensions don't match\n");
        fprintf(stderr, "  matrix [%d x %d] * vector [%d]\n", m, n, $1.n);
        yyerror("Range error");
        delete [] $1.data;
        delete $3;
        YYERROR;
      }
      $$.n=m;
      $$.data = $3->left_mult($1.data);
      delete [] $1.data;
      delete $3;
    }
  | vector_exp PLUS vector_exp {
      if ($1.n!=$3.n) {
        fprintf(stderr, "Vector addition: dimensions don't match\n");
        fprintf(stderr, "  vector [%d] + vector [%d]\n", $1.n, $3.n);
        yyerror("Range error");
        delete [] $1.data;
        delete [] $3.data;
        YYERROR;
      }
      $$.n=$1.n;
      $$.data=new float[$$.n];
      for (integer i=0; i<$$.n; i++) $$.data[i]=$1.data[i]+$3.data[i];
      delete [] $1.data;
      delete [] $3.data;
    }
  | scalar_exp PROD vector_exp {
      $$.n=$3.n;
      $$.data=new float[$$.n];
      for (integer i=0; i<$$.n; i++) $$.data[i]=$1*$3.data[i];
      delete [] $3.data;
    }
  | vector_exp PROD scalar_exp {
      $$.n=$1.n;
      $$.data=new float[$$.n];
      for (integer i=0; i<$$.n; i++) $$.data[i]=$3*$1.data[i];
      delete [] $1.data;
    }
  | vector_exp MOD scalar_exp {
      //this is important for generating subscripts:
      int quotient;
      for (integer i=0; i<$$.n; i++) {
        quotient=$1.data[i]/$3;
        $$.data[i]=$1.data[i]-$3*quotient;
      }
    }
  | scalar_exp PLUS vector_exp {
      $$.n=$3.n;
      $$.data=new float[$$.n];
      for (integer i=0; i<$$.n; i++) $$.data[i]=$1+$3.data[i];
      delete [] $3.data;
    }
  | vector_exp PLUS scalar_exp {
      $$.n=$1.n;
      $$.data=new float[$$.n];
      for (integer i=0; i<$$.n; i++) $$.data[i]=$3+$1.data[i];
      delete [] $1.data;
    }
  | vector_exp MINUS vector_exp {
      if ($1.n!=$3.n) {
        fprintf(stderr, "Vector addition: dimensions don't match\n");
        fprintf(stderr, "  vector [%d] + vector [%d]\n", $1.n, $3.n);
        yyerror("Range error");
        delete [] $1.data;
        delete [] $3.data;
        YYERROR;
      }
      $$.n=$1.n;
      $$.data=new float[$$.n];
      for (integer i=0; i<$$.n; i++) $$.data[i]=$1.data[i]-$3.data[i];
      delete [] $1.data;
      delete [] $3.data;
    }
  | scalar_exp MINUS vector_exp {
      $$.n=$3.n;
      $$.data=new float[$$.n];
      for (integer i=0; i<$$.n; i++) $$.data[i]=$1-$3.data[i];
      delete [] $3.data;
    }
  | vector_exp MINUS scalar_exp {
      $$.n=$1.n;
      $$.data=new float[$$.n];
      for (integer i=0; i<$$.n; i++) $$.data[i]=$1.data[i]-$3;
      delete [] $1.data;
    }
  | matrix_exp LEFT_SUB scalar_exp RIGHT_SUB {
      real *vec;
      integer m, n;

      $1->dimensions(m, n);
      if ($3 >= m || $3 < 0) {
        fprintf(stderr, "Subscript, %g, out-of-range [0, %d)\n", $3, m);
        yyerror("Range error");
        delete $1;
        YYERROR;
      }
      $$.data=(*$1)((int) $3);
      $$.n=n;
      delete $1;
    }
  | SOLVE LEFT_BRAC matrix_exp COMMA vector_exp COMMA LITERAL RIGHT_BRAC {
      integer m, n;
      $3->dimensions(m, n);
      if (n!=$5.n) {
        fprintf(stderr, "Solver: outer dimensions don't match\n");
        fprintf(stderr, "  matrix [%d x %d]^(-1) * vector [%d]\n", m, n, $5.n);
        yyerror("Range error");
        delete $3;
        delete [] $5.data;
        delete [] $7;
        YYERROR;
      }
      $$.n=m;
      $$.data=scalc_solve($3, $5.data, $7);
      delete $3;
      delete [] $5.data;
      delete [] $7;
    }
  | INVERT LEFT_BRAC vector_exp RIGHT_BRAC {
      $$.n=$3.n;
      $$.data=$3.data;
      for (integer i=0; i<$$.n; i++) $$.data[i]=1/$$.data[i];
      delete [] $3.data;
    }
  | EMPTY LEFT_BRAC scalar_exp RIGHT_BRAC {
      $$.n=(integer) $3;
      $$.data=new real[$$.n];
      for (integer i=0; i<$$.n; i++) $$.data[i]=0;
    }
  | READ LEFT_BRAC VECTOR_T COMMA LITERAL RIGHT_BRAC {
      if (read_vector($5, &$$, ASCII)!=0) {
        yyerror("I/O error");
        delete [] $5;
        YYERROR;
      }
      delete [] $5;
    }
  | scalar_exp RANGE scalar_exp {
      integer n;
      n=(integer) $3 - (integer) $1;
      if (n<0) {
        $$.n=1-n;
        $$.data=new real[$$.n];
        for (integer i=0; i<$$.n; i++) $$.data[i]=$1-i;
      } else {
        $$.n=n+1;
        $$.data=new real[$$.n];
        for (integer i=0; i<$$.n; i++) $$.data[i]=i+$1;
      }
    }
  | vector_exp LEFT_SUB vector_exp RIGHT_SUB {
      $$.n=$3.n;
      $$.data=new real[$$.n];
      for (integer i=0; i<$3.n; i++) {
        if ($3.data[i] < 0 || $3.data[i] >= $1.n) {
          fprintf(stderr, "Subscript, %7.0g, out-of-range [%d, %d)\n", $3.data[i], 0, $1.n);
          yyerror("Range error");
          delete [] $1.data;
          delete [] $3.data;
          YYERROR;
        }
        $$.data[i]=$1.data[(integer) $3.data[i]];
      }
      delete [] $1.data;
      delete [] $3.data;
    }   
  | matrix_exp LEFT_SUB vector_exp COMMA vector_exp RIGHT_SUB {
      integer m, n;
      $1->dimensions(m, n);
      if ($3.n != $5.n) {
        fprintf(stderr, "Subscripts must have same dimensions\n");
        yyerror("Range error");
        delete $1;
        delete [] $3.data;
        delete [] $5.data;
        YYERROR;
      }
      $$.n=$3.n;
      $$.data=new real[$$.n];
      for (integer i=0; i<$3.n; i++) {
        if ($3.data[i] < 0 || $3.data[i] >= m || $5.data[i] < 0 || $5.data[i] >= n) {
          fprintf(stderr, "Subscript, [%7.0g, %7.0g], out-of-range\n", $3.data[i], $5.data[i]);
          yyerror("Range error");
          delete $1;
          delete [] $3.data;
          delete [] $5.data;
          YYERROR;
        }
        $$.data[i]=(*$1)($3.data[i], $5.data[i]);
      }
      delete $1;
      delete [] $3.data;
      delete [] $5.data;
    }
  | VECTOR_T LEFT_BRAC SYMBOL RIGHT_BRAC {
      char *fname;
      long id;
     
      fname=add_path($3);
      if (read_vector(fname, &$$)!=0) {
        yyerror("I/O error");
        delete [] $3;
        delete [] fname;
        YYERROR;
      }
      delete [] fname;
      id=symtab.lookup($3);
      if (id < 0) {
        id=symtab.add($3);
      }
      vartyp[id]=VEC_T;
      //delflag.off(id);
      delflag[id]=0;
      delete [] $3;
    }
  | SIZE LEFT_BRAC scalar_exp RIGHT_BRAC {
      $$.n=5;
      $$.data=new real[5];
      $$.data[0]=1;
      $$.data[1]=1;
      $$.data[2]=1;
      $$.data[3]=SCALAR_DATA;
      $$.data[4]=1;
    }
  | SIZE LEFT_BRAC vector_exp RIGHT_BRAC {
      $$.n=5;
      $$.data=new real[5];
      $$.data[0]=$3.n;
      $$.data[1]=1;
      $$.data[2]=$3.n;
      $$.data[3]=VEC_T;
      $$.data[4]=$3.n;;
      delete [] $3.data;
    }
  | SIZE LEFT_BRAC matrix_exp RIGHT_BRAC {
      integer m, n;
      integer mt;
      $3->dimensions(m, n);
      $$.n=5;
      $$.data=new real[5];
      $$.data[0]=m*n;
      $$.data[1]=m;
      $$.data[2]=n;
      mt=matrix_type($3);
      $$.data[3]=mt;
      switch (mt) {
        case (SPARSE_T):
          $$.data[4]=((sparse_t *) $3)->size();
          break;
        case (SPARSE_ARRAY_T):
          $$.data[4]=((sparse_array_t *) $3)->nel();
          break;
        default:
          $$.data[4]=m*n;
      }
      delete $3;
    };

matrix_exp: MATRIX {
      if ($$==NULL) {
        yyclearin;
        scalc_istream=fmemopen(sc_stringbuf, sc_buflen, "r");
        sc_firstpass=2;
      } else {
        $$=$1;
      }
    }
  | LEFT_BRAC matrix_exp RIGHT_BRAC {$$=$2;}
  | MINUS matrix_exp %prec UMINUS {
      $2->scal_mult(-1);
      $$=$2;
    }
  | matrix_exp PROD matrix_exp {
      $$=$1->mat_mult($3);
      delete $1;
      delete $3;
      if ($$==NULL) {
        yyerror("Range error");
        YYERROR;
      }
    }
  | matrix_exp PROD scalar_exp {
      $1->scal_mult($3);
      $$=$1;
    }
  | scalar_exp PROD matrix_exp {
      $3->scal_mult($1);
      $$=$3;
    }
  | matrix_exp PLUS matrix_exp {
      $$=$1->add($3);
      delete $1;
      delete $3;
      if ($$==NULL) {
        yyerror("Range error");
        YYERROR;
      }
    }
  | matrix_exp MINUS matrix_exp {
      $3->scal_mult(-1);
      $$=$1->add($3);
      delete $1;
      delete $3;
      if ($$==NULL) {
        yyerror("Range error");
        YYERROR;
      }
    }
  | matrix_exp PLUS scalar_exp {
      integer m, n;
      full_t *idmat;
      $1->dimensions(m, n);
      idmat=new full_t(m, n);
      idmat->ones();
      idmat->scal_mult($3);
      $$=$1->add(idmat);
      delete $1;
      delete idmat;
    }
  | scalar_exp PLUS matrix_exp {
      integer m, n;
      full_t *idmat;
      $3->dimensions(m, n);
      idmat=new full_t(m, n);
      idmat->ones();
      idmat->scal_mult($1);
      $$=$3->add(idmat);
      delete $3;
      delete idmat;
    }
  | matrix_exp MINUS scalar_exp {
      integer m, n;
      full_t *idmat;
      $1->dimensions(m, n);
      idmat=new full_t(m, n);
      idmat->ones();
      idmat->scal_mult(-$3);
      $$=$1->add(idmat);
      delete $1;
      delete idmat;
    }
  | scalar_exp MINUS matrix_exp {
      integer m, n;
      full_t *idmat;
      $3->dimensions(m, n);
      $3->scal_mult(-1);
      idmat=new full_t(m, n);
      idmat->ones();
      idmat->scal_mult($1);
      $$=$3->add(idmat);
      delete $3;
      delete idmat;
    }
  | matrix_exp PLUS vector_exp {
      integer m, n, n1;
      sparse_t *idmat;
      $1->dimensions(m, n);
      if (m < n) n1=m; else n1=n;
      if ($3.n != n1) {
        yyerror("Range error");
        delete $1;
        delete [] $3.data;
        YYERROR;
      }
      idmat=new sparse_t(m, n);
      idmat->identity();
      for (integer i=0; i<n1; i++) idmat->cel($3.data[i], i, i);
      $$=$1->add(idmat);
      delete idmat;
      delete $1;
      delete [] $3.data;
    }
  | vector_exp PLUS matrix_exp {
      integer m, n, n1;
      sparse_t *idmat;
      $3->dimensions(m, n);
      if (m < n1) n1=m; else n1=n;
      if ($1.n != n1) {
        yyerror("Range error");
        delete $3;
        delete [] $1.data;
        YYERROR;
      }
      idmat=new sparse_t(m, n);
      idmat->identity();
      for (integer i=0; i<n1; i++) idmat->cel($1.data[i], i, i);
      $$=$3->add(idmat);
      delete idmat;
      delete $3;
      delete [] $1.data;
    }
  | matrix_exp MINUS vector_exp {
      integer m, n, n1;
      sparse_t *idmat;
      $1->dimensions(m, n);
      if (m < n) n1=m; else n1=n;
      if ($3.n != n1) {
        yyerror("Range error");
        delete $1;
        delete [] $3.data;
        YYERROR;
      }
      idmat=new sparse_t(m, n);
      for (integer i=0; i<n1; i++) idmat->cel(-$3.data[i], i, i);
      $$=$1->add(idmat);
      delete idmat;
      delete $1;
      delete [] $3.data;
    }
  | vector_exp MINUS matrix_exp {
      integer m, n, n1;
      sparse_t *idmat;
      $3->dimensions(m, n);
      $3->scal_mult(-1);
      if (m < n1) n1=m; else n1=n;
      if ($1.n != n1) {
        yyerror("Range error");
        delete $3;
        delete [] $1.data;
        YYERROR;
      }
      idmat=new sparse_t(m, n);
      for (integer i=0; i<n1; i++) idmat->cel($1.data[i], i, i);
      $$=$3->add(idmat);
      delete idmat;
      delete $3;
      delete [] $1.data;
    }
  | TRANSPOSE LEFT_BRAC matrix_exp RIGHT_BRAC {
      $3->transpose();
      $$=$3;
    }
  | IDENTITY LEFT_BRAC scalar_exp COMMA scalar_exp RIGHT_BRAC {
      sparse_t *id;
      id=new sparse_t($3, $5);
      id->identity();
      $$=id;
    }
  | EVD LEFT_BRAC matrix_exp COMMA scalar_exp COMMA LITERAL RIGHT_BRAC {
      sparse_solver_parm<matrix_t> parm;
      set_solver_parm(&parm);
      $$=eig_de($3, $5, parm.ncv, 0, $7, parm.maxiter, parm.tol);
      if ($$==NULL) {
        yyerror("Solver error");
        YYERROR;
      }
      delete $3;
    }
  | SVD LEFT_BRAC matrix_exp COMMA scalar_exp COMMA LITERAL RIGHT_BRAC {
      sparse_solver_parm<matrix_t> parm;
      set_solver_parm(&parm);
      $$=eig_de($3, $5, parm.ncv, 1, $7, parm.maxiter, parm.tol);
      if ($$==NULL) {
        yyerror("Solver error");
        YYERROR;
      }
      delete $3;
    }
  | EVD LEFT_BRAC matrix_exp COMMA scalar_exp RIGHT_BRAC {
      sparse_solver_parm<matrix_t> parm;
      set_solver_parm(&parm);
      $$=eig_de($3, $5, parm.ncv, 0, "LM", parm.maxiter, parm.tol);
      if ($$==NULL) {
        yyerror("Solver error");
        YYERROR;
      }
      delete $3;
    }
  | SVD LEFT_BRAC matrix_exp COMMA scalar_exp RIGHT_BRAC {
      sparse_solver_parm<matrix_t> parm;
      set_solver_parm(&parm);
      $$=eig_de($3, $5, parm.ncv, 1, "LM", parm.maxiter, parm.tol);
      if ($$==NULL) {
        yyerror("Solver error");
        YYERROR;
      }
      delete $3;
    }
  | EVD LEFT_BRAC matrix_exp RIGHT_BRAC {
      sparse_solver_parm<matrix_t> parm;
      set_solver_parm(&parm);
      $$=eig_de($3, parm.nev, parm.ncv, 0, "LM", parm.maxiter, parm.tol);
      if ($$==NULL) {
        delete $3;
        YYERROR;
      }
      delete $3;
    }
  | SVD LEFT_BRAC matrix_exp RIGHT_BRAC {
      sparse_solver_parm<matrix_t> parm;
      set_solver_parm(&parm);
      $$=eig_de($3, parm.nev, parm.ncv, 1, "LM", parm.maxiter, parm.tol);
      if ($$==NULL) {
        delete $3;
        YYERROR;
      }
      delete $3;
    }
  | matrix_exp CPROD vector_exp {
      sparse_array_t *sa;
      integer m, n;
      $1->dimensions(m, n);
      if (n!=$3.n) {
        fprintf(stderr, "Vector multiply: inner dimensions don't match\n");
        fprintf(stderr, "  matrix [%d x %d] # vector [%d]\n", m, n, $3.n);
        yyerror("Range error");
        delete $1;
        delete [] $3.data;
        YYERROR;
      }
      if (matrix_type($1) != SPARSE_ARRAY_T) {
        sa=new sparse_array_t($1);
        delete $1;
      } else {
        sa=(sparse_array_t *) $1;
      }
      $$=sa->cmult($3.data);
      delete sa;
      delete [] $3.data;
    }
  | vector_exp CPROD matrix_exp {
      sparse_array_t *sa;
      integer m, n;
      $3->dimensions(m, n);
      if (n!=$1.n) {
        fprintf(stderr, "Vector multiply: inner dimensions don't match\n");
        fprintf(stderr, "  vector [%d] # matrix [%d x %d]\n", $1.n, m, n);
        yyerror("Range error");
        delete $3;
        delete [] $1.data;
        YYERROR;
      }
      if (matrix_type($3) != SPARSE_ARRAY_T) {
        sa=new sparse_array_t($3);
        delete $3;
      } else {
        sa=(sparse_array_t *) $3;
      }
      $$=sa->left_cmult($1.data);
      delete sa;
      delete [] $1.data;
    }
  | INVERT LEFT_BRAC matrix_exp COMMA LITERAL RIGHT_BRAC {
      $$=scalc_invert($3, $5);
      delete $3;
      delete [] $5;
      if ($$==NULL) YYERROR;
    }
  | EMPTY LEFT_BRAC scalar_exp COMMA scalar_exp RIGHT_BRAC {
      integer m, n;
      m=(integer) $3;
      n=(integer) $5;
      $$=new sparse_t(m, n);
    }
  | READ LEFT_BRAC MATRIX_T COMMA LITERAL RIGHT_BRAC {
      $$=read_matrix($5, $3, ASCII);
      if ($$==NULL) {
        yyerror("I/O error");
        delete [] $5;
        YYERROR;
      }
      delete [] $5;
    }
  | MATRIX_T LEFT_BRAC matrix_exp RIGHT_BRAC {
      $$=convert_matrix($3, $1);
      if ($$==NULL) {
        delete $3;
        YYERROR;
      }
      delete $3;
    }
  | MATRIX_T LEFT_BRAC vector_exp RIGHT_BRAC {
      $$=convert_vector($3.data, $3.n, $1);
      if ($$==NULL) {
        delete [] $3.data;
        YYERROR;
      }
      delete [] $3.data;
    }
  | MATRIX_T LEFT_BRAC SYMBOL RIGHT_BRAC {
      char *fname;
      long id;
     
      fname=add_path($3);
      $$=read_matrix(fname, $1);
      delete [] fname;
      if ($$==NULL) {
        yyerror("I/O error");
        delete [] $3;
        YYERROR;
      }
      id=symtab.lookup($3);
      if (id < 0) {
        id=symtab.add($3);
      }
      vartyp[id]=$1;
      //delflag.off(id);
      delflag[id]=0;
      delete [] $3;
    }
  | matrix_exp LEFT_BRAC scalar_exp RIGHT_BRAC {
      if (matrix_type($1)==SPARSE_ARRAY_T) {
        $$=((sparse_array_t *) $1)->get_el($3);
        delete $1;
        if ($$==NULL) {
          yyerror("Range error");
          YYERROR;
        }
      } else {
        yyerror("Syntax error");
        delete $1;
        YYERROR;
      }
    }
  | matrix_exp LEFT_BRAC vector_exp RIGHT_BRAC {
      sparse_t *sel[$3.n];
      if (matrix_type($1)==SPARSE_ARRAY_T) {
        for (integer i=0; i<$3.n; i++) {
          sel[i]=((sparse_array_t *) $1)->get_el($3.data[i]);
          if (sel[i]==NULL) {
            yyerror("Range error");
            delete $1;
            delete [] $3.data;
            //very elegant...
            for (integer j=0; j<i; j++) delete sel[j];
            YYERROR;
          }
        }
        $$=new sparse_array_t(sel, $3.n);
        delete $1;
        delete [] $3.data;
      } else {
        yyerror("Syntax error");
        delete $1;
        delete [] $3.data;
        YYERROR;
      }
    }

scalar_exp: SCALAR {$$=$1;}
  | LEFT_BRAC scalar_exp RIGHT_BRAC {$$=$2;}
  | MINUS scalar_exp %prec UMINUS {
      $$=-$2;
    }
  | vector_exp PROD vector_exp {
      if ($1.n != $3.n) {
        yyerror("Range error");
        delete [] $1.data;
        delete [] $3.data;
        YYERROR;
      }
      $$=0;
      for (integer i=0; i<$1.n; i++) $$+=$1.data[i]*$3.data[i];
      delete [] $1.data;
      delete [] $3.data;
    }
  | vector_exp LEFT_SUB scalar_exp RIGHT_SUB {
      if ($3 >= $1.n) {
        fprintf(stderr, "Subscript, %d, out-of-range [0, %d)\n", (integer) $3, $1.n);
        yyerror("Range error");
        delete [] $1.data;
        YYERROR;
      }
      $$=$1.data[(int) $3];
      delete [] $1.data;
    }
  | scalar_exp PLUS scalar_exp {
      $$=$1+$3;
    }
  | scalar_exp MINUS scalar_exp {
      $$=$1-$3;
    }
  | scalar_exp PROD scalar_exp {
      $$=$1*$3;
    }
  | scalar_exp QUOTIENT scalar_exp {
      $$=$1/$3;
    }
  | scalar_exp MOD scalar_exp {
      int quotient;
      quotient=$1/$3;
      $$=$1-$3*quotient;
    }
  | NORMBRAC scalar_exp NORMBRAC {
      $$=fabs($2);
    }
  | NORMBRAC vector_exp NORMBRAC {
      $$=0;
      for (integer i=0; i<$2.n; i++) $$+=$2.data[i]*$2.data[i];
      $$=sqrt($$);
      delete [] $2.data;
    }
  | NORMBRAC matrix_exp NORMBRAC {
      integer m, n;
      $$=$2->norm();
      delete $2;
    }
  | INVERT LEFT_BRAC scalar_exp RIGHT_BRAC {
      $$=1./$3;
    }
  | SUM LEFT_BRAC vector_exp RIGHT_BRAC {
      //do we need this?  just use dot prod w/ all ones...
      $$=0;
      for (integer i=0; i<$3.n; i++) $$+=$3.data[i];
    }
  | matrix_exp LEFT_SUB scalar_exp COMMA scalar_exp RIGHT_SUB {
      $$=(*$1)($3, $5);
      delete $1;
    }
  | SQRT LEFT_BRAC scalar_exp RIGHT_BRAC {
      $$=sqrt($3);
    }
  | SIN LEFT_BRAC scalar_exp RIGHT_BRAC {
      $$=sin($3);
    }
  | COS LEFT_BRAC scalar_exp RIGHT_BRAC {
      $$=cos($3);
    }
  | TAN LEFT_BRAC scalar_exp RIGHT_BRAC {
      $$=tan($3);
    }
  | scalar_exp POW scalar_exp {
      $$=pow($1, $3);
    };

script:
  SYMBOL DELIM {
    if (run_script($1)!=0) {
      yyerror("I/O error");
    }
  }

assignment:
  SYMBOL LEFT_SUB scalar_exp RIGHT_SUB ASSIGN scalar_exp DELIM {
      //scalar subscript of a vector
      long id;
      int dflag;
      char *fname;
      integer vt;
      int err;
      id=symtab.lookup($1);

      if (id<0) {
        fprintf(stderr, "Symbol, %s, is not defined\n", $1);
        yyerror("Undefined variable");
      } else {
        vt=vartyp[id];
        dflag=delflag[id];
        if (vt==DELETED) {
          fprintf(stderr, "Symbol, %s, has been deleted\n", $1);
          yyerror("Undefined variable");
        } else if (vt==SCALAR_DATA) {
          fprintf(stderr, "Scalar variable, %s, cannot accept subscripts\n", $1);
          yyerror("Type mismatch");
        } else if (vt==VEC_T) {
          vector_t vec;
          fname=add_path($1);
          err=read_vector(fname, &vec, BIN);

          if (err!=0) {
            fprintf(stderr, "File, %s, has become corrupted, deleting %s\n", fname, $1);
            yyerror("I/O error");
            vartyp[id]=DELETED;
          } else {
            if ($3 < 0 || $3 >= vec.n) {
              fprintf(stderr, "Subscript, %7.0g, out-of-range [%d, %d)\n", $3, 0, vec.n);
              yyerror("Range error");
            } else if (sc_firstpass!=1) {
              vec.data[(int) $3]=$6;
              vector_assign($1, vec.data, vec.n);
            }
            delete [] vec.data;
          }
          delete [] fname;
        } else { 
          fprintf(stderr, "Row substitution not allowed--use vector subscripts instead\n");
          yyerror("Type mismatch");
        }
        //if (dflag==0) delflag.off(id);
        if (dflag==0) delflag[id]=0;
      }
      delete [] $1;
    }
  | SYMBOL LEFT_SUB vector_exp RIGHT_SUB ASSIGN vector_exp DELIM {
      //vector subscript of a vector
      long id;
      int dflag;
      char *fname;
      integer vt;
      int err;
      id=symtab.lookup($1);
      if (id<0) {
        fprintf(stderr, "Symbol, %s, is not defined\n", $1);
        yyerror("Undefined variable");
      } else {
        vt=vartyp[id];
        dflag=delflag[id];		//has variable been set to be deleted?
        if (vt==DELETED) {
          fprintf(stderr, "Symbol, %s, has been deleted\n", $1);
          yyerror("Undefined variable");
        } else if (vt==SCALAR_DATA) {
          fprintf(stderr, "Scalar variable, %s, cannot accept subscripts\n", $1);
          yyerror("Type mismatch");
        } else if (vt==VEC_T) {
          vector_t vec;
          fname=add_path($1);
          err=read_vector(fname, &vec, BIN);
          if (err!=0) {
            fprintf(stderr, "File, %s, has become corrupted, deleting %s\n", fname, $1);
            yyerror("I/O error");
            vartyp[id]=DELETED;
          } else {
            if ($3.n != $6.n) {
              fprintf(stderr, "LHS subscript and RHS must have same dimensions: %d vs. %d\n", $3.n, $6.n);
              yyerror("Range error");
            } else if (sc_firstpass!=1) {
              for (integer i=0; i<$3.n; i++) {
                if ($3.data[i] < 0 || $3.data[i] >= vec.n) {
                  fprintf(stderr, "Subscript, %7.0g, out-of-range [%d, %d)\n", $3.data[i], 0, vec.n);
                  fprintf(stderr, "Ignoring\n");
                  continue;
                }
                vec.data[(int) $3.data[i]]=$6.data[i];
              }
              vector_assign($1, vec.data, vec.n);
            }
            delete [] vec.data;
          }
          delete [] fname;
        } else { 
          fprintf(stderr, "Row substitution not allowed--use paired vector subscripts instead\n");
          yyerror("Type mismatch");
        }
        //if (dflag==0) delflag.off(id);
        if (dflag==0) delflag[id]=0;
      }
      delete [] $1;
      delete [] $3.data;
      delete [] $6.data;
    }
  | SYMBOL LEFT_SUB vector_exp RIGHT_SUB ASSIGN scalar_exp DELIM {
      //vector subscript of a vector with a scalar RHS (RHS is duplicated...)
      long id;
      int dflag;
      char *fname;
      integer vt;
      int err;
      id=symtab.lookup($1);
      if (id<0) {
        fprintf(stderr, "Symbol, %s, is not defined\n", $1);
        yyerror("Undefined variable");
      } else {
        vt=vartyp[id];
        dflag=delflag[id];
        if (vt==DELETED) {
          fprintf(stderr, "Symbol, %s, has been deleted\n", $1);
          yyerror("Undefined variable");
        } else if (vt==SCALAR_DATA) {
          fprintf(stderr, "Scalar variable, %s, cannot accept subscripts\n", $1);
          yyerror("Type mismatch");
        } else if (vt==VEC_T) {
          vector_t vec;
          fname=add_path($1);
          err=read_vector(fname, &vec, BIN);
          if (err!=0) {
            fprintf(stderr, "File, %s, has become corrupted, deleting %s\n", fname, $1);
            yyerror("I/O error");
            vartyp[id]=DELETED;
          } else if (sc_firstpass!=1) {
            for (integer i=0; i<$3.n; i++) {
              if ($3.data[i] < 0 || $3.data[i] >= vec.n) {
                fprintf(stderr, "Subscript, %7.0g, out-of-range [%d, %d)\n", $3.data[i], 0, vec.n);
                fprintf(stderr, "Ignoring\n");
                continue;
              }
              vec.data[(int) $3.data[i]]=$6;
            }
            vector_assign($1, vec.data, vec.n);
            delete [] vec.data;
          }
          delete [] fname;
        } else { 
          fprintf(stderr, "Row substitution not allowed--use paired vector subscripts instead\n");
          yyerror("Type mismatch");
        }
        //if (dflag==0) delflag.off(id);
        if (dflag==0) delflag[id]=0;
      }
      delete [] $1;
      delete [] $3.data;
    }
  | SYMBOL LEFT_SUB scalar_exp COMMA scalar_exp RIGHT_SUB ASSIGN scalar_exp DELIM {
      //scalar subscript of a matrix
      long id;
      int dflag;
      char *fname;
      integer vt;
      int err;
      matrix_t *mat;
      integer m, n;

      id=symtab.lookup($1);
      if (id<0) {
        fprintf(stderr, "Symbol, %s, is not defined\n", $1);
        yyerror("Undefined variable");
      } else {
        vt=vartyp[id];
        dflag=delflag[id];
        if (vt==DELETED) {
          fprintf(stderr, "Symbol, %s, has been deleted\n", $1);
          yyerror("Undefined variable");
        } else if (vt==SCALAR_DATA) {
          yyerror("Syntax error");
        } else if (vt==VEC_T) {
          yyerror("Syntax error");
        } else if (vt==SPARSE_ARRAY_T) {
          fprintf(stderr, "Subscripted sparse arrays cannot be assigned to.\n");
          yyerror("Type mismatch");
        } else {
          fname=add_path($1);
          mat=read_matrix(fname, vt, BIN);
          if (mat==NULL) {
            fprintf(stderr, "File, %s, has become corrupted, deleting %s\n", fname, $1);
            yyerror("I/O error");
            vartyp[id]=DELETED;
          } else {
            mat->dimensions(m, n);
            if ($3 < 0 || $3 >= m || $5 < 0 || $5 >= n) {
              yyerror("Range error");
            } else if (sc_firstpass!=1) {
              mat->cel($8, $3, $5);
              matrix_assign($1, mat);
            }
            delete mat;
          }
          delete [] fname;
        }
        //if (dflag==0) delflag.off(id);
        if (dflag==0) delflag[id]=0;
      } 
      delete [] $1;
    }
  | SYMBOL LEFT_SUB vector_exp COMMA vector_exp RIGHT_SUB ASSIGN vector_exp DELIM {
      //vector subscript of a matrix
      long id;
      int dflag;
      char *fname;
      integer vt;
      int err;
      matrix_t *mat;
      integer m, n;

      //structured programming is all well and good, but a goto might actually
      //be more readable...
      id=symtab.lookup($1);
      if (id<0) {
        fprintf(stderr, "Symbol, %s, is not defined\n", $1);
        yyerror("Undefined variable");
      } else {
        vt=vartyp[id];
        dflag=delflag[id];
        if (vt==DELETED) {
          fprintf(stderr, "Symbol, %s, has been deleted\n", $1);
          yyerror("Undefined variable");
        } else if (vt==SCALAR_DATA) {
          yyerror("Syntax error");
        } else if (vt==VEC_T) {
          yyerror("Syntax error");
        } else if (vt==SPARSE_ARRAY_T) {
          fprintf(stderr, "Subscripted sparse arrays cannot be assigned to.\n");
          yyerror("Type mismatch");
        } else {
          fname=add_path($1);
          mat=read_matrix(fname, vt, BIN);
          if (mat==NULL) {
            fprintf(stderr, "File, %s, has become corrupted, deleting %s\n", fname, $1);
            yyerror("I/O error");
            vartyp[id]=DELETED;
          } else if (sc_firstpass!=1) {
            mat->dimensions(m, n);
            if ($3.n != $5.n || $3.n != $8.n) {
              fprintf(stderr, "LHS subscripts and RHS must have same dimensions: %d vs. %d vs. %d\n", $3.n, $5.n, $8.n);
              yyerror("Range error");
            } else {
              for (integer i=0; i<$3.n; i++) {
                if ($3.data[i] < 0 || $3.data[i] >= m || $5.data[i] < 0 || $5.data[i] >= n) {
                  fprintf(stderr, "Subscript, [%7.0g, %7.0g], out-of-range\n", $3.data[i], $5.data[i]);
                  fprintf(stderr, "Ignoring...\n");
                  continue;
                }
                mat->cel($8.data[i], $3.data[i], $5.data[i]);
              }
            }
            matrix_assign($1, mat);
            delete mat;
          }
          delete [] fname;
        }
        //if (dflag==0) delflag.off(id);
        if (dflag==0) delflag[id]=0;
      } 
      //(memory leaks (were) (should be) punishment for making errors...)
      delete [] $1;
      delete [] $3.data;
      delete [] $5.data;
      delete [] $8.data;
    }
  | SYMBOL LEFT_SUB vector_exp COMMA vector_exp RIGHT_SUB ASSIGN scalar_exp DELIM {
      //vector subscript of a matrix with a scalar RHS (RHS is duplicated...)
      long id;
      int dflag;
      char *fname;
      integer vt;
      int err;
      matrix_t *mat;
      integer m, n;

      //structured programming is all well and good, but a goto might actually
      //be more readable...
      id=symtab.lookup($1);
      if (id<0) {
        fprintf(stderr, "Symbol, %s, is not defined\n", $1);
        yyerror("Undefined variable");
      } else {
        vt=vartyp[id];
        dflag=delflag[id];
        if (vt==DELETED) {
          fprintf(stderr, "Symbol, %s, has been deleted\n", $1);
          yyerror("Undefined variable");
        } else if (vt==SCALAR_DATA) {
          yyerror("Syntax error");
        } else if (vt==VEC_T) {
          yyerror("Syntax error");
        } else if (vt==SPARSE_ARRAY_T) {
          fprintf(stderr, "Subscripted sparse arrays cannot be assigned to.\n");
          yyerror("Type mismatch");
        } else {
          fname=add_path($1);
          mat=read_matrix(fname, vt, BIN);
          if (mat==NULL) {
            fprintf(stderr, "File, %s, has become corrupted, deleting %s\n", fname, $1);
            yyerror("I/O error");
            vartyp[id]=DELETED;
          } else if (sc_firstpass!=1) {
            mat->dimensions(m, n);
            for (integer i=0; i<$3.n; i++) {
              if ($3.data[i] < 0 || $3.data[i] >= m || $5.data[i] < 0 || $5.data[i] >= n) {
                fprintf(stderr, "Subscript, [%7.0g, %7.0g], out-of-range\n", $3.data[i], $5.data[i]);
                fprintf(stderr, "Ignoring...\n");
                continue;
              }
              mat->cel($8, $3.data[i], $5.data[i]);
            }
            matrix_assign($1, mat);
          }
          delete mat;
          delete [] fname;
        }
        //if (dflag==0) delflag.off(id);
        if (dflag==0) delflag[id]=0;
      } 
      //(memory leaks (were) (should be) punishment for making errors...)
      delete [] $1;
      delete [] $3.data;
      delete [] $5.data;
      sc_firstpass=0;
    }
  | SYMBOL ASSIGN vector_exp DELIM {
      //simple vector assignment
      int err=0;
      err=vector_assign($1, $3.data, $3.n);
      delete [] $3.data;
      delete [] $1;
      if (err!=0) {
        yyerror("I/O error");
      }
    }
  | SYMBOL ASSIGN matrix_exp DELIM {
      //simple matrix assignment
      int err=0;
      err=matrix_assign($1, $3);
      delete $3;
      delete [] $1;
      if (err!=0) yyerror("I/O error");
    }
  | SYMBOL ASSIGN scalar_exp DELIM {
      //simple scalar assignment
      scalar_assign($1, $3);
      printf("%s = %g\n", $1, $3);
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
        //delflag.off(id);
        delflag[id]=0;
        break;
      case ('M'):
        id=symtab.add(optarg);
        vartyp[id]=FULL_T;
        //delflag.off(id);
        delflag[id]=0;
        break;
      case ('a'):
        id=symtab.add(optarg);
        vartyp[id]=SPARSE_ARRAY_T;
        //delflag.off(id);
        delflag[id]=0;
        break;
      case ('V'):
        id=symtab.add(optarg);
        vartyp[id]=VEC_T;
        //delflag.off(id);
        delflag[id]=0;
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

