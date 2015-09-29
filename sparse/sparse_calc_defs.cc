#include <string.h>
#include <sys/stat.h>
#include <unistd.h>

#include "error_codes.h"

#include "av.h"
#include "full_util.h"
#include "sparse_calc_defs.h"

extern "C" {
  int yyerror(const char *);
}

using namespace libpetey;
using namespace libsparse;

namespace sparse_calc {

symbol_table<string_petey> symtab;
vector_s<integer> vartyp(SC_NSYM);
vector_s<real> scaltab(SC_NSYM);
//bit_array delflag(SC_NSYM);
std::bitset<SC_NSYM> delflag;
char *path;

linked_list<FILE *> sc_stream_stack;
FILE *scalc_istream;

int64_t scalc_session_id;

int sc_scriptflag=0;

int sc_makeflag=0;
int sc_firstpass=0;
time_t sc_mod0;
time_t sc_mod1;

char *sc_stringbuf=NULL;
int sc_buflen=0;

time_t fmoddate(char *fname) {
  struct stat sf; 
  if (stat(fname, &sf)!=0) return -1.;
  return sf.st_mtime;
}

void set_secondpass() {
  scalc_istream=fmemopen(sc_stringbuf, sc_buflen, "r");
  sc_firstpass=2;
}

void main_help_screen() {
  char *sym;
  printf("\nSparse matrix calculator version 0\n\n");
  printf("USAGE\n");
  printf("  sparse_calc [-d] [-p path] [-h] [-L histfile] [-e script] \n");
  printf("              [-A nA] [-v nev] [-t tol] [-I maxiter]\n");
  printf("\n");
  printf("COMMAND LINE OPTIONS\n");
  printf("  -d  debug mode\n");
  printf("  -p  data path\n");
  printf("  -h  print this help screen\n");
  printf("  -A  set variable, NARNOLDI,\n");
  printf("      -- number of Arnoldi vectors in eigenvalue analyses\n");
  printf("  -v  set variable, NEV, -- number of eigenvalues to calculate\n");
  printf("  -t  set variable, TOL, -- tolerance in solve routines\n");
  printf("  -I  set variable, MAXNITER,\n");
  printf("      -- maximum number of iterations in solve routines\n");
  printf("  -L  load history file\n");
  printf("  -e  execute file before starting up\n");
  printf("\n");
  printf("Accepts commands of the following type:\n");
  printf("\n");
  printf("DECLARATION:\n");
  printf("$ <type> <name>\n");
  printf("  where <name> is a filename and <type> is one of the following:\n");
  printf("    full sparse sparse_array vector\n");
  printf("\n");
  printf("ASSIGNMENT:\n");
  printf("$ <lvalue>=<expression>\n");
  printf("  where:\n");
  printf("     <lvalue> is a variable or subscripted variable\n");
  printf("     <expression> is an arithmetic expression comprised of\n");
  printf("       operators, variable names, functions and scalar constants\n");
  printf("  - currently accepted operators are: */+-()[]#|:\n");
  printf("  - built-in functions are:\n");
  printf("    empty identity read\n");
  printf("    transpose solve invert evd svd\n");
  printf("    full sparse sparse_array vector\n");
  printf("\n");
  printf("COMMAND:\n");
  printf("$ <command> [<options>] [<expression>]\n");
  printf("  where:\n");
  printf("    <command> is a built in command\n");
  printf("    <options> is a string literal\n");
  printf("    <expression> is a matrix, vector, scalar or string literal\n");
  printf("  - built in commands are: com delete help path print run sys\n");
  printf("\n");
  printf("Variables used (read from or written to) by built-in functions:\n");
  printf("  MAXNITER NEV TOL NARNOLDI CONDITIONER ERR EVAL\n");
  printf("\n");
  printf("Use:\n");
  printf("$ help <symbol>\n");
  printf("to get help on commands, functions, operators and variables\n");
  printf("\n");
}

//mode:
//1st bit = retrieve by lexical order (0) or by entry order (1)
//2nd bit = print only saved (undeleted) variables
void print_allvar(FILE *fs, int mode) {
  char *sym;
  int id;
  int sflag=mode%2;			//print in order of entry
  int dflag=(mode >> 1)%2;		//print only saved variables
  for (long i=0; i<symtab.entries(sflag==0); i++) {
    if (dflag && delflag[i]) continue;
    if (sflag==0) {
      sym=(char *) symtab.let(i);     //does this produce a memory leak??
      id=symtab.getid(i);
    } else {
      sym=(char *) symtab.get(i);
      id=i;
    }
    if (vartyp[id]==SCALAR_DATA) {
      if ((mode >> 1) % 2 == 0) {
        fprintf(fs,"%s = %g\n", sym, scaltab[id]);
      }
    } else {
      switch (vartyp[id]) {
        case FULL_T: 
          fprintf(fs, "full %s\n", sym);
          break;
        case SPARSE_T: 
          fprintf(fs, "sparse %s\n", sym);
          break;
        case SPARSE_ARRAY_T: 
          fprintf(fs, "sparse_array %s\n", sym);
          break;
        case VEC_T: 
          fprintf(fs, "vector %s\n", sym);
          break;
      }
    }
    delete [] sym;
  }
}

int run_script(const char *filename) {
  if (scalc_istream!=NULL) sc_stream_stack.push(scalc_istream);
  scalc_istream=fopen(filename, "r");
  if (scalc_istream==NULL) {
    //fprintf("sparse_calc: could not open script file, %s\n", filename);
    return -1;
  }
  return 0;
}

void literal_help(char *lit) {
  printf("\n");
  if (strcmp(lit, "cg")==0) {
    printf("\"%s\" ", lit);
    printf("designates the conjugate gradient method in the solver routines.\n");
    printf("  Variables needed:\n");
    printf("    TOL       = desired tolerance [0.001]\n");
    printf("    MAXNITER  = maximum number of iterations [1000]\n");
    printf("    Final error is returned in ERR\n");
  } else if (strcmp(lit, "bcg")==0) {
    printf("\"%s\" ", lit);
    printf("designates the bi-conjugate gradient method in the solver routines.\n");
    printf("  Variables needed:\n");
    printf("    TOL         = desired tolerance [%g]\n", DEF_TOL);
    printf("    MAXNITER    = maximum number of iterations [%d]\n", DEF_MAXNITER);
    printf("    CONDITIONER = conditioning matrix [I]\n");
    printf("    Final error is returned in ERR\n");
  } else if (strcmp(lit, "evd")==0) {
    printf("\"%s\" ", lit);
    printf("designates the eigenvalue decomposition solution method.  See evd (no quotes).\n");
  } else if (strcmp(lit, "svd")==0) {
    printf("\"%s\" ", lit);
    printf("designates the singular value decomposition solution method.\n");
    printf(" See svd (no quotes).\n");
  } else if (strcmp(lit, "vars")==0) {
    print_allvar(stdout);
  } else if (strcmp(lit, "funcs")==0) {
    printf("NAME          PURPOSE                         USAGE\n");
    printf("\n");
    printf("empty         returns empty matrix or vector  empty(<m> [,<n>])\n");
    printf("identity      returns identity matrix         identity(<m>,<n>)\n");
    printf("read          reads from an ASCII file        read(<type> \"<filename>\")\n");
    printf("size          return size and type info       size(<expression>)\n");
    printf("transpose     returns transpose               transpose(<matrix>)\n");
    printf("solve         solves a matrix equation        solve(<matrix>,<vector>,\"method\")\n");
    printf("invert        inverts a matrix                invert(<matrix>,\"method\")\n");
    printf("evd           eigenvalue decomposition        evd(<matrix>,[<nev>])\n");
    printf("svd           singular value decomposition    svd(<sym_mat>,[<nev>])\n");
    printf("full          convert to full matrix          full(<expression>)\n");
    printf("sparse        convert to sparse matrix        sparse(<expression>)\n");
    printf("sparse_array  convert to sparse array         sparse_array(<expression>)\n");
    printf("vector        convert to vector               vector(<expression>)\n");
  } else if (strcmp(lit, "commands")==0) {
    printf("NAME          PURPOSE\n");
    printf("\n");
    printf("help          get help\n");
    printf("run           runs a file\n");
    printf("print         print an expression to stdout or to a file\n");
    printf("com           comment string\n");
    printf("sys           system command\n");
    printf("path          change data path\n");
    printf("delete        deletes a variable\n");
  } else if (strcmp(lit, "operators")==0) {
    printf("SYMBOL   PURPOSE                                  PRECEDENCE\n");
    printf("\n");
    printf("()       override normal order of operations/       2\n");
    printf("           extract elements from a sparse array     2\n");
    printf("[]       subscript operator                         3\n");
    printf("\"\"       quote string literal                       -\n");
    printf("||       vector and matrix norm                     -\n");
    printf(":        range operator                             3\n");
    printf("*        product operator                           4\n");
    printf("#        cumulative product                         4\n");
    printf("/        division                                   4\n");
    printf("+        addition operator                          5\n");
    printf("-        subtraction                                5/1\n");
    printf("=        assignment                                 -\n");
  } else if (strcmp(lit, "types")==0) {
    printf("TYPE NAME        DESCRIPTION\n");
    printf("\n");
    printf("full             full matrix\n");
    printf("sparse           sparse matrix\n");
    printf("sparse_array     array of sparse matrices (treated as single matrix)\n");
    printf("vector           vector type\n");
    printf("scalar           scalar (floating point) type\n");
  } else if (strcmp(lit, "reserved")==0) {
    printf("Here is the full list of reserved words:\n");
    printf("\n");
    printf("- COMMANDS:    com delete help path print run sys\n");
    printf("- TYPES:       full sparse sparse_array vector scalar\n");
    printf("- OPERATORS:   ()[]\"|:*/#+-=\n");
    printf("- FUNCTIONS:   empty identity read transpose solve invert\n");
  } else {
    int id=symtab.lookup(lit);

    if (id >= 0) {
      printf("%s ", lit);
      if (vartyp[id]==SCALAR_DATA) {
        printf("=%g\n", scaltab[id]);
      } else {
        switch (vartyp[id]) {
          case FULL_T: 
            printf(" is defined as a full matrix\n");
            break;
          case SPARSE_T: 
            printf(" is defined as a sparse matrix\n");
            break;
          case SPARSE_ARRAY_T: 
            printf(" is defined as an array of sparse matrices\n");
            break;
          case VEC_T: 
            printf(" is defined as a vector\n");
            break;
        }
      }
    }

    if (strcmp(lit, "EVAL")==0) {
      printf("Used to store returned eigenvalues\n");
    } else if (strcmp(lit, "TOL")==0) {
      printf("Determines the tolerance for iterated solution methods.\n");
      printf("Can be set from the command line using -t\n");
    } else if (strcmp(lit, "ERR")==0) {
      printf("Used to return the error from a solver.\n");
    } else if (strcmp(lit, "MAXNITER")==0) {
      printf("Determines the maximum number of iterations for iterated solution methods\n");
      printf("Can be set from the command line using -I\n");
    } else if (strcmp(lit, "NEV")==0) {
      printf("Determines the number of eigenvalues/eigenvectors extracted.\n");
      printf("Can be set from the command line using -v\n");
    } else if (strcmp(lit, "NARNOLDI")==0) {
      printf("Determines the number of Arnoldi vectors for ARPACK based operations.\n");
      printf("Can be set from the command line using -A\n");
    } else if (strcmp(lit, "CONDITIONER")==0) {
      printf("Can be used to store the conditioning matrix for solution methods that use one.\n");
    } else if (id<0) {
      printf("%s is not in the help database\n", lit);
    }
  }
  printf("\n");
}

void type_help(integer type_id) {
  printf("\n");
  switch (type_id) {
    case FULL_T:
      printf("Full matrix type.  Use in a declaration as follows:\n");
      printf("$ full <variable>\n");
      printf("  - Matrix is read in from a binary file of the same name.\n\n");
      printf("Can also be used for type conversion, as in:\n");
      printf("  <matrix>=full(<expression>)\n");
      printf("  - where <expression> is a vector or matrix.\n");
      printf("\n");
      printf("or as a \'shortcut\' declaration:\n");
      printf("  <matrix>=full(<name>)\n");
      printf("  - simultaneously declares the variable and returns its contents\n");
      printf("- Declared variables persist after exit.\n");
      break;
    case SPARSE_T:
      printf("Sparse matrix type.  Use in a declaration as follows:\n");
      printf("  $ sparse <variable>\n");
      printf("  - Matrix is read in from a binary file of the same name.\n\n");
      printf("Can also be used for type conversion, as in:\n");
      printf("  <matrix>=sparse(<expression>)\n");
      printf("  - where <expression> is a vector or matrix.\n");
      printf("\n");
      printf("or in a \'shortcut\' declaration:\n");
      printf("  <matrix>=sparse(<name>)\n");
      printf("  - simultaneously declares the variable and returns its contents\n");
      printf("- Declared variables persist after exit.\n");
      break;
    case SPARSE_ARRAY_T:
      printf("Array of sparse matrices type.  Use in a declaration as follows:\n");
      printf("  $ sparse_array <variable>\n");
      printf("  - Matrix is read in from a binary file of the same name.\n\n");
      printf("Can also be used for type conversion, as in:\n");
      printf("  <matrix>=sparse_array(<expression>)\n");
      printf("  - where <expression> is a vector or matrix. Must be square.\n");
      printf("\n");
      printf("or in a \'shortcut\' declaration:\n");
      printf("  <matrix>=sparse_array(<name>)\n");
      printf("  - simultaneously declares the variable and returns its contents\n");
      printf("- Declared variables persist after exit.\n");
      break;
  }
  printf("\n");
}

char *add_path(const char *sym) {
  char *fname;
  fname=new char[strlen(sym)+strlen(path)+1];
  sprintf(fname, "%s%s", path, sym);
  return fname;
}

matrix_t * read_matrix(const char *fname, int type, int mode) {
  FILE *fs;
  matrix_t *result;
  int err;

  if (sc_firstpass==1) {
    result=new matrix_t();
  } else {
    switch (type) {
      case(SPARSE_T): 
        result=new sparse_t();
        fprintf(stderr, "Reading file %s as sparse\n", fname);
        break;
      case(SPARSE_ARRAY_T): 
        result=new sparse_array_t();
        fprintf(stderr, "Reading file %s as sparse array\n", fname);
        break;
      case(FULL_T):
        result=new full_t();
        fprintf(stderr, "Reading file %s as full matrix\n", fname);
        break;
      default: 
        fprintf(stderr, "%s is not a matrix\n", fname);
        return NULL;
    }
  }
  fs=fopen(fname, "r");
  if (fs==NULL) {
    fprintf(stderr, "sparse_calc: unable to open file, %s, for reading\n", fname);
    return NULL;
  }
  if (mode==BIN) {
    if (sc_firstpass==1) result->read(fs, type); else result->read(fs);
  } else {
    err=result->scan(fs);
    if (err!=0) {
      fprintf(stderr, "sparse_calc: error reading ASCII file, %s\n", fname);
      yyerror("I/O error");
    }
  }
  fclose(fs);
  return result;
}

int read_vector(const char *fname, vector_t *vec, int mode) {
  FILE *fs;
  char **line;
  long n;
  switch (mode) {
    case(ASCII): {
      fprintf(stderr, "Reading ASCII file %s as vector\n", fname);
      line=read_ascii_all(fname, &n);
      if (line==NULL) {
        fprintf(stderr, "sparse_calc: unable to open file, %s, for reading\n", fname);
        return UNABLE_TO_OPEN_FILE_FOR_READING;
      }
      vec->data=new real[n];
      vec->n=n;
      for (int i=0; i<n; i++) {
        sscanf(line[i], "%g", vec->data+i);
        delete [] line[i];
      }
      delete [] line;
      return 0;
    }
    case(BIN): {
      fprintf(stderr, "Reading binary file %s as vector\n", fname);
      fs=fopen(fname, "r");
      if (fs==NULL) {
        fprintf(stderr, "Unable to open file, %s, for reading\n", fname);
        return UNABLE_TO_OPEN_FILE_FOR_READING;
      }
      fseek(fs, 0, SEEK_END);
      n=(ftell(fs)/sizeof(real));
      fseek(fs, 0, SEEK_SET);
      vec->data=new real[n];
      vec->n=n;
      if (sc_firstpass!=1) fread(vec->data, sizeof(real), n, fs);
      fclose(fs);
      return 0;
    }
    default: {
      fprintf(stderr, "%s mode not recognized\n", fname);
      return -1;
    }
  }

}

int scalar_assign(const char *symbol, real data) {
  int id;
  id=symtab.lookup(symbol);
  if (id<0) {
    id=symtab.add(symbol);
  }
  vartyp[id]=SCALAR_DATA;
  scaltab[id]=data;
  //no need to "delete" scalar variables:
  //delflag.off(id);
  delflag[id]=0;
  return 0;
}

int matrix_assign(const char *symbol, matrix_t *matrix) {
  int id;
  FILE *fs;
  int err=0;
  char *fname;
  char *str;
  integer type;

  id=symtab.lookup(symbol);
  type=matrix_type(matrix);

  if (id < 0) {
    id=symtab.add(symbol);
  }

  vartyp[id]=type;
  //delflag.on(id);
  delflag[id]=1;
  fname=add_path(symbol);

  if (sc_firstpass==1) {
    matrix_t dum;
    integer m, n, m1, n1;
    fs=fopen(fname, "r");
    dum.read(fs, type);
    dum.dimensions(m, n);
    matrix->dimensions(m1, n1);
    if (m!=m1 || n!=n1) {
      //how do we indicate that the statement needs to be executed for real?
      sc_firstpass=2;
    }
  } else { 
    fs=fopen(fname, "w");
    matrix->write(fs);
  }

  fclose(fs);
  delete [] fname;

  return err;
}

matrix_t *convert_matrix(matrix_t *matrix, integer t2) {
  matrix_t *output;

  //if (sc_firstpass==1) return matrix;

  switch (t2) {
    case (SPARSE_T): {
        output = new sparse_t(matrix);
        break;
      }
    case (SPARSE_ARRAY_T): {
        output = new sparse_array_t(matrix);
        break;
      }
    case (FULL_T): {
        output = new full_t(matrix);
        break;
      }
    default: {
        yyerror("Type mismatch");
        output=NULL;
      }
  }

  return output;

}

int vector_assign(const char *symbol, real *data, integer n) {
  int id;
  FILE *fs;
  int err=0;
  char *fname;
  matrix_t *output;
  id=symtab.lookup(symbol);


  if (id < 0) {
    id=symtab.add(symbol);
  }
  vartyp[id]=VEC_T;
  //delflag.on(id);
  delflag[id]=1;

  fname=add_path(symbol);
  fs=fopen(fname, "w");
  fwrite(data, sizeof(real), n, fs);
  fclose(fs);

  delete [] fname;

  return err;
}

matrix_t *convert_vector(real *data, integer n, integer t2) {
  matrix_t *output;
  if (t2==FULL_T) {
    output=new full_t(&data, 1, n);
    output->transpose();
  } else if (t2==SPARSE_T) {
    output=new sparse_t(&data, 1, n);
    //output->print(stdout);
    output->transpose();
    //output->print(stdout);
  } else if (t2==SPARSE_ARRAY_T) {
    sparse_t *s;
    s=new sparse_t(&data, 1, n);
    s->transpose();
    //note: this is not a memory vvv leak (at least, I don't think so...)
    output=new sparse_array_t(&s, 1);
  } else {
    yyerror("Type mismatch/syntax error");
    output=NULL;
  }
  return output;
}

matrix_t * eig_de(matrix_t *mat, 	//matrix on which to perform decomposition
		int32_t nv,		//number of eigenvalues
		int32_t na,		//number of Arnoldi vectors
		int ty,			//type (evd, svd)
		const char which[2],
		int32_t maxiter,
		real tol)
{
  int32_t m, n;
  full_t *v1;
  real **v;
  real *s;
  std::complex<real> **v2;
  std::complex<real> *s1;

  mat->dimensions(m, n);
  if (m!=n) {
    fprintf(stderr, "Eigenvalue decomposition only valid for square matrices\n");
    yyerror("Range error");
    return NULL;
  }
  if (ty==0) {
    matrix_t *ev;
    real **s2;
    //s1=new std::complex<real>[nv];
    //v2=cc_arevd(m, nv, na, s1, mat, which, maxiter, tol);
    //s=new real[na*2];
    //for (integer i=0; i<nv; i++) {
    //  s[i]=std::real(s1[i]);
    //  s[i+nv]=std::imag(s1[i]);
    //}
    v=allocate_matrix<real, integer>(na, n);
    s2=allocate_matrix<real, integer>(2, nv);
    FORTRAN_FUNC(sarevd)(&n, &nv, &na, v[0], s2[0], &maxiter, &tol, which, (void *) mat);
    ev=new full_t(s2, 2, nv);
    ev->transpose();
    matrix_assign("EVAL", ev);
    v1=new full_t(v, nv, n);
    delete ev;
    delete [] s2[0];
    delete [] s2;
    delete [] v[0];
    delete [] v;
    //delete [] s1;
    //delete [] v2[0];
    //delete [] v2;
  } else {
    s=new real[nv];
    v=cc_arsvd(m, nv, na, s, mat, which, maxiter, tol);
    vector_assign("EVAL", s, nv);
    v1=new full_t(v, nv, n);
    delete [] s;
    delete [] v[0];
    delete [] v;
  }

  return v1;
}

void set_solver_parm(sparse_solver_parm<matrix_t> *parm) {

  int id;
  id=symtab.lookup("MAXNITER");
  if (id>=0 && vartyp[id]==SCALAR_DATA) {
    parm->maxiter=scaltab[id];
  } else {
    parm->maxiter=DEF_MAXNITER;
  }
  id=symtab.lookup("NEV");
  if (id>=0 && vartyp[id]==SCALAR_DATA) {
    parm->nev=scaltab[id];
  } else {
    parm->nev=DEF_NEV;
  }
  id=symtab.lookup("NARNOLDI");
  if (id>=0 && vartyp[id]==SCALAR_DATA) {
    parm->ncv=scaltab[id];
  } else {
    parm->ncv=DEF_NARNOLDI;
  }
  id=symtab.lookup("TOL");
  if (id>=0 && vartyp[id]==SCALAR_DATA) {
    parm->tol=scaltab[id];
  } else {
    parm->tol=DEF_TOL;
  }

  //let the calling program do that:
  parm->cond=NULL;

}

matrix_t *get_cond(const char *sym, integer m, integer n) {
  int id;
  char *fname;
  matrix_t *cond;

  id=symtab.lookup(sym);
  if (id >= 0 && (vartyp[id]==FULL_T || vartyp[id]==SPARSE_T || vartyp[id]==SPARSE_ARRAY_T)) {
    fname=add_path(sym);
    cond=read_matrix(fname, vartyp[id]);
    delete [] fname;
  } else {
    yyerror("invert: warning, conditioning matrix not found, using identity");
    cond=new sparse_t(m, n);
    ((sparse_t *) cond)->identity();
  }
  return cond;
}

real *get_x0(const char *sym, matrix_t *a, real *b) {
  integer m, n;
  vector_t x0;
  char *fname;
  int id, err;

  a->dimensions(m, n);
  id=symtab.lookup(sym);
  if (id>=0 && vartyp[id]==VEC_T) {
    fname=add_path(sym);
    err=read_vector(fname, &x0);
    if (err==0 && x0.n==n) return x0.data;
    delete [] fname;
  }

  x0.data=new real[n];
  x0.n=n;
  for (integer i=0; i<m && i<n; i++) {
    real val;
    val=(*a)(i, i);
    if (val!=0) x0.data[i]=b[i]/val; else x0.data[i]=1;
  }

  //vector_assign(sym, x0.data, n);

  return x0.data;
}

matrix_t *scalc_invert(matrix_t *mat, char *type) {
  sparse_solver_parm<matrix_t> parm;
  int (*solver) (matrix_t &, real *, real *, void *);
  int id;
  integer m, n;
  real **result1;
  full_t *result2;
  int err;

  mat->dimensions(m, n);
  if (m!=n) {
    yyerror("Square matrices only supported for inverse function");
    return NULL;
  }

  set_solver_parm(&parm);

  if (strcmp(type, "cg")==0) {
    solver=&conj_grad<integer, real, matrix_t>;
  } else if (strcmp(type, "bcg")==0) {
    solver=&biconj_grad<integer, real, matrix_t>;
    parm.cond=get_cond("CONDITIONER", m, n);
  } else if (strcmp(type, "evd")==0) {
    result1=sparse_evd_invert<integer, real, matrix_t>(*mat, (void *) &parm);
  } else if (strcmp(type, "svd")==0) {
    result1=sparse_svd_invert<integer, real, matrix_t>(*mat, (void *) &parm);
  } else {
    yyerror("invert: error, method not recognized");
    return NULL;
  }
  if (err!=0) yyerror("invert: solver returned non-zero error code");

  if (strcmp(type, "cg")==0 || strcmp(type, "bcg")==0) {
    result1=allocate_matrix<real, integer>(n, m);
    err=sparse_invert<integer, real, matrix_t>(*mat, result1, (void *) &parm, solver);
    scalar_assign("ERR", parm.err);
  }

  result2=new full_t(result1, n, m);

  delete_matrix(result1);
  if (parm.cond!=NULL) delete parm.cond;

  return result2;
  
}

real *scalc_solve(matrix_t *mat, real *b, char *type) {
  sparse_solver_parm<matrix_t> parm;
  int (*solver) (matrix_t &, real *, real *, void *);
  int id;
  integer m, n;
  real *result;
  char *fname;
  int err;

  mat->dimensions(m, n);

  if (m!=n) {
    yyerror("Square matrices only supported for solver");
    return NULL;
  }

  set_solver_parm(&parm);

  if (strcmp(type, "cg")==0 || strcmp(type, "n:cg")==0) {
    solver=&conj_grad<integer, real, matrix_t>;
    result=get_x0("X0", mat, b);
  } else if (strcmp(type, "bcg")==0 || strcmp(type, "n:bcg")==0) {
    solver=&biconj_grad<integer, real, matrix_t>;
    result=get_x0("X0", mat, b);
    parm.cond=get_cond("CONDITIONER", m, n);
  } else if (strcmp(type, "evd")==0 || strcmp(type, "n:evd")==0) {
    solver=&sparse_evd_solver<integer, real, matrix_t>;
    result=new real[m];
  } else if (strcmp(type, "svd")==0 || strcmp(type, "n:svd")==0) {
    solver=&sparse_svd_solver<integer, real, matrix_t>;
    result=new real[m];
  } else {
    yyerror("solve: method not recognized");
    return NULL;
  }

  if (type[0]=='n') {
    err=sparse_normal<integer, real, matrix_t>(*mat, b, result, (void *) &parm, solver);
    if (err!=0) yyerror("invert: solver returned non-zero error code");
  } else {
    err=(*solver) (*mat, b, result, (void *) &parm);
    if (err!=0) yyerror("invert: solver returned non-zero error code");
  }

  if (strcmp(type, "cg")==0 || strcmp(type, "bcg")==0 ||
	strcmp(type, "n:cg")==0 || strcmp(type, "n:bcg")==0) {
    scalar_assign("ERR", parm.err);
  }
  if (parm.cond!=NULL) delete parm.cond;

  return result;
  
}

matrix_t * sub_mat(matrix_t *mat, matrix_t *s1, matrix_t *s2) {
  integer m1, m2, n1, n2;
  //all dimen

}

template <class index_t, class scalar>
integer matrix_type(matrix_base<index_t, scalar> *mat) {
  integer result;
  if (typeid(*mat) == typeid(sparse<index_t, scalar>)) {
    result=SPARSE_T;
  } else if (typeid(*mat) == typeid(sparse_array<index_t, scalar>)) {
    result=SPARSE_ARRAY_T;
  } else if (typeid(*mat) == typeid(full_matrix<index_t, scalar>)) {
    result=FULL_T;
  } else {
    result=-1;
  }
  return result;
}

template integer matrix_type(matrix_base<integer, real> *);

} //end namespace sparse_calc

