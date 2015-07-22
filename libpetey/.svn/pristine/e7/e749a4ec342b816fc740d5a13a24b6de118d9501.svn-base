#include <sys/types.h>

#include "bit_array.h"

#include "linked.cc"
#include "vector_s.h"
#include "string_petey.h"
#include "symbol_table.h"
#include "sparse_array.h"
#include "full_matrix.h"
#include "read_ascii_all.h"

#include "sparse_solve.h"

//initial space for symbols:
#define SC_NSYM 1000

//maximum size of function stack:
#define MAXFSTACK 1000

using namespace libpetey;
using namespace libsparse;

namespace sparse_calc {

typedef float real;
typedef int32_t integer;

#define SPARSE_T 1
#define SPARSE_ARRAY_T 2
#define FULL_T 3
#define VEC_T 4
#define SCALAR_DATA 5
#define DELETED -1

//read modes:
#define BIN 0
#define ASCII 1

typedef matrix_base<integer, real> matrix_t;
typedef sparse<integer, real> sparse_t;
typedef sparse_array<integer, real> sparse_array_t;
typedef full_matrix<integer, real> full_t;

struct vector_t {
  integer n;
  real *data;
};

#define DEF_NARNOLDI 20
#define DEF_NEV 5

extern symbol_table<string_petey> symtab;
extern vector_s<integer> vartyp;
extern vector_s<real> scaltab;
extern bit_array delflag;
extern char *path;
//extern symbol_table<string_petey> funtab;
//extern int (*) (int nmat, matrix **mat, int nscal, real *scal, void &*result) [MAXNSYM];

//input stream stack:
extern linked_list<FILE *> sc_stream_stack;
extern FILE *scalc_istream;
extern int sc_istr_stckptr;

extern int64_t scalc_session_id;

//if we want the program to exit once scalc_istream ends:
extern int sc_scriptflag;

//for make-like dependency checking:
extern int sc_makeflag;
extern int sc_firstpass;
extern time_t sc_mod0;		//mod. date of script or executable
					//--whichever is sooner
extern time_t sc_mod1;		//mod. date of LHS

//in case we want to read from a string buffer:
extern char *sc_stringbuf;
extern int sc_buflen;

time_t fmoddate(char *fname);
void set_secondpass();

//help functions:
void main_help_screen();
void literal_help(char *lit);
void type_help(integer type_id);

//operations on variables:
char *add_path(const char *sym);
matrix_t * read_matrix(const char *fname, int type, int mode=BIN);
int read_vector(const char *fname, vector_t *vec, int mode=BIN);
int matrix_assign(const char *symbol, matrix_t *matrix);
matrix_t *convert_matrix(matrix_t *old, integer t);
int vector_assign(const char *symbol, real *data, integer n);
matrix_t *convert_vector(real *data, integer n, integer t);
int scalar_assign(const char *symbol, real data);

//system operations:
void print_allvar(FILE *fs, int mode=0);
int run_script(const char *filename);

matrix_t * eig_de(matrix_t *mat,        //matrix on which to perform decomposition
                int32_t nv,             //number of eigenvalues
                int32_t na,             //number of Arnoldi vectors
                int ty,                //type (evd, svd)
		const char which[2],
		int32_t maxiter,
		real tol);
void set_solver_parm(sparse_solver_parm<matrix_t> *parm);
matrix_t *get_cond(const char *sym);
matrix_t *scalc_invert(matrix_t *mat, char *type);
real *scalc_solve(matrix_t *mat, real *b, char *type);

//hmmm....  these have been inlined...
real *sub_mat(matrix_t *mat, real *s1, real *s2, integer n);
real *sub_vec(real *vec, real *sub, integer n);

int lval_mat(matrix_t *mat, real *s1, real *s2, real *vec);
int lval_vec(real *v1, real *sub, real *v2, integer n);

template <class index_t, class scalar>
integer matrix_type(matrix_base<index_t, scalar> *mat);

} //end namespace sparse_calc

