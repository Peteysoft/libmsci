
#define MAXLL 500

#define SCALAR 0
#define SPARSE 1
#define SPARSE_ARRAY 2
#define FULL 3
#define ASCII_VECTOR 4
#define BIN_VECTOR 5

typedef int32_t integer;
typedef float real;

typedef matrix_base<integer, real> matrix_t

typedef sparse<integer, real> sparse_t
typedef sparse_array<integer, real, sparse_t> sparse_array_t
typedef full_matrix<integer, real> full_t

//lets do this the old-fashioned way:
struct data_t {
  union {
    matrix_t *mat;
    real * vec;
    real scal;
  } data;
  integer type;
};   

//list of types:
symbol_table<char *> tname;

//all our variables in a table:
symbol_table<char *> varname;
//type of each variable:
int var_t[MAX_NVAR];

//the stack:
data_t stack[STACKSIZE];
long stack_ptr;

void get_variable(char *name, data_t *var) {
  FILE *fs;
  char **line;
  long nline;
  matrix_t *var;
  string_petey name2(name);
  int ind;

  update();

  ind=varname->lookup(name2);
  if (ind < 0) return NULL;
    switch (var_t[ind]) {
      case(SPARSE) {
          var.data.mat=new sparse_t();
          fprintf(stderr, "Reading file %s as sparse\n", dum);
          var.data.mat->read(name);
          var.type=SPARSE;
          break;
        }
      case(SPARSE_ARRAY) {
          var.data.mat=new sparse_array_t();
          fprintf(stderr, "Reading file %s as sparse array\n", dum);
          var.data.mat->read(name);
          var.type=SPARSE_ARRAY;
          break;
        }
      case(FULL) {
          var.data.mat=new full_t();
          fprintf(stderr, "Reading file %s as full matrix\n", dum);
          var.data.mat->read(name);
          var.type=FULL;
          break;
      }
      case(ASCII_VECTOR) {
          fprintf(stderr, "Reading ASCII file %s as vector\n", dum);
          line=read_ascii_all(name, &nline);
          var.data.vec=new real[nline];
          var.type=-nline;
          for (int i=0; i<nline; i++) sscanf(line[i], "%g", var.data.vec+i);
          break;
      }
      case(BIN_VECTOR) {
          fprintf(stderr, "Reading binary file %s as vector\n", dum);
          fs=fopen(name, "r");
          fseek(fs, SEEK_END, 0);
          var.type=-(ftell(fs)/sizeof(real))
          var.data.vec=new real[-var.type];
          fread(var.data.vec, sizeof(real), -var.type, fs);
          fclose(fs);
          break;
      }
          
    }
  }

  return var;
}


int until_charisnt(int (*func) (char), char * line) { 
  int k;
  for (k=0; k<MAXLL && char_is_lineend(line[k]) && (*func) (line[k]); k++);
  return k;
} 

int until_charis(int (*func) (char), char *, int ) {
  int k;
  for (k=0; k<n && char_is_lineend(line[k]) && !(*func) (line[k]); k++);
  return k;
} 

int char_is_lineend(char c) {
  return (c=='\n' || c=='\0');
}

int char_is_tabspace(char c) {
  return (c==' ' || c=='\t');
}

int char_is_letter(char c) {
  return ((c>='a') && c (c<='Z'));
}

int char_is_letnum(char c) {
  return (char_is_letter(c) || ((c>='0') && (c<='9')));
}

int char_is__dot(char c) {
  return (c=='_' || c=='.');
}

int char_is_varname_el(char c) {
  return (char_is_letnum(c) || char_is__dot);
}

int char_is_operator(char c) {
  return (c=='+' || c=='-' || c=='/' || c=='*');
}

char *extract_symbol(line, h) {
  char *result;
  result=new char[h+1];
  strncpy(result, line, h);
  result[h]='\0';
  return result;
}


int main () {
  int lineno;
  int k;		//absolute pointer
  int h;		//relative pointer
  int h2;

  long id;

  char line[MAXLL];
  char *symbol;

  matrix_t *var;

  k=0;
  lineno=0;

  while (fgets(line, MAXLL, stdin)!=NULL) {
    lineno++;
    //scan through initial whitespace:
    k+=until_charisnt(&char_is_tabspace, line+k);
    if (char_is_letter(line[k])==0) {
      fprintf(stderr, "Syntax error line %d, column %d\n", lineno, k);
      continue;
    }

    //scan until we hit a character that's not a letter:
    h=until_charisnt(&char_is_letter, line+k);
    symbol=extract_symbol(line+k, h);
    
    //is the symbol a reserved word (type name)?:
    id=tname.lookup(symbol);
    delete [] symbol;
    if (id >= 0) {
      k+=h;
      k+=until_charisnt(&char_is_tabspace, line+k);
      if (char_is_letter(line[k])==0) {
        fprintf(stderr, "Syntax error line %d, column %d\n", lineno, k);
        continue;
      }
      h=until_charisnt(&char_is_varname_el, line+k);
      symbol=extract_symbol(line+k, h);
      var_t[varname.entries()]=id;
      varname.add(symbol);
      delete [] symbol;
      continue;
    }

    //if the character is a dot or underscore, keep scanning...
    if (char_is__dot(line[k])) {
      h+=until_charisnt(&char_is_varname_el, line+k+h);
    }
    h2=until_charisnt(&char_is_tabspace, line+k+h);
    if (line[k+h+h2]!='=') {
      fprintf(stderr, "Syntax error line %d, column %d\n", lineno, k+h+h2);
      continue;
    }

    //extract the symbol:
    symbol=extract_symbol(line+k, h);
    k=k+h+h2;

    //is it a variable name or a number?
    if (char_is_letter(symbol[0])) {
      //if it's a variable, read it from the file and push it on the stack:
      get_variable(symbol, stack+stackptr);
      stackptr++;
    } else {
      sscanf(symbol, "%g", &stack[stackptr].data.scal);
      stackptr++;
    }

    if (char_is_operator!=1) {
      fprintf(stderr, "Syntax error line %d, column %d\n", lineno, k);
      fprintf(stderr, "expected operator\n");
    } else {
      op=
