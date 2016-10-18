#include <stdio.h>
#include <vector>

#include "string_petey.h"
#include "symbol_table.h"

using namespace std;
using namespace libpetey;

class grammar_t;
class parse_tree;

class expr_t {
  public:
    expr_t();
    virtual ~expr_t();

    virtual parse_tree * parse(char *string, int ver)=0;
    virtual int parse(char *string, parse_tree *pt, int depth=0)=0;

    virtual int nver()=0;
    virtual int ns(int v)=0;
    virtual char *nameof()=0;
};

class primitive: public expr_t {
  protected:
    char val;
  public:
    primitive(char v1);
    virtual ~primitive();

    virtual parse_tree * parse(char *string, int ver);
    virtual int parse(char *string, parse_tree *pt, int depth=0);

    virtual int nver();
    virtual int ns(int v);

    virtual char *nameof();

};

class compound: public expr_t {
  protected:
    grammar_t *grammar;		//each expression knows it's own identity
    long id;
    vector<expr_t **> sub;	//list of sub-expressions 
    vector<int> nsub;		//number of sub-expressions per version
    int ntyp;			//number of versions

  public:
    compound(grammar_t *g, long i);
    virtual ~compound();

    int add(expr_t **neex, int n);
    virtual int nver();
    virtual int ns(int v);

    virtual char *nameof();

    virtual parse_tree * parse(char *string, int ver);
    virtual int parse(char *string, parse_tree *pt, int depth=0);
};

class parse_tree {
  public:
    expr_t *expr;
    int ver;			//which version of the expression
    parse_tree **sub;		//sub-trees
    char *loc;			//location where expression starts
    int len;			//total length of the expression

    parse_tree();
    parse_tree(expr_t *e, int v, char * l, int n, parse_tree **s);

    ~parse_tree();

    void print(int depth=0);
    char * value();

};

class grammar_t {
  friend class compound;
  protected:
    symbol_table<string_petey> name;
    vector<expr_t *> exprlist;
    long last;		//id of last expression defined
  public:
    grammar_t();
    ~grammar_t();

    int add(char *nm, char *def);	//add a new expression
    int load(FILE *fs);		//load a complete grammar
    parse_tree * parse(char *expr_name, char *string);
    parse_tree *parse(char *string);
};

