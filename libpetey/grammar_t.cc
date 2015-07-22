#include <ctype.h>
#include <stdlib.h>
#include <string.h>

#include "linked.h"
#include "grammar_t.h"
#include "read_ascii_all.h"

//for debugging:
parse_tree *trunk;

parse_tree::parse_tree() {
  expr=NULL;
  ver=0;
  loc=NULL;
  len=0;
  sub=NULL;
}

parse_tree::parse_tree(expr_t *e, int v, char *l, int n, parse_tree **s) {
  expr=e;
  ver=v;
  loc=l;
  sub=s;
  len=n;
}

parse_tree::~parse_tree() {
  if (sub!=NULL) {
    for (int i=0; i<expr->ns(ver); i++) {
      delete sub[i];
    }
    delete [] sub;
  }
}

char *parse_tree::value() {
  char *result;
  result=new char[len+1];
  for (int i=0; i<len; i++) result[i]=loc[i];
  result[len]='\0';
  return result;
}

void parse_tree::print(int depth) {
  char *name;
  char *val=value();

  if (expr==NULL) return;
  name=expr->nameof();

  for (int i=0; i<depth; i++) printf("  ");
  printf("%s %d %s\n", name, ver, val);
  delete [] name;
  delete [] val;

  if (sub!=NULL) {
    for (int i=0; i<expr->ns(ver); i++) {
      if (sub[i]!=NULL) sub[i]->print(depth+1);
    }
  }
}


expr_t::expr_t() {}

expr_t::~expr_t() {}

primitive::primitive(char v1) {
  val=v1;
}

primitive::~primitive() {
}

int primitive::nver() {
  return 1;
}

int primitive::ns(int v) {
  return 0;
}

char * primitive::nameof() {
  char *result=new char[3];
  result[0]='\'';
  result[1]=val;
  result[2]='\0';
  return result;
}

parse_tree *primitive::parse(char *string, int ver) {
  parse_tree *result=NULL;
  char *nm=nameof();

  //printf("Parsing primitive %s with string %s\n", nm, string);

  if (string[0]==val && ver==0) {
    printf("Primitive %s succeeded; returning\n", nm);
    result=new parse_tree(this, 0, string, 1, NULL);
    delete [] nm;
  }
  printf("Primitive %s failed; returning\n", nm);
  delete [] nm;
  return result;
}

int primitive::parse(char *string, parse_tree *pt, int depth) {
  int err=1;
  char *nm=nameof();
  char *dstring;

  dstring=new char[depth*2+1];
  for (int i=0; i<depth*2; i++) dstring[i]=' ';
  dstring[depth*2]='\0';


  printf("%s%s: with string %s\n", dstring, nm, string);
  //printf("\n");

  pt->sub=NULL;
  if (pt->loc==NULL) {
    //printf("Parsing primitive %s with string %s\n", nm, string);
    if (string[0]==val) {
      pt->loc=string;
      pt->len=1;
      pt->ver=0;
      pt->expr=this;
      err=0;
    }
  } else if (string==pt->loc) {
    //printf("Re-parsing primitive %s with string %s (failure)\n", nm, string);
    err=-1;
  } else {
    //printf("Re-parsing primitive %s with string %s\n", nm, string);
    if (string[0]==val) {
      pt->loc=string;
      pt->len=1;
      pt->ver=0;
      pt->expr=this;
      err=0;
    }
  }
  if (err==0) {
    printf("%sPrimitive %s succeeded; returning\n", dstring, nm);
  } else {
    printf("%sPrimitive %s failed; returning\n", dstring, nm);
  }

  trunk->print();

  delete [] nm;
  return err;

}

compound::compound(grammar_t *g, long i) {
  grammar=g;
  ntyp=0;
  id=i;
}

compound::~compound() {
}

int compound::add(expr_t **neex, int n) {
  sub[ntyp]=neex;
  nsub[ntyp]=n;
  ntyp++;

  return ntyp;
}

int compound::nver() {
  return ntyp;
}

int compound::ns(int v) {
  return nsub[v];
}

char * compound::nameof() {
  return grammar->name.get(id);
}

parse_tree * compound::parse(char *string, int ver) {
  parse_tree **sub_tree;
  int v1;
  int i;
  int len;
  char *lnew;
  char *nm=nameof();

  printf("Parsing type %s version %d with string %s\n", nm, ver, string);
  printf("  there are %d sub-expr.\n", nsub[ver]);

  sub_tree=new parse_tree*[nsub[ver]];
  for (int j=0; j<nsub[ver]; j++) sub_tree[j]=NULL;

  for (i=0; i<nsub[ver]; i++) {
    printf("Parsing the %dth sub-expression\n", i);
    if (sub_tree[i]==NULL) {
      v1=0;
    } else {
      v1=sub_tree[i]->ver+1;
      delete sub_tree[i];
      sub_tree[i]=NULL;
    }
    for (int j=v1; j<sub[ver][i]->nver(); j++) {
      if (i>0) {
        lnew=sub_tree[i-1]->loc+sub_tree[i-1]->len;
        printf("version %d; sub_tree[%d]->loc=%s; sub_tree[%d]->len=%d\n", j, i-1, sub_tree[i-1]->loc, i-1, sub_tree[i-1]->len);
      } else {
        lnew=string;
      }
      sub_tree[i]=sub[ver][i]->parse(lnew, j);
      if (sub_tree[i]!=NULL) break;
    }
    if (sub_tree[i]==NULL) {
      i-=2;
      if (i<-1) break;
    }
  }

  if (i<nsub[ver]) {
    printf("Type %s version %d failed; returning\n", nm, ver);
    delete [] sub_tree;
    delete [] nm;
    return NULL;
  }

  printf("Type %s version %d succeeded; returning\n", nm, ver);
  delete [] nm;
  len=0;
  for (int j=0; j<nsub[ver]; j++) len+=sub_tree[j]->len;

  return new parse_tree(this, ver, string, len, sub_tree);
}

int compound::parse(char *string, parse_tree *pt, int depth) {
  int err=0;
  char *lnew;
  int j;
  char *dstring;

  dstring=new char[depth*2+1];
  for (int i=0; i<depth*2; i++) dstring[i]=' ';
  dstring[depth*2]='\0';

  char *nm=nameof();

  printf("%s%s: with string %s\n", dstring, nm, string);
  //printf("\n");

  if (pt->ver >= ntyp) {
    if (pt->sub!=NULL) {
      for (int k=0; k<nsub[pt->ver]; k++) delete pt->sub[k];
      delete [] pt->sub;
    }
    pt->ver=0;
    printf("Type %s failed; returning\n", nm);
    return -2;
  }

  if (pt->sub==NULL) {
    for (; pt->ver < ntyp; pt->ver++) {
      printf("%sParsing version %d with %d sub-expressions\n", dstring, pt->ver, nsub[pt->ver]);
      pt->sub=new parse_tree*[nsub[pt->ver]];
      for (j=0; j<nsub[pt->ver]; j++) pt->sub[j]=new parse_tree();
      for (j=0; j<nsub[pt->ver]; j++) {
        if (j>0) {
          lnew=pt->sub[j-1]->loc+pt->sub[j-1]->len;
        } else {
          lnew=string;
        }
        err=sub[pt->ver][j]->parse(lnew, pt->sub[j], depth+1);
        if (err!=0) {
          j-=2;
          if (j<-1) break;
        }
      }
      if (j<-1) {
        for (int k=0; k<nsub[pt->ver]; k++) delete pt->sub[k];
        delete [] pt->sub;
        pt->sub=NULL;
      } else {
        break;
      }
    }
    if (pt->sub==NULL) {
      err=-1;
    }
  } else {
    printf("%sRe-parsing version %d with %d sub-expressions\n", dstring, pt->ver, nsub[pt->ver]);
    for (j=nsub[pt->ver]-1; j<nsub[pt->ver]; j++) {
      pt->ver=0;
      if (j>0) {
        lnew=pt->sub[j-1]->loc+pt->sub[j-1]->len;
      } else {
        lnew=string;
      }
      err=sub[pt->ver][j]->parse(lnew, pt->sub[j], depth+1);
      if (err!=0) {
        j-=2;
        if (j<-1) break;
      }
    }
    if (j<-1) {
      for (int k=0; k<nsub[pt->ver]; k++) delete pt->sub[k];
      delete [] pt->sub;
      pt->sub=NULL;
      pt->ver++;
      err=parse(string, pt, depth);
    }
  }

  if (err==0) {
    printf("%sType %s succeeded; returning\n", dstring, nm);
    pt->len=0;
    for (int j=0; j<nsub[pt->ver]; j++) pt->len+=pt->sub[j]->len;
    pt->loc=string;
    pt->expr=this;
  } else {
    printf("%sType %s failed; returning\n", dstring, nm);
  }

  pt->expr=this;
  trunk->print();

  delete [] nm;

  return err;
}

grammar_t::grammar_t() {
}

grammar_t::~grammar_t() {
  for (int i=0; i<name.entries(); i++) delete exprlist[i];
}

int grammar_t::add(char *name1, char *def) {
  linked_list<int> start;
  linked_list<int> end;
  char oldc;
  int *ind1=NULL;
  int *ind2=NULL;
  long n1, n2;
  char *name2;
  expr_t **neex=NULL;
  long id, id1;
  int i;
  int err=0;

  printf("%s:%s\n", name1, def);

  //do nothing but find all the starts and ends of words (space-separated character strings):
  oldc=def[0];
  if (isalpha(oldc)) start.add(0);
  for (i=1; def[i]!='\0'; i++) {
    if (isblank(oldc) && isblank(def[i])!=1) {
      start.add(i);
    }
    if (isblank(def[i]) && isblank(oldc)!=1) {
      end.add(i);
    }
    oldc=def[i];
  }
  if (isblank(def[i-1])!=1) end.add(i);

  ind1=start.make_array(n1);
  ind2=end.make_array(n2);

  assert(n1==n2);

  neex=new expr_t *[n1];

  id1=name.lookup(name1);
  if (id1<0) {
    id1=name.add(name1);
    exprlist[id1]=new compound(this, id1);
  }

  for (i=0; i<n1; i++) {
    def[ind2[i]]='\0';
    name2=def+ind1[i];
    printf("%d %s\\\n", i, name2);

    if (name2[0]=='\'') {
      if (name2[1]=='\0') {
        fprintf(stderr, "grammar_t: syntax error, subexpr. %d\n", i);
        fprintf(stderr, "        --single quote must be followed by primitive\n");
        err=-1;
        goto finish;
      }
      if (name2[2]!='\0' && name2[2]!='\'') {
        fprintf(stderr, "grammar_t: syntax error, subexpr. %d\n", i);
        fprintf(stderr, "        --primitives must be enclosed by single quotes or a quote and a space\n");
        err=-1;
        goto finish;
      }
      if (name2[2]=='\'') name2[2]='\0';
      id=name.lookup(name2);
      if (id<0) {
        id=name.add(name2);
        exprlist[id]=new primitive(name2[1]);
      }
    } else if (isalpha(name2[0])==0) {
      fprintf(stderr, "grammar_t: syntax error, subexpr. %d\n", i);
      fprintf(stderr, "        --expression names must begin with a letter\n");
      err=-1;
      goto finish;
    } else {
      for (int j=ind1[i]+1; j<ind2[i]; j++) {
        if (isalnum(def[j])==0 && def[j]!='_') {
          fprintf(stderr, "grammar_t: syntax error, subexpr. %d\n", i);
          fprintf(stderr, "        --expression names must be composed of letters, number or underscores\n");
          err=-1;
          goto finish;
        }
      }
      id=name.lookup(name2);
      if (id<0) {
        fprintf(stderr, "grammar_t: symbol, %s, is undefined\n", name2);
        err=-2;
        goto finish;
      }
    }
    neex[i]=exprlist[id];
  }

  //printf("\n");
  //name.print();
  //printf("\n");

  ((compound *) exprlist[id1])->add(neex, n1);
  last=id1;

  finish:
    if (ind1!=NULL) delete [] ind1;
    if (ind2!=NULL) delete [] ind2;

  return err;
} 

int grammar_t::load(FILE *fs) {
  char **line;
  long nlines;
  int ptr1, ptr2;
  int len;
  char *name;
  char *def;
  int nerr=0;

  line=read_ascii_all(fs, &nlines);

  for (long i=0; i<nlines; i++) {
    for (ptr1=0; isalpha(line[i][ptr1])==0 && line[i][ptr1]!='\n'; ptr1++);
    if (line[i][ptr1]=='\n') continue;
    for (ptr2=ptr1; isblank(line[i][ptr2])==0 && line[i][ptr2]!='\n' && line[i][ptr2]!=':'; ptr2++);
    if (line[i][ptr2]=='\n') {
      fprintf(stderr, "(1) syntax error line %d; skipped\n", i+1);
      nerr++;
      goto skip;
    }

    name=line[i]+ptr1;
    if (line[i][ptr2]!=':') {
      name[ptr2]='\0';
      for (ptr2++; line[i][ptr2]!=':' && line[i][ptr2]!='\n'; ptr2++);
      if (line[i][ptr2]=='\n' || line[i][ptr2]!=':') {
        fprintf(stderr, "(2) syntax error line %d; skipped\n", i+1);
        nerr++;
        goto skip;
      }
    } else {
      name[ptr2]='\0';
    }

    def=line[i]+ptr2+1;
    len=strlen(def);
    if (def[len-1]=='\n') def[len-1]='\0';

    if (add(name, def)!=0) {
      fprintf(stderr, "(3) syntax error line %d; skipped\n", i+1);
      nerr++;
    }

    skip:
      delete [] line[i];
  }

  return nerr;

}

parse_tree * grammar_t::parse(char *expr_name, char *string) {
  long id;
  int nver;
  parse_tree *result;

  id=name.lookup(string);
  if (id<0) return NULL;

  nver=exprlist[id]->nver();
  for (int i=0; i<nver; i++) {
    result=exprlist[id]->parse(string, i);
    if (result!=NULL) return result;
  }

  return result;
}

parse_tree * grammar_t::parse(char *string) {
  long id;
  int nver;
  parse_tree *result;
  int err;

  id=last;
  if (id<0) return NULL;

  trunk=new parse_tree();

  err=exprlist[id]->parse(string, trunk);
  result=trunk;
  if (err!=0) {
    delete result;
    result=NULL;
  }

  return result;
}



