%{
#include "linked_list.h"
#include "multi_class.h"

extern "C" {
  int yylex(void);
  int yyerror(const char *);
}

%}

%union {
  char *string;
  int integer;
  classifier_obj<real_a, cls_ta> *classifier;
  multiclass_hier<real_a, cls_ta> *multi;
  linked<*multiclass_hier<real_a, cls_ta> > *blist;
  linked<cls_ta> *clist;
}

%token <string> FNAME
%token <integer> CLS
%token <multi> branch
%token <classifier> model
%token <blist> branch_list
%token <clist> class_list

branch: model '{' branch_list '}' | CLS ';' {
    $$=new multi_class_hier<real_a, cls_ta>($1);
};

class_list: CLS {
    $$=new linked<cls_ta>();
    $$->add($1);
} | CLS class_list {
    $$->add($1);
};

branch_list: branch {
    $$=new linked<*multiclass_hier<real_a, cls_ta> >();
    $$->add($1);
} | branch branch_list {
    $$->add($1);
};

partition_list: partition | partition partition_list;

model: FNAME {
    $$=new agf2class<real_a, cls_ta>($1);
} | partition_list;

partition: FNAME class_list '|' class_list ';';


