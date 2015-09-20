#ifndef SYMBOL_TABLE_H_INCLUDED
#define SYMBOL_TABLE_H_INCLUDED

#include <stdio.h>

namespace libpetey {

  template <class sym_t>
  class symbol_table {
    protected:
      long n;			//number of entries
      long nunq;		//number of unique symbols
      long array_size;		//size of array
      sym_t *sym;		//the symbols
      long *sind;		//indexes the variables in terms of the symbols

      int update_flag;

      //for thread-safe code, update must be atomic
      void update();
      
    public:
      symbol_table();
      ~symbol_table();

      long * collect_garbage(long &ndup);

      //adds a new symbol to the table--returns its id
      long add(sym_t name);

      //returns a unique, pervasive (until the death of the object) ID:
      long lookup(sym_t name);

      //looks up symbol based on ID:
      sym_t get(long id);

      //looks up symbol based on lexical ordering:
      sym_t let(long ind);
      //looks up id based on lexical ordering:
      long getid(long ind);

      void print();

      long entries(int uflag=0);

  };

}

#endif

