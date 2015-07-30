//"bit array" class:

#ifndef BIT_ARRAY_H_INCLUDED
#define BIT_ARRAY_H_INCLUDED 1

namespace libpetey {

  typedef unsigned int word;

  class bit_array {
    protected:
      long nbits;		//number of bits stored
      long nwords;	//number of words in data array
      word *data;		//the data

    public:
      //constructors:
      bit_array();
      bit_array(long n);
      bit_array(long n, char value);
      bit_array(char *d, long n);
      bit_array(word *d, long nw, long nb);

      //copy constructor:
      bit_array(const bit_array &other);

      //equals operator:
      bit_array operator = (const bit_array &other);

      //type conversion:
      operator char *();
      operator short *();
      operator long *();

      //manipulate the array:
      void resize(long n);

      //destructor:
      ~bit_array();
    
      //change the values in the array:
      char on(long ind);
      char off(long ind);
      char flip(long ind);
      char change(char value, long ind);

      //access an element of the array:
      char operator [] (long ind);
      //number of nonzero elements up to index value:
      long nnonzero (long ind);

      //how big is it?
      long nel();

      void print();

  };

  inline long bit_array::nel() {
    return nbits;
  }

  inline char bit_array::change(char value, long ind) {
    if (value == 0) return off(ind); else return on(ind);
  }

} //end namespace libpetey

#endif

