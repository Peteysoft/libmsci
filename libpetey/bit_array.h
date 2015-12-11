//"bit array" class:

#ifndef BIT_ARRAY_H_INCLUDED
#define BIT_ARRAY_H_INCLUDED 1

namespace libpetey {

  typedef unsigned int word;

  //even thought the STL contains a "bitset" class, we've decided to keep this 
  //since the functionality doesn't quite overlap
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
      bit_array & operator = (const bit_array &other);

      //type conversion:
      //operator char *();
      //operator short *();
      //operator long *();
      template <class type>
      operator type *();

      //manipulate the array:
      void resize(long n);

      //destructor:
      ~bit_array();
    
      //change the values in the array:
      char on(long ind);
      char off(long ind);
      char flip(long ind);
      char set(char value, long ind);

      //access an element of the array:
      char operator [] (long ind);
      //number of nonzero elements up to index value:
      long nnonzero (long ind);

      //how big is it?
      long nel();

      void print();

      //set bits randomly (need to call ran_init() before using)
      void random();

      //test the bit array:
      int test(int ntrial);		//number of trials

      //can be treated as numbers:
      int operator == (const bit_array &other);
      int operator > (const bit_array &other);
      int operator < (const bit_array &other);
      int operator >= (const bit_array &other);
      int operator <= (const bit_array &other);
      int operator != (const bit_array &other);

  };

  inline long bit_array::nel() {
    return nbits;
  }

  inline char bit_array::set(char value, long ind) {
    if (value == 0) return off(ind); else return on(ind);
  }

  int test_bit_array(int size,		//size of bit array
		  int ntrial);		//number of trials

} //end namespace libpetey

#endif

