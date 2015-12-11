#include <stdio.h>
#include <assert.h>

#include "bit_array.h"
#include "randomize.h"

namespace libpetey {

bit_array::bit_array() {
  nbits=0;
  nwords=1;
  data=new word[1];
  data[0]=0;
}

bit_array::bit_array(long n) {
  nbits=n;
  nwords=(n-1)/(sizeof(word)*8)+1;
  data=new word[nwords];
  //all unused space is to be padded with zeros:
  for (long i=0; i<nwords; i++) data[i]=0;
}

bit_array::bit_array(long n, char value) {
  nbits=n;
  nwords=(n-1)/(sizeof(word)*8)+1;
  data=new word[nwords];
  if (value == 0) {
    for (long i=0; i<nwords; i++) data[i]=0;
  } else {
    for (long i=0; i<nwords-1; i++) data[i]=255;
    data[nwords-1]=0;		//all extra space should be padded with 0s
    for (long i=0; i<nwords*sizeof(word)*8-n; i++) on(i+(nwords-1)*sizeof(word));
  }
}

bit_array::bit_array(char *d, long n) {
  nbits=n;
  nwords=(n-1)/(sizeof(word)*8)+1;
  data=new word[nwords];
  data[nwords-1]=0;	//all extra space must be padded with 0s
  for (long i=0; i<n; i++) {
    if (d[i] > 0) on(i); else off(i);
  }
}

bit_array::bit_array(word *d, long nw, long nb) {
  assert(nb<=sizeof(word)*8);
  nwords=nw;
  nbits=nb;
  data=new word[nw];
  for (long i=0; i<nw; i++) data[i]=d[nw-i-1];
}

bit_array::bit_array(const bit_array &other) {
  nwords=(other.nbits-1)/sizeof(word)/8+1;
  if (nwords<=0) nwords=1;
  data=new word[nwords];
  nbits=other.nbits;
  for (long i=0; i<nwords; i++) data[i]=other.data[i];
}

bit_array & bit_array::operator = (const bit_array &other) {
  if (nwords>=0) delete [] data;
  nwords=(other.nbits-1)/sizeof(word)/8+1;
  if (nwords<=0) nwords=1;
  data=new word[nwords];
  nbits=other.nbits;
  for (long i=0; i<nwords; i++) data[i]=other.data[i];
  return *this;
}

bit_array::~bit_array() {
  delete [] data;
}

template <class type>
bit_array::operator type * () {
  type *result=new type[nbits];
  for (long i=0; i<nbits; i++) result[i]=(*this)[i];
  return result;
}

template bit_array::operator char *();
template bit_array::operator short *();
template bit_array::operator int *();
template bit_array::operator float *();
template bit_array::operator double *();

void bit_array::resize (long n) {
  long k;
  long nwn;
  word *new_data;

  //calculate new number of words:
  nwn=(n-1)/(sizeof(word)*8)+1;
  if (nwn>nwords && nwn<nwords*2) nwn=nwords*2;

  if (nwn > nwords) {
    //create extended array and copy data to it:
    new_data=new word[nwn];
    for (long i=0; i<nwords; i++) new_data[i]=data[i];
    for (long i=nwords; i<nwn; i++) new_data[i]=0;
    nwords=nwn;
    delete [] data;
    data=new_data;
  }

  nbits=n;

}

char bit_array::on(long ind) {
  long wordind;
  long bitind;
  word tmpl;
 
  //check for out of bounds subscripts:
  if (ind < 0) {
    printf("bit_array::on: negative subscript encountered\n");
    return -1;
  }
  if (ind >= nbits) resize(ind+1);

  //calculate the subscript for the word and subscript for its bit:
  wordind=ind/(sizeof(word)*8);
  bitind=ind % (sizeof(word)*8);

  //left shift a bit to the correct position and do a bitwise or:
  tmpl=1 << bitind;
  data[wordind]=data[wordind] | tmpl;

  return 1;
}

char bit_array::off(long ind) {
  long wordind;
  long bitind;
  word tmpl;
  
  //check for out of bounds subscripts:
  if (ind < 0) {
    printf("bit_array::off: negative subscript encountered\n");
    return -1;
  }
  if (ind >= nbits) resize(ind+1);

  //calculate the subscript for the word and subscript for its bit:
  wordind=ind/(sizeof(word)*8);
  bitind=ind % (sizeof(word)*8);

  //left shift a bit to the correct position and do a bitwise and of its complement:
  tmpl=1 << bitind;
  data[wordind]=data[wordind] & ~tmpl;

  return 0;
}

char bit_array::flip(long ind) {
  long wordind;
  long bitind;
  word tmpl, tmpl2;
  
  //check for out of bounds subscripts:
  if (ind < 0) {
    printf("bit_array::flip: negative subscript encountered\n");
    return -1;
  }
  if (ind >= nbits) resize(ind+1);

  //calculate the subscript for the word and subscript for its bit:
  wordind=ind/(sizeof(word)*8);
  bitind=ind % (sizeof(word)*8);

  tmpl=1 << bitind;
  data[wordind]=data[wordind] ^ tmpl;
}

char bit_array::operator [] (long ind) {
  long wordind;
  long bitind;
  word result1;
  
  //check for out of bounds subscripts:
  if (ind < 0) {
    fprintf(stderr, "bit_array::[]: negative subscript encountered\n");
    return -1;
  }
  if (ind >= nbits) {
    fprintf(stderr, "bit_array::[]: out-of-bounds subscript encountered\n");
    return 0;
  }

  //find the location of the bit:
  wordind=ind/(sizeof(word)*8);
  bitind=ind % (sizeof(word)*8);

  //right shift it until it is in the lowest order position:
  result1=data[wordind] >> bitind;

  //taking the modolo wrt 2 will get it:
  return char (result1 % 2);

}

long bit_array::nnonzero(long ind) {
  long wordind;
  long bitind;
  word tw;
  long n;

  //check for out of bounds subscripts:
  if (ind < 0) {
    printf("Warning: negative subscript encountered\n");
    return -1;
  }
  if (ind >= nbits) {
    printf("Warning: out-of-bounds subscript encountered\n");
    return -1;
  }

  //find the location of the bit:
  wordind=ind/(sizeof(word)*8);
  bitind=ind % (sizeof(word)*8);

  n=0;
  for (long i=0; i<wordind; i++) {
    tw=data[i];
    for (int j=1; j<sizeof(word)*8; j++) {
      n+=tw % 2;
      tw = tw >> 1;
    }
    n+=tw % 2;
  }

  tw=data[wordind];
  for (long i=0; i<bitind; i++) {
    n+=tw % 2;
    tw = tw >> 1;
  }
  n+=tw % 2;

  return n;
}

void bit_array::print() {
  for (long i=0; i<nbits; i++) printf("%d ", (*this)[i]);
  //for (long i=0; i<nwords; i++) printf("%b", data[i]);
  printf("\n");
}

void bit_array::random() {
  for (long i=0; i<nbits; i++) {
    if (2*ranu() < 1) {
      off(i);
    } else {
      on(i);
    }
  }
}

int bit_array::test(int ntrial) {
  char bitstr[nbits];		//bitstring as char array
  for (int i=0; i<ntrial; i++) {
    for (long j=0; j<nbits; j++) {
      if (2*ranu() < 1) {
        off(j);
	bitstr[j]=0;
      } else {
        on(j);
	bitstr[j]=1;
      }
    }
    for (long j=0; j<nbits; j++) {
      flip(j);
      flip(j);
    }
    for (long j=0; j<nbits; j++) {
      if (bitstr[j]!=(*this)[j]) {
        fprintf(stderr, "bit_array::test: failed at bit %d; trial %d\n", j, i);
	fprintf(stderr, "Bit-array: ");
	for (int k=0; k<nbits; k++) fprintf(stderr, "%hhd", bitstr[k]);
	fprintf(stderr, "\n");
	return -1;
      }
    }
  }
  return 0;
}

int bit_array::operator == (const bit_array &other) {
  int nwd;		//number of words containing data
  //once again, do we allow comparison between different sized arrays?
  assert(nbits==other.nbits);

  nwd=(nbits-1)/(sizeof(word)*8)+1;
  //because we've padded extra space with zeros (we hope) comparisons are
  //simpler:
  for (int i=0; i<nwd; i++) {
    if (data[i]!=other.data[i]) return 0;
  }
  return 1;
}

int bit_array::operator > (const bit_array &other) {
  int nwd;		//number of words containing data
  //once again, do we allow comparison between different sized arrays?
  assert(nbits==other.nbits);

  nwd=(nbits-1)/(sizeof(word)*8)+1;
  //because we've padded extra space with zeros (we hope) comparisons are
  //simpler:
  for (int i=0; i<nwd; i++) {
    if (data[i]>other.data[i]) return 1;
    if (data[i]<other.data[i]) return 0;
  }
  return 0;
}

int bit_array::operator < (const bit_array &other) {
  int nwd;		//number of words containing data
  //once again, do we allow comparison between different sized arrays?
  assert(nbits==other.nbits);

  nwd=(nbits-1)/(sizeof(word)*8)+1;
  //because we've padded extra space with zeros (we hope) comparisons are
  //simpler:
  for (int i=0; i<nwd; i++) {
    if (data[i]<other.data[i]) return 1;
    if (data[i]>other.data[i]) return 1;
  }
  return 0;
}

int bit_array::operator >= (const bit_array &other) {
  int nwd;		//number of words containing data
  //once again, do we allow comparison between different sized arrays?
  assert(nbits==other.nbits);

  nwd=(nbits-1)/(sizeof(word)*8)+1;
  //because we've padded extra space with zeros (we hope) comparisons are
  //simpler:
  for (int i=0; i<nwd; i++) {
    if (data[i]<other.data[i]) return 0;
    if (data[i]>other.data[i]) return 1;
  }
  return 1;
}

int bit_array::operator <= (const bit_array &other) {
  int nwd;		//number of words containing data
  //once again, do we allow comparison between different sized arrays?
  assert(nbits==other.nbits);

  nwd=(nbits-1)/(sizeof(word)*8)+1;
  //because we've padded extra space with zeros (we hope) comparisons are
  //simpler:
  for (int i=0; i<nwd; i++) {
    if (data[i]>other.data[i]) return 0;
    if (data[i]<other.data[i]) return 1;
  }
  return 1;
}

int bit_array::operator != (const bit_array &other) {
  int nwd;		//number of words containing data
  //once again, do we allow comparison between different sized arrays?
  assert(nbits==other.nbits);

  nwd=(nbits-1)/(sizeof(word)*8)+1;
  //because we've padded extra space with zeros (we hope) comparisons are
  //simpler:
  for (int i=0; i<nwd; i++) {
    if (data[i]!=other.data[i]) return 1;
  }
  return 0;
}

//test function:
int test_bit_array(int size, int ntrial) {
  bit_array a1(size);
  bit_array *a2;
  bit_array *a3;
  char *s1;
  char *s2;
  int err;
  err=a1.test(ntrial);
  if (err!=0) return err;

  for (int i=0; i<ntrial; i++) {
    a1.random();
    s1=(char *) a1;
    a2=new bit_array(s1, size);
    a1=*a2;
    a3=new bit_array();
    for (int j=0; j<size; j++) a3->set(a1[j], j);
    s2=(char *) *a3;
    for (int j=0; j<size; j++) {
      if (s1[j]!=s2[j]) {
        fprintf(stderr, "test_bit_array: failed at bit %d; trial %d\n", j, i);
	fprintf(stderr, "Bit-array 1: ");
	for (int k=0; k<size; k++) fprintf(stderr, "%hhd", s1[k]);
	fprintf(stderr, "\n");
	fprintf(stderr, "Bit-array 2: ");
	for (int k=0; k<size; k++) fprintf(stderr, "%hhd", s2[k]);
	fprintf(stderr, "\n");
	return -1;
      }
    }
    delete [] s1;
    delete [] s2;
    delete a2;
    delete a3;
  }
  return 0;
}

}
