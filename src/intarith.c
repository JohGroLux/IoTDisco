#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "config.h"
#include "intarith.h"


typedef uint16_t Word;   // single-length word
typedef uint32_t DWord;  // double-length word
typedef int32_t SDWord;  // signed double-length word


// Bitlength of a Word
#define WSIZE (8*sizeof(Word))


/*------check whether a multiprecision integer is 0------*/
int int_is0(const Word *a, int len)
{
  Word word = 0;
  int i;
  
  for (i = len - 1; i >= 0; i--) word |= a[i];
  return (word == 0);
}


/*------check whether a multiprecision integer is 1------*/
int int_is1(const Word *a, int len)
{
  Word word = 0;
  int i;
  
  for (i = len - 1; i > 0; i--) word |= a[i];
  return ((word == 0) & (a[0] == 1));
}


/*------set a multiprecision integer to a word------*/
void int_set(Word *r, Word a, int len)
{
  int i;
  
  for (i = len - 1; i > 0; i--) r[i] = 0;
  r[0] = a;
}


/*------comparison of multi-precision integers------*/
int int_cmp(const Word *a, const Word *b, int len)
{
  int i;
  
  for (i = len - 1; i >= 0; i--) {
    // printf("a[i] = %04x, b[i] = %04x\n", a[i], b[i]);
    if (a[i] > b[i]) return 1;   // a > b
    if (a[i] < b[i]) return -1;  // a < b
  }
  return 0;  // a = b
}


/*------copy a multi-precision integer------*/
void int_copy(Word *r, const Word *a, int len)
{
  int i;
  
  for (i = len - 1; i >= 0; i--) r[i] = a[i];
}


/*------printing of a multi-precision integer------*/
void int_print(const char *c, const Word *a, int len)
{
  int i;
  
  if ((c != NULL) && (strlen(c) > 0)) printf("%s", c);
  for (i = len - 1; i >= 0; i--) printf("%04x", a[i]);
  printf("\n");
}


/*------multiprecision right-shift------*/
int int_shr_c99(Word *r, const Word *a, int len)
{
  int i, retval;
  
  retval = a[0] & 1;  // return value
  for (i = 0; i < len - 1; i++) r[i] = (a[i+1] << (WSIZE - 1)) | (a[i] >> 1);
  r[len-1] = a[len-1] >> 1;
  
  return retval;
}


/*------multiprecision addition------*/
int int_add_c99(Word *r, const Word *a, const Word *b, int len)
{
  DWord sum;
  int i;
  
  sum = 0;
  for (i = 0; i < len; i++) {
    sum += (DWord) a[i] + b[i];
    r[i] = (Word) sum;
    sum >>= WSIZE;
  }
  
  return (int) sum;  // carry bit
}


/*------Multiprecision subtraction------*/
int int_sub_c99(Word *r, const Word *a, const Word *b, int len)
{
  DWord dif;
  int i;
  
  dif = 1;
  for (i = 0; i < len; i++) {
    dif += (DWord) a[i] + (~b[i]);
    r[i] = (Word) dif;
    dif >>= WSIZE;
  }
  
  return (1 - ((int) dif));  // borrow bit
}


/*------Multiprecision multiplication------*/
void int_mul_c99(Word *r, const Word *a, const Word *b, int len)
{
  DWord prod = 0;
  int i, j;
  
  // multiplication of A by b[0]
  for (j = 0; j < len; j++) {
    prod += (DWord) a[j]*b[0];
    r[j] = (Word) prod;
    prod >>= WSIZE;
  }
  r[j] = (Word) prod;
  
  // multiplication of A by b[i] for 1 <= i < len
  for (i = 1; i < len; i++) {
    prod = 0;
    for (j = 0; j < len; j++) {
      prod += (DWord) a[j]*b[i];
      prod += r[i+j];
      r[i+j] = (Word) prod;
      prod >>= WSIZE;
    }
    r[i+j] = (Word) prod;
  }
}


/*------Multiplication by 32-bit integer------*/
void int_mul32_c99(Word *r, const Word *a, const Word *b, int len)
{
  DWord prod = 0;
  int i, j;
  
  // multiplication of A by b[0]
  for (j = 0; j < len; j++) {
    prod += (DWord) a[j]*b[0];
    r[j] = (Word) prod;
    prod >>= WSIZE;
  }
  r[j] = (Word) prod;
  
  // multiplication of A by b[i] for 1 <= i < 32/WSIZE
  for (i = 1; i < 32/WSIZE; i++) {
    prod = 0;
    for (j = 0; j < len; j++) {
      prod += (DWord) a[j]*b[i];
      prod += r[i+j];
      r[i+j] = (Word) prod;
      prod >>= WSIZE;
    }
    r[i+j] = (Word) prod;
  }
}


/*------Multiprecision squaring------*/
void int_sqr_c99(Word *r, const Word *a, int len)
{
  DWord prod = 0, sum = 0;
  int i, j;
  
  // compute A[1,...,len-1]*a[0]
  r[0] = 0;
  for (j = 1; j < len; j++) {
    prod += (DWord) a[j]*a[0];
    r[j] = (Word) prod;
    prod >>= WSIZE;
  }
  r[j] = (Word) prod;
  
  // compute A[i+1,...,len-1]*a[i] for 1 <= i < len
  for (i = 1; i < len; i++) {
    prod = 0;
    for (j = i + 1; j < len; j++) {
      prod += (DWord) a[j]*a[i];
      prod += r[i+j];
      r[i+j] = (Word) prod;
      prod >>= WSIZE;
    }
    r[i+j] = (Word) prod;
  }
  
  // double the result obtained so far and
  // add the partial products a[i]*a[i]
  for (i = 0; i < len; i++) {
    prod = (DWord) a[i]*a[i];
    sum += (Word) prod;
    sum += (DWord) r[2*i] + r[2*i];
    r[2*i] = (Word) sum;
    sum >>= WSIZE;
    sum += ((Word) (prod >> WSIZE));
    sum += (DWord) r[2*i+1] + r[2*i+1];
    r[2*i+1] = (Word) sum;
    sum >>= WSIZE;
  }
}
