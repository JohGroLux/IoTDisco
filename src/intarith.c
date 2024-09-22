#include <assert.h>  // assert()
#include <stdint.h>  // uint32_t
#include <stdio.h>   // printf()
#include <stdlib.h>  // ??
#include <string.h>  // strlen()
#include "intarith.h"


// printf() formatstring
#if (WSIZE == 16)
#define FORMATSTR "%04x"
#else  // WSIZE == 32
#define FORMATSTR "%08x"
#endif


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
  Word agb = 0, asb = 0;
  int i, r = 0;
  
  // assert(WSIZE >= len);
  
  for (i = len - 1; i >= 0; i--) {
    agb = (agb << 1) | (a[i] > b[i]);
    asb = (asb << 1) | (a[i] < b[i]);
  }
  r += (agb > asb);  // r = +1 if a is greater than b
  r -= (agb < asb);  // r = -1 of a is smaller than b
  
  return r;  // r = 0 if a equals b
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
  for (i = len - 1; i >= 0; i--) printf(FORMATSTR, a[i]);
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
  DWord sum = 0;
  int i;
  
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
  DWord dif = 1;
  int i;
  
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
  
  // assert((r != a) && (r != b));
  
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
  
  // assert(r != a);
  
  // multiplication of A[1,...,len-1] by a[0]
  r[0] = 0;
  for (j = 1; j < len; j++) {
    prod += (DWord) a[j]*a[0];
    r[j] = (Word) prod;
    prod >>= WSIZE;
  }
  r[j] = (Word) prod;
  
  // multiplication of A[i+1,...,len-1] by a[i] for 1 <= i < len
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
  
  // double existing result, add squares a[i]^2 for 0 <= i < len
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
