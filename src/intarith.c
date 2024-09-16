#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "config.h"
#include "intarith.h"


/*------check whether a multiprecision integer is 0------*/
int int_is0(const UINT16 *a, int len)
{
  int i;
  UINT16 word = 0;
  
  for (i = len-1; i >= 0; i--) word |= a[i];
  return (word == 0);
}


/*------check whether a multiprecision integer is 1------*/
int int_is1(const UINT16 *a, int len)
{
  int i;
  UINT16 word = 0;
  
  for (i = len-1; i > 0; i--) word |= a[i];
  return ((word == 0) & (a[0] == 1));
}


/*------set a multiprecision integer to a word------*/
void int_set(UINT16 *r, UINT16 a, int len)
{
  int i;
  
  for (i = len-1; i > 0; i--) r[i] = 0;
  r[0] = a;
}


/*------comparison of multi-precision integers------*/
int int_cmp(const UINT16 *a, const UINT16 *b, int len)
{
  int i;
  
  for (i = len-1; i >= 0; i--)
  {
    // printf("a[i] = %04x, b[i] = %04x\n", a[i], b[i]);
    if (a[i] > b[i]) return 1;   // a > b	
    if (a[i] < b[i]) return -1;  // a < b				
  }
  return 0;  // a = b
}


/*------copy a multi-precision integer------*/
void int_copy(UINT16 *r, const UINT16 *a, int len)
{
  int i;
  
  for (i = len-1; i >= 0; i--) r[i] = a[i];
}


/*------printing of a multi-precision integer------*/
void int_print(const char *c, const UINT16 *a, int len)
{
  int i;
  
  if ((c != NULL) && (strlen(c) > 0)) printf("%s", c);
  for (i = len - 1; i >= 0; i--) printf("%04x", a[i]);
  printf("\n");
}


/*------multiprecision right-shift------*/
int int_shr_c99(UINT16 *r, const UINT16 *a, int len)
{
  int i, retval;
  
  retval = a[0] & 1;  // return value
  for(i = 0; i < len-1; i++) r[i] = (a[i+1] << 15) | (a[i] >> 1);
  r[len-1] = a[len-1] >> 1;
  
  return retval;
}


/*------multiprecision addition------*/
int int_add_c99(UINT16 *r, const UINT16 *a, const UINT16 *b, int len)
{
  int i;
  UINT32 sum;
  
  sum = 0;
  for (i = 0; i < len; i++)
  {
    sum += (UINT32) a[i] + b[i];
    r[i] = (UINT16) sum;
    sum >>= 16;
  }
  
  return (int) sum;  // carry bit
}


/*------Multiprecision subtraction------*/
int int_sub_c99(UINT16 *r, const UINT16 *a, const UINT16 *b, int len)
{
  int i;
  UINT32 dif;
  
  dif = 1;
  for (i = 0; i < len; i++)
  {
    dif += (UINT32) a[i] + (~b[i]);
    r[i] = (UINT16) dif;
    dif >>= 16;
  }
  
  return (1 - ((int) dif));  // borrow bit
}


/*------Multiprecision multiplication------*/
void int_mul_c99(UINT16 *r, const UINT16 *a, const UINT16 *b, int len)
{
  int i, j;
  UINT32 prod = 0;
  
  /* multiplication of A by b[0] */
  for(j = 0; j < len; j ++) 
  {
    prod += (UINT32) a[j]*b[0];
    r[j] = (UINT16) prod;
    prod >>= 16;
  }
  r[j] = (UINT16) prod;
  
  /* multiplication of A by b[i] for 1 <= i < len */
  for(i = 1; i < len; i ++) 
  {
    prod = 0;
    for(j = 0; j < len; j ++) 
    {
      prod += (UINT32) a[j]*b[i];
      prod += r[i+j];
      r[i+j] = (UINT16) prod;
      prod >>= 16;
    }
    r[i+j] = (UINT16) prod;
  }
}


/*------Multiplication by 32-bit integer------*/
void int_mul32_c99(UINT16 *r, const UINT16 *a, const UINT16 *b, int len)
{
  int i;
  UINT32 prod;
  
  /* multiplication of A by b[0] */
  prod = 0;
  for (i = 0; i < len; i++)  
  { 
    prod += (UINT32) a[i]*b[0];
    r[i] = (UINT16) prod;
    prod >>= 16;
  }
  r[len] = (UINT16) prod;
  
  /* multiplication of A by b[1] */
  prod = 0;
  for (i = 0; i < len; i++)  
  { 
    prod += (UINT32) a[i]*b[1] + r[i+1];
    r[i+1] = (UINT16) prod;
    prod >>= 16;
  }
  r[len+1] = (UINT16) prod;
}


/*------Multiprecision squaring------*/
void int_sqr_c99(UINT16 *r, const UINT16 *a, int len)
{
  int i, j;
  UINT32 prod = 0, sum = 0;
  
  /* compute A[1,...,len-1]*a[0] */
  r[0] = 0;
  for(j = 1; j < len; j ++) 
  {
    prod += (UINT32) a[j]*a[0];
    r[j] = (UINT16) prod;
    prod >>= 16;
  }
  r[j] = (UINT16) prod;
  
  /* compute A[i+1,...,len-1]*a[i] for 1 <= i < len */
  for(i = 1; i < len; i ++) 
  {
    prod = 0;
    for(j = i+1; j < len; j ++) 
    {
      prod += (UINT32) a[j]*a[i];
      prod += r[i+j];
      r[i+j] = (UINT16) prod;
      prod >>= 16;
    }
    r[i+j] = (UINT16) prod;
  }
  
  /* double the result obtained so far and */
  /* add the partial products a[i]*a[i]    */
  for (i = 0; i < len; i ++)
  {
    prod = (UINT32) a[i]*a[i];
    sum += (UINT16) prod;
    sum += (UINT32) r[2*i] + r[2*i];
    r[2*i] = (UINT16) sum;
    sum >>= 16;
    sum += ((UINT16) (prod>>16));
    sum += (UINT32) r[2*i+1] + r[2*i+1];
    r[2*i+1] = (UINT16) sum;
    sum >>= 16;
  }
}
