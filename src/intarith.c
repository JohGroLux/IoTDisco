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
  int i, j, k;
  UINT64 accu;
  
  /* first iteration of 1st outer loop */
  accu = (UINT32) a[0]*b[0];
  r[0] = (UINT16) accu;
  accu >>= 16;
  /* len-1 iterations of 1st outer loop */
  for (i = 1; i < len; i++)  
  { 
    j = 0; k = i;
    while (k >= 0) accu += (UINT32) a[j++]*b[k--];
    r[i] = (UINT16) accu;
    accu >>= 16;
  }
  /* len-2 iterations of 2nd outer loop */
  for (i = len; i < 2*len-2; i++)
  { 
    k = len-1; j = i-k;
    while (j < len) accu += (UINT32) a[j++]*b[k--];
    r[i] = (UINT16) accu;
    accu >>= 16;
  }
  /* last iteration of 2nd outer loop */
  accu += (UINT32) a[len-1]*b[len-1];
  r[2*len-2] = (UINT16) accu;
  r[2*len-1] = (UINT16) (accu >> 16);
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
  int i, j, k;
  UINT64 accu = 0;
  
  r[0] = 0;
  /* len-1 iterations of 1st outer loop */  
  for(i = 1; i < len; i++) 
  {
    j = 0; k = i;
    while (j < k) accu += (UINT32) a[j++]*a[k--];
    r[i] = (UINT16) accu;
    accu >>= 16;
  }
  /* len-2 iterations of 2nd outer loop */
  for(i = len; i < 2*len-2; i++) 
  {
    k = len-1; j = i-k;
    while (j < k) accu += (UINT32) a[j++]*a[k--];
    r[i] = (UINT16) accu;
    accu >>= 16;
  }
  r[i] = (UINT16) accu;
  r[i+1] = 0;
  /* double the result obtained so far and add all a[i]*a[i] */
  accu = 0;
  for(i = 0; i < 2*len; i += 2) 
  {
    k = i>>1;
    accu += (UINT32) a[k]*a[k] + r[i] + r[i];
    r[i] = (UINT16) accu;
    accu >>= 16; 
    accu += (UINT32) r[i+1] + r[i+1];
    r[i+1] = (UINT16) accu;
    accu >>= 16;
  }
}
