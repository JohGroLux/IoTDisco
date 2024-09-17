#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "config.h"
#include "intarith.h"
#include "gfparith.h"


#ifdef MSPECC_USE_VLA  // use variable-length arrays for field elements
#define _len len       // requires "Allow VLA" in C/C++ Compiler Options
#else
#define _len MSPECC_MAX_LEN
#endif


// arithmetic right-shift of a signed 32-bit integer
#if ((((INT32) -1) >> 1) == ((INT32) -1))
#define ars32(x,y) ((x) >> (y))
#else
#define ars32(x,y) (((x) >> (y)) | ((~(((x) >> 31) - 1)) << (32 - (y))))
#endif


/*------Set r to p = 2^(16*len-1) - c------*/
void gfp_set(UINT16 *r, UINT16 c, int len)
{
  int i;
  
  r[len-1] = 0x7FFF;
  for (i = len-2; i > 0; i--) r[i] = 0xFFFF;
  r[0] = -c;
}


/*------compare a gfp element with p------*/
int gfp_isp(const UINT16 *a, UINT16 c, int len)
{
  int i;
  
  if (a[len-1] != 0x7FFF) return 0;  // a != p
  for (i = len-2; i > 0; i--) { if (a[i] != 0xFFFF) return 0; }  // a != p
  if (a[0] != -c) return 0;  // a != p
  
  return 1;  // a == p
}


/*------Modular addition------*/
void gfp_add_c99(UINT16 *r, const UINT16 *a, const UINT16 *b, UINT16 c, int len)
{
  int i;
  UINT16 msw;
  UINT32 sum;
  
  sum = (UINT32) a[len-1] + b[len-1];
  msw = ((UINT16) sum) & 0x7FFF;
  sum = (UINT32) c*((UINT16) (sum >> 15));
  
  for (i = 0; i < len-1; i++)
  {
    sum += (UINT32) a[i] + b[i];
    r[i] = (UINT16) sum;
    sum >>= 16;  // sum is <= 2
  }
  r[len-1] = msw + ((UINT16) sum);
}


/*------Modular subtraction------*/
void gfp_sub_c99(UINT16 *r, const UINT16 *a, const UINT16 *b, UINT16 c, int len)
{
  int i;
  UINT16 msw;
  INT32 sum;
  
  /* we compute r = 4*p + a - b mod p = 2^(16*len+1) + a - b - 4*c mod p */
  sum = (INT32) 0x1FFFC + a[len-1] - b[len-1];
  msw = ((UINT16) sum) & 0x7FFF;
  sum = (INT32) c*((UINT16) (sum >> 15));
  sum = sum - (c<<1) - (c<<1);  // (c<<1) can be up to 16 bits long!
  
  for (i = 0; i < len-1; i++)
  {
    sum += (INT32) a[i] - b[i];
    r[i] = (UINT16) sum;
    sum = ars32(sum, 16);  // arithmetic right-shift; now sum is in [-2,1]
  }
  r[len-1] = msw + ((UINT16) sum) + 4;  // 0x1FFFC+4 = 0x20000 = 4*p>>16*(len-1)
}


/*------Modular subtraction (2nd variant)------*/
void gfp_sub_c99_v2(UINT16 *r, const UINT16 *a, const UINT16 *b, UINT16 c, int len)
{
  int i;
  UINT16 msw;
  UINT32 sum;
  
  /* we compute r = 4*p + a - b */
  sum = (UINT32) 0x1FFFC + a[len-1] - b[len-1];
  msw = ((UINT16) sum) & 0x7FFF;
  sum = (UINT32) c*((UINT16) (sum >> 15));
  sum = sum - (c<<1) - (c<<1) + 4;  // (c<<1) can be up to 16 bits long!
  
  for (i = 0; i < len-1; i++)
  {
    sum += (UINT32) 0x3FFFC + a[i] - b[i];
    r[i] = (UINT16) sum;
    sum >>= 16;
  }
  r[len-1] = msw + ((UINT16) sum);
}


/*------Conditional negation------*/
void gfp_cneg_c99(UINT16 *r, const UINT16 *a, UINT16 c, int neg, int len)
{
  int i;
  UINT16 msw, mask;
  INT32 sum;
  
  mask = ~(((UINT16) (neg&1)) - 1);
  
  // t = 2^(k+1) - a - 1 is the (k+1)-bit one's complement of a. When the LSB
  // of neg = 1 then we compute r = 4*p - a mod p = 2*(2^(k+1) - 2*c) - a mod p
  // = 2^(k+1) + 2^(k+1) - 4*c - a - 1 + 1 mod p = 2^{k+1) + t - 4*c + 1 mod p.
  // The one's complement t can be obtained by xoring each word of a with an
  // all-1 mask. On the other hand, when the LSB of neg = 0, we simply compute
  // r = a mod p.
  
  sum = (INT32) (mask & 0xFFFC) + (mask ^ a[len-1]);
  msw = ((UINT16) sum) & 0x7FFF;
  sum = (INT32) c*((UINT16) (sum >> 15));
  sum = sum - (mask & (c<<1)) - (mask & (c<<1)) + (mask & 1);
  
  for (i = 0; i < len-1; i++)
  {
    sum += (INT32) (mask ^ a[i]);
    r[i] = (UINT16) sum;
    sum = ars32(sum, 16);  // arithmetic right-shift; now sum is in [-1,1]
  }
  r[len-1] = msw + ((UINT16) sum) + (mask & 4);
}


/*------Modular halving------*/
void gfp_hlv_c99(UINT16 *r, const UINT16 *a, UINT16 c, int len)
{
  int i;
  UINT16 tmp, mask;
  INT32 sum;
  
  // masked addition of prime p to a
  mask = ~((a[0] & 1) - 1);  // 0 or 0xFFFF
  sum = (INT32) a[0] - (c & mask);
  tmp = (UINT16) sum;
  sum = ars32(sum, 16);  // arithmetic right-shift; now sum is in [-1,0]
  
  for (i = 1; i < len-1; i ++)
  {
    sum += (INT32) a[i];
    r[i-1] = (((UINT16) sum) << 15) | (tmp >> 1);
    tmp = (UINT16) sum;
    sum = ars32(sum, 16);  // arithmetic right-shift; now sum is in [-1,0]
  }
  sum += (INT32) a[len-1] + (0x8000 & mask);
  r[len-2] = (((UINT16) sum) << 15) | (tmp >> 1);
  r[len-1] = (UINT16) (sum >> 1);
}


/*------Modular halving (2nd variant)------*/
void gfp_hlv_c99_v2(UINT16 *r, const UINT16 *a, UINT16 c, int len)
{
  int i;
  UINT16 tmp, mask;
  UINT32 sum;
  
  // masked addition of prime p to a
  mask = ~((a[0] & 1) - 1);  // 0 or 0xFFFF
  // mask = 0 - (a[0] & 1);  // 0 or 0xFFFF
  sum = (UINT32) a[0] + ((0-c) & mask);
  tmp = (UINT16) sum;
  sum >>= 16;
  for (i = 1; i < len-1; i ++)
  {
    sum += (UINT32) a[i] + mask;
    r[i-1] = (((UINT16) sum) << 15) | (tmp >> 1);
    tmp = (UINT16) sum;
    sum >>= 16;
  }
  sum += (UINT32) a[len-1] + (0x7FFF & mask);
  r[len-2] = (((UINT16) sum) << 15) | (tmp >> 1);
  r[len-1] = (UINT16) (sum >> 1);
}


/*------Reduction by a pseudo-Mersenne prime p=2^n-c------*/
void gfp_red_c99(UINT16 *r, const UINT16 *a, UINT16 c, int len)
{
  int i;
  UINT16 msw, d = (c<<1);
  UINT32 prod, sum;
  
  /* first round */
  prod = 0;
  for (i = 0; i < len-1; i++)
  {
    prod += (UINT32) a[i+len]*d + a[i];
    r[i] = (UINT16) prod;
    prod >>= 16;
  }
  prod += (UINT32) a[2*len-1]*d + a[len-1];
  
  /* second round */
  msw = ((UINT16) prod) & 0x7FFF;
  sum = (UINT32) c*(prod >>= 15);  // sum is max 32 bits if c is max 15 bits
  for (i = 0; i < len-1; i++)
  {
    sum += r[i];
    r[i] = (UINT16) sum;
    sum >>= 16; 
  }
  r[len-1] = msw + ((UINT16) sum);
}


/*------Modular multiplication------*/
void gfp_mul_c99(UINT16 *r, const UINT16 *a, const UINT16 *b, UINT16 c, int len)
{
  int i, j;
  UINT16 msw, d = (c<<1);
  UINT32 prod = 0;
  UINT16 t[2*MSPECC_MAX_LEN];
  
  /* multiplication of A by b[0] */
  for(j = 0; j < len; j ++) 
  {
    prod += (UINT32) a[j]*b[0];
    t[j] = (UINT16) prod;
    prod >>= 16;
  }
  t[j] = (UINT16) prod;
  
 /* multiplication of A by b[i] for 1 <= i < len */
  for(i = 1; i < len; i ++) 
  {
    prod = 0;
    for(j = 0; j < len; j ++) 
    {
      prod += (UINT32) a[j]*b[i];
      prod += t[i+j];
      t[i+j] = (UINT16) prod;
      prod >>= 16;
    }
    t[i+j] = (UINT16) prod;
  }
  
  /* first round of modular reduction */
  prod = 0;
  for (i = 0; i < len-1; i ++)
  {
    prod += (UINT32) t[i+len]*d + t[i];
    t[i] = (UINT16) prod;
    prod >>= 16;
  }
  prod += (UINT32) t[2*len-1]*d + t[len-1];
  
  /* second round of modular reduction */
  msw = ((UINT16) prod) & 0x7FFF;
  prod = (UINT32) c*(prod >> 15);  // prod is max 32 bits if c is max 15 bits
  for (i = 0; i < len-1; i ++)
  {
    prod += t[i];
    r[i] = (UINT16) prod;
    prod >>= 16; 
  }
  r[len-1] = msw + ((UINT16) prod);
}


/*------Modular squaring------*/
void gfp_sqr_c99(UINT16 *r, const UINT16 *a, UINT16 c, int len)
{
  int i, j;
  UINT16 msw, d = (c<<1);
  UINT32 prod = 0, sum = 0;
  UINT16 t[2*MSPECC_MAX_LEN];
  
  /* compute A[1,...,len-1]*a[0] */
  t[0] = 0;
  for(j = 1; j < len; j ++) 
  {
    prod += (UINT32) a[j]*a[0];
    t[j] = (UINT16) prod;
    prod >>= 16;
  }
  t[j] = (UINT16) prod;
  
  /* compute A[i+1,...,len-1]*a[i] for 1 <= i < len */
  for(i = 1; i < len; i ++) 
  {
    prod = 0;
    for(j = i+1; j < len; j ++) 
    {
      prod += (UINT32) a[j]*a[i];
      prod += t[i+j];
      t[i+j] = (UINT16) prod;
      prod >>= 16;
    }
    t[i+j] = (UINT16) prod;
  }
  
  /* double the result obtained so far and */
  /* add the partial products a[i]*a[i]    */
  for (i = 0; i < len; i ++)
  {
    prod = (UINT32) a[i]*a[i];
    sum += (UINT16) prod;
    sum += (UINT32) t[2*i] + t[2*i];
    t[2*i] = (UINT16) sum;
    sum >>= 16;
    sum += ((UINT16) (prod>>16));
    sum += (UINT32) t[2*i+1] + t[2*i+1];
    t[2*i+1] = (UINT16) sum;
    sum >>= 16;
  }
  
  /* first round of modular reduction */
  prod = 0;
  for (i = 0; i < len-1; i ++)
  {
    prod += (UINT32) t[i+len]*d + t[i];
    t[i] = (UINT16) prod;
    prod >>= 16;
  }
  prod += (UINT32) t[2*len-1]*d + t[len-1];
  
  /* second round of modular reduction */
  msw = ((UINT16) prod) & 0x7FFF;
  prod = (UINT32) c*(prod >> 15);  // prod is max 32 bits if c is max 15 bits
  for (i = 0; i < len-1; i ++)
  {
    prod += t[i];
    r[i] = (UINT16) prod;
    prod >>= 16; 
  }
  r[len-1] = msw + ((UINT16) prod);
}


/*------Reduction of a (16*len+32)-bit integer------*/
void gfp_red32_c99(UINT16 *r, const UINT16 *a, UINT16 c, int len)
{
  int i;
  UINT16 msw, word;
  UINT32 prod;
  
  prod = 0;
  msw = a[len-1] & 0x7FFF;
  
  // compute words r[0] and r[1]
  for (i = 0; i < 2; i ++)
  {
    word = (a[i+len] << 1) | (a[i+len-1] >> 15);
    prod += (UINT32) word*c + a[i];
    r[i] = (UINT16) prod;
    prod >>= 16;
  }
  
  // compute word r[2]
  word = -(a[len+1] >> 15);  // either 0 or 0xFFFF
  prod += (UINT32) (word&c) + a[2];
  r[2] = (UINT16) prod;
  prod >>= 16;
  
  // compute r[i] = a[i] + carry
  for (i = 3; i < len-1; i++)
  {
    prod += (UINT32) a[i];
    r[i] = (UINT16) prod;
    prod >>= 16;
  }
  r[len-1] = ((UINT16) prod) + msw;
}


/*------Modular Multiplication by 32-bit integer------*/
void gfp_mul32_c99(UINT16 *r, const UINT16 *a, const UINT16 *b, UINT16 c, int len)
{
  int i;
  UINT16 msw, word;
  UINT32 prod = 0;
  UINT16 t[MSPECC_MAX_LEN+2];
  
  /* multiplication of A by b[0] */
  for (i = 0; i < len; i++)  
  { 
    prod += (UINT32) a[i]*b[0];
    t[i] = (UINT16) prod;
    prod >>= 16;
  }
  t[len] = (UINT16) prod;
  
  /* multiplication of A by b[1] */
  prod = 0;
  for (i = 0; i < len; i++)  
  { 
    prod += (UINT32) a[i]*b[1] + t[i+1];
    t[i+1] = (UINT16) prod;
    prod >>= 16;
  }
  t[len+1] = (UINT16) prod;
  
  prod = 0;
  msw = t[len-1] & 0x7FFF;
  
  // compute words r[0] and r[1]
  for (i = 0; i < 2; i ++)
  {
    word = (t[i+len] << 1) | (t[i+len-1] >> 15);
    prod += (UINT32) word*c + t[i];
    r[i] = (UINT16) prod;
    prod >>= 16;
  }
  
  // compute word r[2]
  word = -(t[len+1] >> 15);  // either 0 or 0xFFFF
  prod += (UINT32) (word&c) + t[2];
  r[2] = (UINT16) prod;
  prod >>= 16;
  
  // compute r[i] = a[i] + carry
  for (i = 3; i < len-1; i++)
  {
    prod += (UINT32) t[i];
    r[i] = (UINT16) prod;
    prod >>= 16;
  }
  r[len-1] = ((UINT16) prod) + msw;
}


/*------Least non-negative residue------*/
void gfp_lnr(UINT16 *r, const UINT16 *a, UINT16 c, int len)
{
  int i;
  UINT16 mask;
  UINT32 sum = c;
  
  // we first compute r = a - p by adding the two's complement of p
  for (i = 0; i < len-1; i++)
  {
    sum += (UINT32) a[i];
    r[i] = (UINT16) sum;
    sum >>= 16;
  }
  sum += (UINT32) a[len-1] + 0x8000;
  r[len-1] = (UINT16) sum;
  
  // mask is 0 when the addition produced a carry, 0xFFFF otherwise
  mask = ((UINT16) (sum>>16)) - 1;
  
  // perform masked addition of p (i.e. if r < 0 compute r = r + p)
  sum = (UINT32) r[0] + ((-c) & mask);
  r[0] = (UINT16) sum;
  sum >>= 16;
  for (i = 1; i < len-1; i++)
  {
    sum += (UINT32) r[i] + mask;
    r[i] = (UINT16) sum;
    sum >>= 16;
  }
  sum += (UINT32) r[len-1] + (0x7FFF & mask);
  r[len-1] = (UINT16) sum;
}


/*------compare two gfp elements that may be incompletely reduced------*/
int gfp_cmp(UINT16 *a, UINT16 *b, UINT16 c, int len)
{
  int i;
  UINT16 diff = 0;
  
  gfp_lnr(a, a, c, len);
  gfp_lnr(b, b, c, len);
  // compare a and b in constant time
  for (i = len-1; i >= 0; i--) diff |= (a[i] ^ b[i]);
  
  return (diff != 0);
}


/*****************************************************************************/
/* Inversion r = a^(-1) mod m for a special 'low-weight' modulus of the form */
/* m = u*2^k + 1 such as used for OPFs. The operand 'u' specifies the 16     */
/* most-significant bits of m and 'n' refers to the number of 32-bit words   */
/* the operand 'a' consists of. The result 'r' may not be fully reduced, but */
/* it is at most 'n' words long.                                             */
/*****************************************************************************/

int gfp_inv(UINT16 *r, const UINT16 *a, UINT16 c, int len)
{
  int uvlen = len;
  UINT16 tmp[3*_len];  // temporary space for three gfp elements
  UINT16 *ux = tmp, *vx = &tmp[_len], *x1 = &tmp[2*_len], *x2 = r;
  
  int_copy(ux, a, len);  // set ux = a
  gfp_set(vx, c, len);   // set vx = p
  int_set(x1, 1, len);   // set x1 = 1
  int_set(x2, 0, len);   // set x2 = 0
  
  while (int_cmp(ux, vx, len) >= 0) int_sub(ux, ux, vx, len);
  if (int_is0(ux, len)) return MSPECC_ERR_INVERSION_ZERO;
  
  while((!int_is1(ux, uvlen)) && (!int_is1(vx, uvlen)))
  {
    while((ux[0]&1) == 0)   // ux is even
    {
      int_shr(ux, ux, uvlen);
      gfp_hlv(x1, x1, c, len);
    }
    while((vx[0]&1) == 0)   // vx is even
    {
      int_shr(vx, vx, uvlen);
      gfp_hlv(x2, x2, c, len);
    }
    if (int_cmp(ux, vx, uvlen) >= 0)
    {
      int_sub(ux, ux, vx, uvlen);
      gfp_sub(x1, x1, x2, c, len);
    }
    else
    {
      int_sub(vx, vx, ux, uvlen);
      gfp_sub(x2, x2, x1, c, len);
    }
    if ((ux[uvlen-1] == 0) && (vx[uvlen-1] == 0)) uvlen--;
  }
  
  if (int_is1(ux, len)) int_copy(r, x1, len);
  return MSPECC_NO_ERROR;
}
