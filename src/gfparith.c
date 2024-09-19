#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "config.h"
#include "intarith.h"
#include "gfparith.h"


typedef uint16_t Word;   // single-length word
typedef uint32_t DWord;  // double-length word
typedef int32_t SDWord;  // signed double-length word


// Bitlength of a Word
#define WSIZE (8*sizeof(Word))
// All-1 Mask: 0xFF..FF
#define ALL1MASK ((Word) -1L)
// MSB-1 Mask: 0x80..00
#define MSB1MASK ((Word) (1L << (WSIZE - 1)))
// MSB-0 Mask: 0x7F..FF
#define MSB0MASK (ALL1MASK >> 1)
// Minus-4 Mask: 0xFF..FC
#define MIN4MASK  ((Word) -4L)
// 4*p[len-1] (WSIZE+1 bits long)
#define FOURXPHI (((DWord) MSB0MASK) << 2)
// 4*p[i], 1 <= i < len-1 (WSIZE+2 bits long)
#define FOURXPMI (((DWord) ALL1MASK) << 2)


#ifdef MSPECC_USE_VLA  // use variable-length arrays for field elements
#define _len len     // requires "Allow VLA" in C/C++ Compiler Options
#else
#define _len MSPECC_MAX_LEN
#endif


// Arithmetic Right-Shift (ARS) of a SDWord
#if ((((SDWord) -1) >> 1) == ((SDWord) -1))
#define SDWARS(x,y) ((x) >> (y))
#else
// SDMASK is all-1 if x is negative, 0 otherwise
#define SDMASK(x) (~(((x) >> (2*WSIZE - 1)) - 1))
#define SDWARS(x,y) (((x) >> (y)) | (SDMASK(x) << (2*WSIZE - (y))))
#endif


/*------Set r to p = 2^(w*len-1) - c------*/
void gfp_set(Word *r, Word c, int len)
{
  int i;
  
  r[len-1] = MSB0MASK;  // 0x7F..FF
  for (i = len - 2; i > 0; i--) r[i] = ALL1MASK;  // 0xFF..FF
  r[0] = 0 - c;
}


/*------Compare a gfp element with p------*/
int gfp_isp(const Word *a, Word c, int len)
{
  int i;
  
  if (a[len-1] != MSB0MASK) return 0;  // 0x7F..FF
  for (i = len - 2; i > 0; i--) {
    if (a[i] != ALL1MASK) return 0;    // 0xFF..FF
  }
  if (a[0] != (0 - c)) return 0;       // a != p
  
  return 1;  // a == p
}


/*------Modular addition------*/
void gfp_add_c99(Word *r, const Word *a, const Word *b, Word c, int len)
{
  DWord sum;
  Word msw;
  int i;
  
  sum = (DWord) a[len-1] + b[len-1];
  msw = ((Word) sum) & MSB0MASK;  // 0x7F..FF
  sum = (DWord) c*((Word) (sum >> (WSIZE - 1)));
  
  for (i = 0; i < len - 1; i++) {
    sum += (DWord) a[i] + b[i];
    r[i] = (Word) sum;
    sum >>= WSIZE;  // sum is <= 2
  }
  r[len-1] = msw + ((Word) sum);
}


/*------Modular subtraction------*/
void gfp_sub_c99(Word *r, const Word *a, const Word *b, Word c, int len)
{
  SDWord sum;  // signed!
  Word msw;
  int i;
  
  // we compute r = 4*p + a - b mod p = 2^(k+2) + a - b - 4*c mod p
  sum = (SDWord) FOURXPHI + a[len-1] - b[len-1];  // 0x1FF..FC
  msw = ((Word) sum) & MSB0MASK;  // 0x7F..FF
  sum = (SDWord) c*((Word) (sum >> (WSIZE - 1)));
  sum = sum - (c << 1) - (c << 1);  // (c << 1) can be up to WSIZE bits long!
  
  for (i = 0; i < len - 1; i++) {
    sum += (SDWord) a[i] - b[i];
    r[i] = (Word) sum;
    sum = SDWARS(sum, WSIZE);  // arithmetic right-shift: sum is now in [-2,1]
  }
  r[len-1] = msw + ((Word) sum) + 4;  // 0x1FF..FC+4=0x200..00 = MSW of 2^(k+2)
}


/*------Modular subtraction (2nd variant)------*/
void gfp_sub_c99_v2(Word *r, const Word *a, const Word *b, Word c, int len)
{
  DWord sum;
  Word msw;
  int i;
  
  // we compute r = 4*p + a - b
  sum = (DWord) FOURXPHI + a[len-1] - b[len-1];  // 0x1FF..FC
  msw = ((Word) sum) & MSB0MASK;  // 0x7F..FF
  sum = (DWord) c*((Word) (sum >> (WSIZE - 1)));
  sum = sum - (c << 1) - (c << 1) + 4;  // (c<<1) can be up to WSIZE bits long!
  
  for (i = 0; i < len - 1; i++) {
    sum += (DWord) FOURXPMI + a[i] - b[i];  // 0x3FF..FC
    r[i] = (Word) sum;
    sum >>= WSIZE;
  }
  r[len-1] = msw + ((Word) sum);
}


/*------Conditional negation------*/
void gfp_cneg_c99(Word *r, const Word *a, Word c, int neg, int len)
{
  SDWord sum;  // signed!
  Word msw, mask;
  int i;
  
  mask = ~(((Word) (neg & 1)) - 1);  // 0 or all-1
  
  // t = 2^(k+1) - a - 1 is the (k+1)-bit one's complement of a. When the LSB
  // of neg = 1 then we compute r = 4*p - a mod p = 2*(2^(k+1) - 2*c) - a mod p
  // = 2^(k+1) + 2^(k+1) - 4*c - a - 1 + 1 mod p = 2^(k+1) + t - 4*c + 1 mod p.
  // The one's complement t can be obtained by xoring each word of a with an
  // all-1 mask. On the other hand, when the LSB of neg = 0, we simply compute
  // r = a mod p.
  
  sum = (SDWord) (mask & MIN4MASK) + (mask ^ a[len-1]);  // 0xFF..FC
  msw = ((Word) sum) & MSB0MASK;  // 0x7F..FF
  sum = (SDWord) c*((Word) (sum >> (WSIZE - 1)));
  sum = sum - (mask & (c << 1)) - (mask & (c << 1)) + (mask & 1);
  
  for (i = 0; i < len - 1; i++) {
    sum += (SDWord) (mask ^ a[i]);
    r[i] = (Word) sum;
    sum = SDWARS(sum, WSIZE);  // arithmetic right-shift: sum is now in [-1,1]
  }
  r[len-1] = msw + ((Word) sum) + (mask & 4);
}


/*------Modular halving------*/
void gfp_hlv_c99(Word *r, const Word *a, Word c, int len)
{
  SDWord sum;  // signed!
  Word tmp, mask;
  int i;
  
  // masked addition of prime p to a
  mask = ~((a[0] & 1) - 1);  // 0 or all-1
  sum = (SDWord) a[0] - (c & mask);
  tmp = (Word) sum;
  sum = SDWARS(sum, WSIZE);  // arithmetic right-shift: sum is now in [-1,0]
  
  for (i = 1; i < len - 1; i++) {
    sum += (SDWord) a[i];
    r[i-1] = (((Word) sum) << (WSIZE - 1)) | (tmp >> 1);
    tmp = (Word) sum;
    sum = SDWARS(sum, WSIZE);  // arithmetic right-shift: sum is now in [-1,0]
  }
  sum += (SDWord) a[len-1] + (MSB1MASK & mask);  // 0x80..00
  r[len-2] = (((Word) sum) << (WSIZE - 1)) | (tmp >> 1);
  r[len-1] = (Word) (sum >> 1);
}


/*------Modular halving (2nd variant)------*/
void gfp_hlv_c99_v2(Word *r, const Word *a, Word c, int len)
{
  DWord sum;
  Word tmp, mask;
  int i;
  
  // masked addition of prime p to a
  mask = ~((a[0] & 1) - 1);  // 0 or all-1
  // mask = 0 - (a[0] & 1);  // 0 or all-1
  sum = (DWord) a[0] + ((0 - c) & mask);
  tmp = (Word) sum;
  sum >>= WSIZE;
  
  for (i = 1; i < len - 1; i++) {
    sum += (DWord) a[i] + mask;
    r[i-1] = (((Word) sum) << (WSIZE - 1)) | (tmp >> 1);
    tmp = (Word) sum;
    sum >>= WSIZE;
  }
  sum += (DWord) a[len-1] + (mask >> 1);
  r[len-2] = (((Word) sum) << (WSIZE - 1)) | (tmp >> 1);
  r[len-1] = (Word) (sum >> 1);
}


/*------Reduction by a pseudo-Mersenne prime p=2^n-c------*/
void gfp_red_c99(Word *r, const Word *a, Word c, int len)
{
  DWord prod = 0, sum;
  Word msw, d = (c << 1);
  int i;
  
  // first step
  for (i = 0; i < len - 1; i++) {
    prod += (DWord) a[i+len]*d + a[i];
    r[i] = (Word) prod;
    prod >>= WSIZE;
  }
  prod += (DWord) a[2*len-1]*d + a[len-1];
  
  // second step
  msw = ((Word) prod) & MSB0MASK;  // 0x7F..FF
  sum = (DWord) c*(prod >> (WSIZE - 1));  // sum is max 2*WSIZE bits long!
  for (i = 0; i < len - 1; i++) {
    sum += r[i];
    r[i] = (Word) sum;
    sum >>= WSIZE;
  }
  r[len-1] = msw + ((Word) sum);
}


/*------Modular multiplication------*/
void gfp_mul_c99(Word *r, const Word *a, const Word *b, Word c, int len)
{
  Word t[2*_len];
  DWord prod = 0;
  Word msw, d = (c << 1);
  int i, j;
  
  // multiplication of A by b[0]
  for (j = 0; j < len; j++) {
    prod += (DWord) a[j]*b[0];
    t[j] = (Word) prod;
    prod >>= WSIZE;
  }
  t[j] = (Word) prod;
  
  // multiplication of A by b[i] for 1 <= i < len
  for (i = 1; i < len; i++) {
    prod = 0;
    for (j = 0; j < len; j++) {
      prod += (DWord) a[j]*b[i];
      prod += t[i+j];
      t[i+j] = (Word) prod;
      prod >>= WSIZE;
    }
    t[i+j] = (Word) prod;
  }
  
  // first step of modular reduction
  prod = 0;
  for (i = 0; i < len - 1; i++) {
    prod += (DWord) t[i+len]*d + t[i];
    t[i] = (Word) prod;
    prod >>= WSIZE;
  }
  prod += (DWord) t[2*len-1]*d + t[len-1];
  
  // second step of modular reduction
  msw = ((Word) prod) & MSB0MASK;  // 0x7F..FF
  prod = (DWord) c*(prod >> (WSIZE - 1));  // prod is max 2*WSIZE bits long!
  for (i = 0; i < len - 1; i++) {
    prod += t[i];
    r[i] = (Word) prod;
    prod >>= WSIZE;
  }
  r[len-1] = msw + ((Word) prod);
}


/*------Modular squaring------*/
void gfp_sqr_c99(Word *r, const Word *a, Word c, int len)
{
  Word t[2*_len];
  DWord prod = 0, sum = 0;
  Word msw, d = (c << 1);
  int i, j;
  
  // compute A[1,...,len-1]*a[0]
  t[0] = 0;
  for (j = 1; j < len; j++) {
    prod += (DWord) a[j]*a[0];
    t[j] = (Word) prod;
    prod >>= WSIZE;
  }
  t[j] = (Word) prod;
  
  // compute A[i+1,...,len-1]*a[i] for 1 <= i < len
  for (i = 1; i < len; i++) {
    prod = 0;
    for (j = i + 1; j < len; j++) {
      prod += (DWord) a[j]*a[i];
      prod += t[i+j];
      t[i+j] = (Word) prod;
      prod >>= WSIZE;
    }
    t[i+j] = (Word) prod;
  }
  
  // double the result obtained so far and
  // add the partial products a[i]*a[i]
  for (i = 0; i < len; i++) {
    prod = (DWord) a[i]*a[i];
    sum += (Word) prod;
    sum += (DWord) t[2*i] + t[2*i];
    t[2*i] = (Word) sum;
    sum >>= WSIZE;
    sum += ((Word) (prod >> WSIZE));
    sum += (DWord) t[2*i+1] + t[2*i+1];
    t[2*i+1] = (Word) sum;
    sum >>= WSIZE;
  }
  
  // first step of modular reduction
  prod = 0;
  for (i = 0; i < len - 1; i++) {
    prod += (DWord) t[i+len]*d + t[i];
    t[i] = (Word) prod;
    prod >>= WSIZE;
  }
  prod += (DWord) t[2*len-1]*d + t[len-1];
  
  // second step of modular reduction
  msw = ((Word) prod) & MSB0MASK;  // 0x7F..FF
  prod = (DWord) c*(prod >> (WSIZE - 1));  // prod is max 2*WSIZE bits long!
  for (i = 0; i < len - 1; i++) {
    prod += t[i];
    r[i] = (Word) prod;
    prod >>= WSIZE;
  }
  r[len-1] = msw + ((Word) prod);
}


/*------Reduction of a (w*len+32)-bit integer------*/
void gfp_red32_c99(Word *r, const Word *a, Word c, int len)
{
  DWord prod;
  Word msw, d = (c << 1);
  int i;
  
  prod = (DWord) c*(a[len-1] >> (WSIZE - 1));
  msw = a[len-1] & MSB0MASK;  // 0x7F..FF
  
  // compute first 32 bits of result
  for (i = 0; i < 32/WSIZE; i++) {
    prod += (DWord) a[i+len]*d + a[i];
    r[i] = (Word) prod;
    prod >>= WSIZE;
  }
    
  // compute r[i] = t[i] + carry
  for (i = 32/WSIZE; i < len - 1; i++) {
    prod += (DWord) a[i];
    r[i] = (Word) prod;
    prod >>= WSIZE;
  }
  r[len-1] = ((Word) prod) + msw;
}


/*------Modular Multiplication by 32-bit integer------*/
void gfp_mul32_c99(Word *r, const Word *a, const Word *b, Word c, int len)
{
  Word t[_len+(32/WSIZE)];
  DWord prod = 0;
  Word msw, d = (c << 1);
  int i, j;
  
  // multiplication of A by b[0]
  for (j = 0; j < len; j++) {
    prod += (DWord) a[j]*b[0];
    t[j] = (Word) prod;
    prod >>= WSIZE;
  }
  t[j] = (Word) prod;
  
  // multiplication of A by b[i] for 1 <= i < 32/WSIZE
  for (i = 1; i < 32/WSIZE; i++) {
    prod = 0;
    for (j = 0; j < len; j++) {
      prod += (DWord) a[j]*b[i];
      prod += t[i+j];
      t[i+j] = (Word) prod;
      prod >>= WSIZE;
    }
    t[i+j] = (Word) prod;
  }
  
  prod = (DWord) c*(t[len-1] >> (WSIZE - 1));
  msw = t[len-1] & MSB0MASK;  // 0x7F..FF
  
  // compute first 32 bits of result
  for (i = 0; i < 32/WSIZE; i++) {
    prod += (DWord) t[i+len]*d + t[i];
    r[i] = (Word) prod;
    prod >>= WSIZE;
  }
    
  // compute r[i] = t[i] + carry
  for (i = 32/WSIZE; i < len - 1; i++) {
    prod += (DWord) t[i];
    r[i] = (Word) prod;
    prod >>= WSIZE;
  }
  r[len-1] = ((Word) prod) + msw;
}


/*------Least non-negative residue------*/
void gfp_lnr(Word *r, const Word *a, Word c, int len)
{
  DWord sum = c;
  Word mask;
  int i;
  
  // we first compute r = a - p by adding the two's complement of p
  for (i = 0; i < len - 1; i++) {
    sum += (DWord) a[i];
    r[i] = (Word) sum;
    sum >>= WSIZE;
  }
  sum += (DWord) a[len-1] + MSB1MASK;  // 0x80..00
  r[len-1] = (Word) sum;
  
  // mask is 0 when the addition produced a carry, all-1 otherwise
  mask = ((Word) (sum >> WSIZE)) - 1;
  
  // perform masked addition of p (i.e. if r < 0 compute r = r + p)
  sum = (DWord) r[0] + ((0 - c) & mask);
  r[0] = (Word) sum;
  sum >>= WSIZE;
  for (i = 1; i < len - 1; i++) {
    sum += (DWord) r[i] + mask;
    r[i] = (Word) sum;
    sum >>= WSIZE;
  }
  sum += (DWord) r[len-1] + (mask >> 1);
  r[len-1] = (Word) sum;
}


/*------Compare two gfp elements that may be incompletely reduced------*/
int gfp_cmp(Word *a, Word *b, Word c, int len)
{
  Word diff = 0;
  int i;
  
  gfp_lnr(a, a, c, len);
  gfp_lnr(b, b, c, len);
  // compare a and b in constant time
  for (i = len - 1; i >= 0; i--) diff |= (a[i] ^ b[i]);
  
  return (diff != 0);
}


/*------Inversion r = a^(-1) mod m for a modulus of the form m = 2^(w*len-1) - c------*/
int gfp_inv(Word *r, const Word *a, Word c, int len)
{
  Word tmp[3*_len];  // temporary space for three gfp elements
  Word *ux = tmp, *vx = &tmp[_len], *x1 = &tmp[2*_len], *x2 = r;
  int uvlen = len;
  
  int_copy(ux, a, len);  // set ux = a
  gfp_set(vx, c, len);   // set vx = p
  int_set(x1, 1, len);   // set x1 = 1
  int_set(x2, 0, len);   // set x2 = 0
  
  while (int_cmp(ux, vx, len) >= 0) int_sub(ux, ux, vx, len);
  if (int_is0(ux, len)) return MSPECC_ERR_INVERSION_ZERO;
  
  while((!int_is1(ux, uvlen)) && (!int_is1(vx, uvlen))) {
    while((ux[0] & 1) == 0) {  // ux is even
      int_shr(ux, ux, uvlen);
      gfp_hlv(x1, x1, c, len);
    }
    while((vx[0] & 1) == 0) {  // vx is even
      int_shr(vx, vx, uvlen);
      gfp_hlv(x2, x2, c, len);
    }
    if (int_cmp(ux, vx, uvlen) >= 0) {
      int_sub(ux, ux, vx, uvlen);
      gfp_sub(x1, x1, x2, c, len);
    } else {
      int_sub(vx, vx, ux, uvlen);
      gfp_sub(x2, x2, x1, c, len);
    }
    if ((ux[uvlen-1] == 0) && (vx[uvlen-1] == 0)) uvlen--;
  }
  
  if (int_is1(ux, len)) int_copy(r, x1, len);
  return MSPECC_NO_ERROR;
}
