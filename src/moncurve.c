///////////////////////////////////////////////////////////////////////////////
// moncurve.c: Point arithmetic and scalar mul. on Montgomery curves.        //
// This file is part of SECC430, a Scalable ECC implementation for MSP430.   //
// Version 0.7.0 (2023-06-24), see <http://www.cryptolux.org/> for updates.  //
// License: GPLv3 (see LICENSE file), other licenses available upon request. //
// ------------------------------------------------------------------------- //
// This program is free software: you can redistribute it and/or modify it   //
// under the terms of the GNU General Public License as published by the     //
// Free Software Foundation, either version 3 of the License, or (at your    //
// option) any later version. This program is distributed in the hope that   //
// it will be useful, but WITHOUT ANY WARRANTY; without even the implied     //
// warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the  //
// GNU General Public License for more details. You should have received a   //
// copy of the GNU General Public License along with this program. If not,   //
// see <http://www.gnu.org/licenses/>.                                       //
///////////////////////////////////////////////////////////////////////////////


#include <stdio.h>
#include <stdlib.h>
#include "intarith.h"
#include "gfparith.h"
#include "moncurve.h"
#include "tedcurve.h"
#include "ecdparam.h"


#ifdef MSPECC_USE_VLA  // use variable-length arrays for field elements
#define _len len       // requires "Allow VLA" in C/C++ Compiler Options
#else  // use maximum-length arrays
#define _len (MSPECC_MAX_LEN/WSIZE)
#endif


#if (WSIZE == 16)
#define GET_BIT(k, i) ((k[((i)>>4)] >> ((i) & 0x0F)) & 1)
#else  // 32 bits
#define GET_BIT(k, i) ((k[((i)>>5)] >> ((i) & 0x1F)) & 1)
#endif


#if (WSIZE == 16)
static const Word SECC_INV_MASK[16] = {                           \
  0x5F58, 0xE072, 0x28DB, 0x1703, 0xBC96, 0x22E6, 0x97C4, 0xA158, \
  0x646A, 0xCED0, 0x2D36, 0xE628, 0x9A79, 0x4908, 0x4D46, 0x76F9 };
#if ((MSPECC_MAX_LEN/WSIZE) > 16)
#error "Multiplicative mask for secure inversion must be extended!"
#endif // MSPECC_MAX_LEN/WSIZE) > 16)
#else  // 32 bits
static const Word SECC_INV_MASK[8] = {            \
  0xE0725F58, 0x170328DB, 0x22E6BC96, 0xA15897C4, \
  0xCED0646A, 0xE6282D36, 0x49089A79, 0x76F94D46  };
#if ((MSPECC_MAX_LEN/WSIZE) > 8)
#error "Multiplicative mask for secure inversion must be extended!"
#endif // MSPECC_MAX_LEN/WSIZE) > 8)
#endif


/*****************************************************************************/
/* Copy the coordinates of a projective point.                               */
/*****************************************************************************/

void mon_copy(PROPOINT *r, const PROPOINT *p, int len)
{
  int_copy(r->x, p->x, len);
  if ((r->y != NULL) && (p->y != NULL)) int_copy(r->y, p->y, len);
  if ((r->z != NULL) && (p->z != NULL)) int_copy(r->z, p->z, len);
  if ((r->z != NULL) && (p->z == NULL)) int_set(r->z, 1, len);
}


/*****************************************************************************/
/* Differential point addition P = P + Q on a Montgomery curve. The points P */
/* and Q are expected to be given in standard projective coordinates of the  */
/* form (X,Z), i.e. the Y-coordinate is not needed. The parameter 'xd' is    */
/* the x-coordinate of the difference D = P - Q, which is normally the same  */
/* as the x-coordinate of the "original" base point when performing a scalar */
/* multiplication according to the so-called Montgomery ladder.              */
/*****************************************************************************/

void mon_add(PROPOINT *p, const PROPOINT *q, const Word *xd,
             const ECDPARAM *m)
{
  int len = m->len; Word c = m->c;
  Word *xp = p->x, *zp = p->z, *xq = q->x, *zq = q->z;
  Word *t1 = p->y, *t2 = p->slack, *prod = &(p->slack[len]);
  (void) prod;  // to silence a warning
  
  gfp_add(t1, xp, zp, c, len);          // t1 := xp+zp;
  gfp_sub(t2, xp, zp, c, len);          // t2 := xp-zp;
  gfp_sub(xp, xq, zq, c, len);          // xr := xq-zq;
  gfp_mul(zp, t1, xp, c, len);          // zr := t1*xr;
  gfp_add(t1, xq, zq, c, len);          // t1 := xq+zq;
  gfp_mul(xp, t1, t2, c, len);          // xr := t1*t2;
  gfp_sub(t1, xp, zp, c, len);          // t1 := xr-zr;
  gfp_add(t2, xp, zp, c, len);          // t2 := xr+zr;
  gfp_sqr(xp, t2, c, len);              // xr := t2*t2;
  gfp_sqr(t2, t1, c, len);              // t2 := t1*t1;
  gfp_mul(zp, xd, t2, c, len);          // zr := xd*t2;
}


/*****************************************************************************/
/* Point doubling P = 2*P on a Montgomery curve. The point P is expected to  */
/* be given in standard projective coordinates of the form (X,Z), i.e. the   */
/* Y-coordinate is not needed.                                               */
/*****************************************************************************/

void mon_double(PROPOINT *p, const ECDPARAM *m)
{
  int len = m->len; Word c = m->c;
  Word *xp = p->x, *zp = p->z;
  Word *t1 = p->y, *t2 = p->slack, *prod = &(p->slack[len]);
  (void) prod;  // to silence a warning
  
  gfp_add(t1, xp, zp, c, len);          // t1 := xp+zp;
  gfp_sqr(t2, t1, c, len);              // t2 := t1*t1;
  gfp_sub(t1, xp, zp, c, len);          // t1 := xp-zp;
  gfp_sqr(zp, t1, c, len);              // zr := t1*t1;
  gfp_mul(xp, t2, zp, c, len);          // xr := t2*zr;
  gfp_sub(t1, t2, zp, c, len);          // t1 := t2-zr;
  gfp_mul32(t2, t1, m->a24, c, len);    // t2 := t1*a24;
  gfp_add(t2, t2, zp, c, len);          // t2 := t2+zr;
  gfp_mul(zp, t1, t2, c, len);          // zr := t1*t2;
}


/*****************************************************************************/
/* Check whether a point P on a Montgomery curve has low order by computing  */
/* the (projective) product R = 8*P and checking whether its Z-coordinate is */
/* 0. Only the (affine) x-coordinate of the base point P is needed. This     */
/* check can help to prevent certain kinds of side-channel attack, including */
/* the one described in the paper "To Infinity and Beyond: Combined Attack   */
/* on ECC Using Points of Low Order" (Proceedings of CHES 2011).             */
/*****************************************************************************/

int mon_check_order(PROPOINT *r, const Word *xp, const ECDPARAM *m)
{
  int len = m->len;
  Word c = m->c;
  
  // initialize point R with x-coordinate of P
  int_copy(r->x, xp, len); int_set(r->z, 1, len);
  
  // co-factor multiplication (three doublings)
  mon_double(r, m);
  mon_double(r, m);
  mon_double(r, m);
  
  // check whether Z-coordinate equals 0 or prime p (we do not need to check
  // whether Z is 2*p since the result of a field-arithmetic operation is
  // always smaller than 2*p)
  if (int_is0(r->z, len) || gfp_isp(r->z, c, len))
    return MSPECC_ERR_INVALID_POINT;
  
  return MSPECC_NO_ERROR;
}


/*****************************************************************************/
/* Montgomery ladder for scalar multiplication R = k*P on a Montgomery       */
/* curve. This implementation is optimized for secret scalars as specified   */
/* in the paper "Curve25519: New Diffie-Hellman Speed Records" (Proceedings  */
/* of PKC 2005), which means the three least-significant bits are always 0   */
/* and, thus, three doublings are sufficient to process these bits. Only the */
/* (affine) x-coordinate of the base point P is needed. The result R is      */
/* given in projective coordinates of the form (X,Z), i.e. the Y-coordinate  */
/* is not computed.                                                          */
/*****************************************************************************/

void mon_mul_ladder(PROPOINT *r, const Word *k, const Word *xp,
                    const ECDPARAM *m)
{
  int ki, len = m->len, i = WSIZE*len - 1;
  Word tmp[3*_len]; // temporary space for three gfp elements
  PROPOINT q = { tmp, &tmp[len], &tmp[2*len], NULL, r->slack };
  PROPOINT *t[2] = { r, &q };
  
  // find position of first non-zero bit in 'k'. This does not introduce a
  // vulnerability to side-channel attacks since, normally, the position of the
  // leading "1" in a secret scalar is fixed. For example, a secret scalar in
  // Curve25519 consists of 32 bytes, whereby the highest bit of the most
  // significant byte is always "0" and the second-highest bit is always "1".
  while ((GET_BIT(k, i) == 0) && (i > 0)) i--;
  
  // initialize T[0] with (X,Z) = (xp,1)  
  int_copy(t[0]->x, xp, len);
  int_set(t[0]->z, 1, len);
  // initialize T[1] with 2*T[0]
  mon_copy(t[1], t[0], len);
  mon_double(t[1], m);
  
  // simple left-to-right Montgomery ladder
  for(i = i - 1; i >= 0; i--) {
    ki = GET_BIT(k, i);
    mon_add(t[1-ki], t[ki], xp, m);
    mon_double(t[ki], m);
  }
  
  // At the end of the ladder, the point T[0] contains the X and Z coordinate
  // of R = k*P, while T[1] contains the X and Z coordinate of T[0]+P, i.e. we
  // have T[1] = k*P + P = (k+1)*P. These four coordinates, along with the x
  // and y-coordinate of the base point P, allow a recovery of the Y-coordinate
  // of R. Hence, we store the X and Z-coordinate of T[1] in the Word arrays
  // addressed by the pointers 'y' and 'slack' of the PROPOINT structure for R.
  int_copy(r->y, q.x, len);
  int_copy(r->slack, q.z, len);
}


/*****************************************************************************/
/* Montgomery ladder for scalar multiplication R = k*P on a Montgomery       */
/* curve. Only the (affine) x-coordinate of the base point P is needed. The  */
/* result R is given in projective coordinates of the form (X,Z), i.e. the   */
/* Y-coordinate is not computed. However, besides the X and Z-coordinate of  */
/* R = k*P, also the X and Z-coordinate of the point Q = R + P = (k+1)*P is  */
/* computed by the Montgomery latter. These two coordinates are stored in    */
/* the Word arrays referenced by the two elements 'y' and 'slack' of the     */
/* PROPOINT structure representing the result R. The X and Z-coordinate of   */
/* Q, along with the X and Z-coordinate of R and the affine x and            */
/* y-coordinate of the base point P, can be used to recover the Y coordinate */
/* of R.                                                                     */
/*****************************************************************************/

void mon_mul_ladder_consttime(PROPOINT *r, const Word *k, const Word *xp,
                              const ECDPARAM *m)
{
  int i, ki, len = m->len;
  Word tmp[3*_len]; // temporary space for three gfp elements
  PROPOINT q = { tmp, &tmp[len], &tmp[2*len], NULL, r->slack };
  PROPOINT *t[2] = { r, &q };
  
  // A naive implementation of the Montgomery ladder first scans k to find the
  // leading "1" and then initializes T[0] with P and T[1] with 2*P. However,
  // this approach leaks the bitlength of k through the execution time unless
  // the leading "1" is always at a fixed position. In order to achieve a
  // constant execution time independent of k, T[0] has to be initialized with
  // the point at infinity (by setting (X,Z) to (1,0)) and T[1] with P. In this
  // case, leading "0" bits in k do not change T[0] and T[1], which means the
  // overall number of loop iterations is always constant, namely 16*'len'.
  
  // initialize T[0] with (X,Z) = (1,0)  
  int_set(t[0]->x, 1, len);
  int_set(t[0]->z, 0, len);
  // initialize T[1] with (X,Z) = (xp,1)  
  int_copy(t[1]->x, xp, len);
  int_set(t[1]->z, 1, len);
  
  // simple left-to-right Montgomery ladder
  for(i = WSIZE*len - 1; i >= 0; i--) {
    ki = GET_BIT(k, i);
    mon_add(t[1-ki], t[ki], xp, m);
    mon_double(t[ki], m);
  }
  
  // At the end of the ladder, the point T[0] contains the X and Z coordinate
  // of R = k*P, while T[1] contains the X and Z coordinate of T[0]+P, i.e. we
  // have T[1] = k*P + P = (k+1)*P. These four coordinates, along with the x
  // and y-coordinate of the base point P, allow a recovery of the Y-coordinate
  // of R. Hence, we store the X and Z-coordinate of T[1] in the Word arrays
  // addressed by the pointers 'y' and 'slack' of the PROPOINT structure for R.
  int_copy(r->y, q.x, len);
  int_copy(r->slack, q.z, len);
}


/*****************************************************************************/
/* Conversion of a point P given in standard projective coordinates of the   */
/* form (X,Z) into a point R in standard affine coordinates, whereby only    */
/* x-coordinate is computed. This conversion requires the inversion of Z,    */
/* which can leak information about the secret scalar when implemented in a  */
/* straightforward way using e.g. the extended Euclidean algorithm; see the  */
/* paper "Projective Coordinates Leak" (Proceedings of EUROCRYPT 2004). To   */
/* prevent this leakage, a simple masking technique is applied, which means  */
/* that, instead of inverting Z "directly", it is first multiplied by a      */
/* field element 'u' that is unknown to the attacker, and then the product   */
/* Z*u is inverted. Finally, the inverse (Z*u)^-1 is multiplied by 'u' to    */
/* obtain Z^-1. This approach for secure inversion is described in Section 6 */
/* of the paper "SPA Vulnerabilities of the Binary Extended Euclidean        */
/* Algorithm" (Journal of Cryptographic Engineering, 2016).                  */
/*****************************************************************************/

int mon_proj_affine(PROPOINT *r, const PROPOINT *p, const ECDPARAM *m)
{
  int err, len = m->len; Word c = m->c;
  Word *xp = p->x, *yp = p->y, *zp = p->z;
  Word *xr = r->x, *yr = r->y, *zr = r->z;
  Word *t1 = r->slack, *prod = &(r->slack[len]);
  (void) prod;  // to silence a warning
  
  // "masked" inversion of Z to thwart timing attacks
  gfp_mul(t1, zp, SECC_INV_MASK, c, len);
  err = gfp_inv(t1, t1, c, len);
  if (err != MSPECC_NO_ERROR) return err;
  gfp_mul(zr, t1, SECC_INV_MASK, c, len);
  
  // get least non-negative residue of x = X*(1/Z)
  gfp_mul(t1, xp, zr, c, len);
  gfp_lnr(xr, t1, c, len);
  
  if ((yp != NULL) && (yr != NULL)) {
    // get least non-negative residue of x = X*(1/Z)
    gfp_mul(t1, yp, zr, c, len);
    gfp_lnr(yr, t1, c, len);
  }
  
  // overwrite Z with 1
  int_set(zr, 1, len);
  
  return MSPECC_NO_ERROR;
}


/*****************************************************************************/
/* Recover the Y-coordinate of a point Q obtained by the Montgomery ladder.  */
/* Besides the X and Z-coordinate of Q = k*P, also the X and Z-coordinate of */
/* the point Q + P = (k+1)*P as well as the affine x and y-coordinate of the */
/* base point P are needed. The X and Z-coordinate of Q + P are expected to  */
/* be contained in Word arrays referenced by the elements 'y' and 't' of the */
/* PROPOINT structure representing the point Q. This implementation of the   */
/* Y-coordinate recovery is optimized for Montgomery curves with B = 1.      */
/*****************************************************************************/

void mon_recover_y(PROPOINT *r, const PROPOINT *q, const PROPOINT *p,
                   const ECDPARAM *m)
{
  int len = m->len; Word c = m->c;
  Word tmp[3*_len]; // temporary space for three gfp elements
  Word *t1 = tmp, *t2 = &tmp[len], *t3 = &tmp[2*len];
  Word *x1 = q->x, *z1 = q->z, *x2 = q->y, *z2 = q->slack;
  Word *xr = r->x, *yr = r->y, *zr = r->z, *xp = p->x, *yp = p->y;
  Word *prod = &(r->slack[len]);
  (void) prod;  // to silence a warning
  
  gfp_mul(t1, xp, x1, c, len);          // t1 := xp*x1;
  gfp_sub(t1, t1, z1, c, len);          // t1 := t1-z1;
  gfp_mul(t2, z1, xp, c, len);          // t2 := z1*xp;
  gfp_sub(t2, x1, t2, c, len);          // t2 := x1-t2;
  gfp_mul(t3, z2, t1, c, len);          // t3 := z2*t1;
  gfp_mul(t1, x2, t2, c, len);          // t1 := x2*t2;
  gfp_add(t2, t3, t1, c, len);          // t2 := t3+t1;
  gfp_sub(t3, t3, t1, c, len);          // t3 := t3-t1;
  gfp_mul(t1, x2, yp, c, len);          // t1 := x2*yp;
  gfp_mul(yr, t2, t3, c, len);          // yr := t2*t3;
  gfp_add(t3, z2, z2, c, len);          // t3 := 4*b;
  gfp_add(t2, t3, t3, c, len);          // t2 := t3*z2;
  gfp_mul(t3, t2, t1, c, len);          // t3 := t2*t1;  
  gfp_mul(t2, t3, z1, c, len);          // t2 := t3*z1;
  gfp_mul(zr, t2, z1, c, len);          // zr := t2*z1;
  gfp_mul(xr, t2, x1, c, len);          // xr := t2*x1;
}


/*****************************************************************************/
/* Variable-base scalar multiplication R = k*P on a Montgomery curve,        */
/* including some tests to ensure the validity of inputs and outputs. Only   */
/* the (affine) x-coordinate of the base point P is required. The result R   */
/* given in affine coordinates and also consists of only the x-coordinate,   */
/* i.e. the y-coordinate is not computed.                                    */
/*****************************************************************************/

int mon_mul_varbase(Word *r, const Word *k, const Word *xp,
                    const ECDPARAM *m)
{
  int err, len = m->len;
  Word tmp[6*_len]; // temporary space for six gfp elements
  PROPOINT q = { tmp, &tmp[len], &tmp[2*len], NULL, &tmp[3*len] };  
  
  // set r to 0 when k is 0 (should normally never happen)
  if (int_is0(k, len)) { int_set(r, 0, len); return MSPECC_ERR_INVALID_SCALAR; }
  
  // check the order of P to prevent the attack described in "To Infinity and
  // Beyond: Combined Attack on ECC Using Points of Low Order" (CHES 2011)
  // err = mon_check_order(&q, xp, m);
  // if (err != MSPECC_NO_ERROR) { int_set(r, 0, len); return err; }
  
  // perform the scalar multiplication (Montgomery ladder)
  mon_mul_ladder(&q, k, xp, m);
  
  // convert result from projective to affine coordinates
  err = mon_proj_affine(&q, &q, m);
  if (err != MSPECC_NO_ERROR) { int_set(r, 0, len); return err; }
  
  // assign x-coordinate of affine point to output
  int_copy(r, q.x, len);
  
  return MSPECC_NO_ERROR;
}


int mon_mul_fixbase(Word *r, const Word *k, const ECDPARAM *m)
{
  int err, len = m->len, c = m->c;
  Word tmp[8*_len];
  PROPOINT q = { tmp, &tmp[len], &tmp[2*len], &tmp[3*len], &tmp[5*len] };
  Word *prod = &(q.slack[len]);
  (void) prod;  // to silence a warning
  
  // set r to 0 when k is 0 (should normally never happen)
  if (int_is0(k, len)) { int_set(r, 0, len); return MSPECC_ERR_INVALID_SCALAR; }
  
  // perform scalar multiplication via fixed-base comb method
  ted_mul_comb4b(&q, k, m);
  
  // from twisted Edwards curve to Montgomery curve u = (Z+Y)/(Z-Y)
  gfp_sub(q.extra, q.z, q.y, c, len);
  gfp_add(q.slack, q.z, q.y, c, len);
  
  // "masked" inversion of Z-Y to thwart timing attacks
  gfp_mul(q.x, q.extra, INV_MASK, c, len);
  err = gfp_inv(q.x, q.x, c, len);
  if (err != MSPECC_NO_ERROR) return err;
  gfp_mul(q.extra, q.x, INV_MASK, c, len);
   
  // get least non-negative residue of u = (Z+Y)/(Z-Y)
  gfp_mul(q.x, q.slack, q.extra, c, len);
  gfp_lnr(r, q.x, c, len);
  
  return MSPECC_NO_ERROR;
}



/*****************************************************************************/
/* Convert a point P on a Montgomery curve to the corresponding point R on   */
/* the birationally equivalent twisted Edwards curve. The point P is         */
/* expected to to be given in standard projective coordinates and the result */
/* R is also given in standard projective coordinates. This conversion       */
/* requires a pre-computed constant c = SquareRoot(-(A+2)/B) where A, B are  */
/* the curve parameters of the Montgomery curve. Note that (A+2)/B equals    */
/* the parameter a of the birationally-equivalent twisted Edwards curve.     */
/*****************************************************************************/

void mon_to_ted(PROPOINT *r, const PROPOINT *p, const ECDPARAM *m)
{
  int len = m->len; Word c = m->c;
  Word tmp[2*_len]; // temporary space for two gfp elements
  Word *t1 = tmp, *t2 = &tmp[len], *t3 = r->slack, *prod = &(r->slack[len]);
  Word *xm = p->x, *ym = p->y, *zm = p->z;
  Word *xt = r->x, *yt = r->y, *zt = r->z; 
  (void) prod;  // to silence a warning
  
  gfp_add(t1, xm, zm, c, len);
  gfp_sub(t2, xm, zm, c, len);
  gfp_mul(t3, xm, m->rma, c, len);
  gfp_mul(xt, t3, t1, c, len);          // xt := c*xm*(xm + zm);
  gfp_mul(zt, ym, t1, c, len);          // zt := ym*(xm + zm);
  gfp_mul(yt, ym, t2, c, len);          // yt := ym*(xm - zm);
}


void mon_test25519(void)
{
  #if (WSIZE == 16)
  // Word x[256/WSIZE] = {                                             \
  //   0x0009, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, \
  //   0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000 };
  // Word k[256/WSIZE] = {                                             \
  //   0xAB58, 0x7E08, 0x4A62, 0x4B8A, 0xE179, 0x8B7F, 0x8083, 0xE60E, \
  //   0x3B6F, 0x29B1, 0x1826, 0xFDB6, 0x2F1C, 0x278B, 0x88FF, 0x6BE0 };
  // Test vector from https://tools.ietf.org/html/draft-irtf-cfrg-curves-10
  // r = 0x4F2B886F147EFCAD4D67785BC843833F3735E4ECC2615BD3B4C17D7B7DDB9EDE
  Word x[256/WSIZE] = {                                             \
    0xDBE6, 0x6768, 0x3058, 0xDB30, 0x9435, 0xA4C1, 0xB124, 0x7C5F, \
    0x6672, 0xEC24, 0xB326, 0x3B35, 0xA910, 0xA603, 0xABD0, 0x4C1C };
  Word k[256/WSIZE] = {                                             \
    0x46A0, 0x6BE3, 0x52F0, 0x9D7C, 0x163B, 0x4B15, 0x4682, 0xDD5E, \
    0x1462, 0x0A4C, 0xFCC1, 0x185A, 0x6A50, 0x4422, 0x44BA, 0x449A };
  #else  // WSIZE == 32
  // Word x[256/WSIZE] = {                             \
  //   0x00000009, 0x00000000, 0x00000000, 0x00000000, \
  //   0x00000000, 0x00000000, 0x00000000, 0x00000000 };
  // Word k[256/WSIZE]  = {                            \
  //   0x7E08AB58, 0x4B8A4A62, 0x8B7FE179, 0xE60E8083, \
  //   0x29B13B6F, 0xFDB61826, 0x278B2F1C, 0x6BE088FF };
  // r = 0x4F2B886F147EFCAD4D67785BC843833F3735E4ECC2615BD3B4C17D7B7DDB9EDE
  // Test vector from https://tools.ietf.org/html/draft-irtf-cfrg-curves-10
  Word x[256/WSIZE] = {                             \
    0x6768DBE6, 0xDB303058, 0xA4C19435, 0x7C5FB124, \
    0xEC246672, 0x3B35B326, 0xA603A910, 0x4C1CABD0 };
  Word k[256/WSIZE] = {                             \
    0x6BE346A0, 0x9D7C52F0, 0x4B15163B, 0xDD5E4682, \
    0x0A4C1462, 0x185AFCC1, 0x44226A50, 0x449A44BA };
  #endif
  
  Word r[256/WSIZE];
  int len = 256/WSIZE;
  
  // pruning: make sure that scalar k is a valid scalar
  k[len-1] &= (((Word) -1L) >> 1);            // 0x7F..FF 
  k[len-1] |= ((Word) (1UL << (WSIZE - 2)));  // 0x40..00 
  k[0]     &= ((Word) -8L);                   // 0xFF..F8
  
  #ifdef MSPECC_DEBUG_PRINT
  int_print("x = ", x, len);
  int_print("k = ", k, len);
  #endif
  
  mon_mul_varbase(r, k, x, &CURVE25519);
  // F1611: 12,191,721 cycles, FR5969: 10,878,075
  
  #ifdef MSPECC_DEBUG_PRINT
  int_print("r = ", r, len);
  #endif
  // correct result:
  // r = 0x5285A2775507B454F7711C4903CFEC324F088DF24DEA948E90C6E99D3755DAC3
}
