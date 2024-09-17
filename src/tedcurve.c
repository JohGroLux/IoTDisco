///////////////////////////////////////////////////////////////////////////////
// tedcurve.c: Point arithmetic and scalar mul. on twisted Edwards curves.   //
// This file is part of SECC430, a Scalable ECC implementation for MSP430.   //
// Version 1.0.1 (2023-06-24), see <http://www.cryptolux.org/> for updates.  //
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
#include "config.h"
#include "intarith.h"
#include "gfparith.h"
#include "ecdparam.h"
#include "moncurve.h"
#include "tedcurve.h"


#ifdef MSPECC_USE_VLA  // use variable-length arrays for field elements
#define _len len       // requires "Allow VLA" in C/C++ Compiler Options
#else
#define _len MSPECC_MAX_LEN
#endif


#define get_bit(k, i) ((k[((i)>>4)] >> ((i)&0xF)) & 1)


#if (MSPECC_MAX_LEN <= 16)
static const UINT16 INV_MASK[16] = { 0x5F58, 0xE072, 0x28DB, 0x1703, 0xBC96, \
  0x22E6, 0x97C4, 0xA158, 0x646A, 0xCED0, 0x2D36, 0xE628, 0x9A79, 0x4908,    \
  0x4D46, 0x76F9 };
#endif


/*****************************************************************************/
/* Print a point P given in projective coordinates. The parameter 'num' can  */
/* be either 2 or 3 and specifies the number of coordinates to be printed.   */
/* However, the x and y-coordinate are printed in any case, regardless of    */
/* the actual value of 'num'. If 'num' is 3 then the z-coordinate is also    */
/* printed.                                                                  */
/*****************************************************************************/

void ted_print(const PROPOINT *p, int len, int num)
{
  int_print("x = ", p->x, len);
  int_print("y = ", p->y, len);
  if ((num == 3) && (p->z != NULL)) int_print("z = ", p->z, len);
}


/*****************************************************************************/
/* Set an affine point P to the neutral element (0,1).                       */
/*****************************************************************************/

void ted_set0_aff(AFFPOINT *p, int len)
{
  int_set(p->x, 0, len);
  int_set(p->y, 1, len);
}


/*****************************************************************************/
/* Set a projective point P to the neutral element (0,1,1).                  */
/*****************************************************************************/

void ted_set0_pro(PROPOINT *p, int len)
{
  int_set(p->x, 0, len);
  int_set(p->y, 1, len);
  int_set(p->z, 1, len);
}


/*****************************************************************************/
/* Copy certain coordinates of a source point P to a destination point R.    */
/* The parameter 'num' can be either 2, 3, or 5 and specifies the number of  */
/* coordinates to be copied. However, the x and y-coordinate are copied in   */
/* any case, regardless of the actual value of 'num'. If 'num' is 2 then the */
/* z-coordinate of R is set to 1. On the other hand, if 'num' is 3 or 5, the */
/* z-coordinate is also copied. Furthermore, if 'num' is 5, it is assumed    */
/* that P and R are points in extended projective coordinates and, thus, the */
/* coordinates referenced by the element 'extra' of the PROPOINT structure   */
/* are copied as well.                                                       */
/*****************************************************************************/

void ted_copy(PROPOINT *r, const PROPOINT *p, int len, int num)
{
  int_copy(r->x, p->x, len);
  int_copy(r->y, p->y, len);
  if (num == 2) int_set(r->z, 1, len);
  if ((num == 3) || (num == 5)) int_copy(r->z, p->z, len);
  if (num == 5) int_copy(r->extra, p->extra, 2*len);
}


/*****************************************************************************/
/* Mixed point addition P = P + Q on a twisted Edwards curve. The point P is */
/* expected to be given in extended projective coordinates of the form       */
/* (X,Y,Z,E,H) with E*H = T = X*Y/Z. Note that the coordinates E and H can   */
/* only be accessed through the element 'extra' of the PROPOINT structure.   */
/* The point Q is expected to be given in extended affine coordinates of the */
/* form (u,v,w) where u = (x+y)/2, v = (y-x)/2, and w = d*x*y.               */
/*****************************************************************************/

void ted_add(PROPOINT *p, const PROPOINT *q, const ECDPARAM *m)
{
  int len = m->len; UINT16 c = m->c;
  UINT16 *x1 = p->x, *y1 = p->y, *z1 = p->z;
  UINT16 *e1 = p->extra, *h1 = &(p->extra[len]);
  UINT16 *t1 = p->slack, *prod = &(p->slack[len]);
  const UINT16 *u2 = q->x, *v2 = q->y, *w2 = q->z;
  
  gfp_mul(t1, e1, h1, c, len);          // t1 := e1*h1;
  gfp_sub(e1, y1, x1, c, len);          // e3 := y1-x1;
  gfp_add(h1, y1, x1, c, len);          // h3 := y1+x1;
  gfp_mul(x1, e1, v2, c, len);          // x3 := e3*v2;
  gfp_mul(y1, h1, u2, c, len);          // y3 := h3*u2;
  gfp_sub(e1, y1, x1, c, len);          // e3 := y3-x3;
  gfp_add(h1, y1, x1, c, len);          // h3 := y3+x3;
  gfp_mul(x1, t1, w2, c, len);          // x3 := t1*w2;
  gfp_sub(t1, z1, x1, c, len);          // t1 := z1-x3;
  gfp_add(x1, z1, x1, c, len);          // x3 := z1+x3;
  gfp_mul(z1, t1, x1, c, len);          // z3 := t1*x3;
  gfp_mul(y1, x1, h1, c, len);          // y3 := x3*h3;
  gfp_mul(x1, e1, t1, c, len);          // x3 := e3*t1;
}


/*****************************************************************************/
/* Extended point doubling P = 2*P on a twisted Edwards curve. The point P   */
/* is expected to be given in extended projective coordinates of the form    */
/* (X,Y,Z,E,H) with E*H = T = X*Y/Z. Note that the coordinates E and H can   */
/* only be accessed through the element 'extra' of the PROPOINT structure.   */
/*****************************************************************************/

void ted_double(PROPOINT *p, const ECDPARAM *m)
{
  int len = m->len; UINT16 c = m->c;
  UINT16 *x1 = p->x, *y1 = p->y, *z1 = p->z;
  UINT16 *e1 = p->extra, *h1 = &(p->extra[len]);
  UINT16 *t1 = p->slack, *prod = &(p->slack[len]);
  
  gfp_sqr(e1, x1, c, len);              // e3 := x1*x1;
  gfp_sqr(h1, y1, c, len);              // h3 := y1*y1;
  gfp_sub(t1, e1, h1, c, len);          // t1 := e3-h3;
  gfp_add(h1, e1, h1, c, len);          // h3 := e3+h3;
  gfp_add(x1, x1, y1, c, len);          // x3 := x1+y1;
  gfp_sqr(e1, x1, c, len);              // e3 := x3*x3;
  gfp_sub(e1, h1, e1, c, len);          // e3 := h3-e3;
  gfp_sqr(y1, z1, c, len);              // y3 := z1*z1;
  gfp_add(y1, y1, y1, c, len);          // y3 := 2*y3;
  gfp_add(y1, t1, y1, c, len);          // y3 := t1+y3;
  gfp_mul(x1, e1, y1, c, len);          // x3 := e3*y3;
  gfp_mul(z1, y1, t1, c, len);          // z3 := y3*t1;
  gfp_mul(y1, t1, h1, c, len);          // y3 := t1*h3;
}


/*****************************************************************************/
/* Projective point addition R = R + P on a twisted Edwards curve. Both R    */
/* and P are expected to be given in standard projective coordinates of the  */
/* form (X,Y,Z).                                                             */
/*****************************************************************************/

void ted_add_pro(PROPOINT *r, const PROPOINT *p, const ECDPARAM *m)
{
  int len = m->len; UINT16 c = m->c;
  UINT16 *x1 = r->x, *y1 = r->y, *z1 = r->z;  
  UINT16 *t1 = r->extra, *t2 = &(r->extra[len]);
  UINT16 *t3 = r->slack, *prod = &(r->slack[len]);
  UINT16 *x2 = p->x, *y2 = p->y, *z2 = p->z;
  
  gfp_add(t1, x1, y1, c, len);          // t1 := x1+y1;
  gfp_add(t2, x2, y2, c, len);          // t2 := x2+y2;
  gfp_mul(t3, t1, t2, c, len);          // t3 := t1*t2;
  gfp_mul(t1, z1, z2, c, len);          // t1 := z1*z2;
  gfp_mul(z1, x1, x2, c, len);          // z3 := x1*x2;
  gfp_mul(x1, y1, y2, c, len);          // x3 := y1*y2;
  gfp_add(y1, z1, x1, c, len);          // y3 := z3+x3;
  gfp_mul(t2, z1, x1, c, len);          // t2 := z3*x3;
  gfp_mul(x1, t2, m->d, c, len);        // x3 := d*t2;
  gfp_sqr(t2, t1, c, len);              // t2 := t1^2;
  gfp_sub(z1, t3, y1, c, len);          // z3 := t3-y3;
  gfp_sub(t3, t2, x1, c, len);          // t3 := t2-x3;
  gfp_add(x1, t2, x1, c, len);          // x3 := t2+x3;
  gfp_mul(t2, x1, y1, c, len);          // t2 := x3*y3;
  gfp_mul(y1, t1, t2, c, len);          // y3 := t1*t2;
  gfp_mul(t2, t3, z1, c, len);          // t2 := t3*z3;
  gfp_mul(z1, x1, t3, c, len);          // z3 := x3*t3;
  gfp_mul(x1, t1, t2, c, len);          // x3 := t1*t2;
}


/*****************************************************************************/
/* Conversion of a point P given in standard affine coordinates (x,y) into a */
/* point R in extended affine coordinates of the form (u,v,w) where u =      */
/* (x+y)/2, v = (y-x)/2, and w = d*x*y.                                      */
/*****************************************************************************/

void ted_affine_extaff(PROPOINT *r, const AFFPOINT *p, const ECDPARAM *m)
{
  int len = m->len; UINT16 c = m->c;
  UINT16 *x = p->x, *y = p->y;
  UINT16 *u = r->x, *v = r->y, *w = r->z;
  UINT16 *t1 = r->slack, *prod = &(r->slack[len]);
  
  gfp_add(t1, x, y, c, len);
  gfp_hlv(u, t1, c, len);
  gfp_sub(t1, y, x, c, len);
  gfp_hlv(v, t1, c, len);
  gfp_mul(t1, x, y, c, len);
  gfp_mul(w, t1, m->d, c, len);
}


/*****************************************************************************/
/* Conversion of a point P given in extended affine coordinates of the form  */
/* (u,v,w) where u = (x+y)/2, v = (y-x)/2, and w = d*x*y into a point R in   */
/* extended projective coordinates of the form (X,Y,Z,E,H) with E*H = T =    */
/* X*Y/Z.                                                                    */
/*****************************************************************************/

void ted_extaff_extpro(PROPOINT *r, const PROPOINT *p, const ECDPARAM *m)
{
  int len = m->len; UINT16 c = m->c;
  UINT16 *u = p->x, *v = p->y;
  UINT16 *x = r->x, *y = r->y, *z = r->z;
  
  gfp_add(x, u, u, c, len);             // x := 2*u;
  gfp_add(y, v, v, c, len);             // y := 2*v;
  gfp_add(y, y, x, c, len);             // y := y+x;
  gfp_hlv(y, y, c, len);                // y := y/2;
  gfp_sub(x, x, y, c, len);             // x := x-y;
  int_set(z, 1, len);
  if (r->extra != NULL)
  {
    int_copy(r->extra, x, len);
    int_copy(&(r->extra[len]), y, len);
  }
}


/*****************************************************************************/
/* Conversion of a point P given in standard affine coordinates (x,y) into a */
/* point R in (extended) projective coordinates of the form (X,Y,Z).         */
/*****************************************************************************/

void ted_aff_to_pro(PROPOINT *r, const AFFPOINT *p, const ECDPARAM *m)
{
  int len = m->len;
  UINT16 *xp = p->x, *yp = p->y;
  UINT16 *xr = r->x, *yr = r->y, *zr = r->z;
  
  int_copy(xr, xp, len);
  int_copy(yr, yp, len);
  int_set(zr, 1, len);
  if (r->extra != NULL)
  {
    int_copy(r->extra, xr, len);
    int_copy(&(r->extra[len]), yr, len);
  }
}


/*****************************************************************************/
/* Validate a point P given in projective coordinates (i.e. check whether    */
/* X, Y, and Z coordinate of P satisfy the projective curve equation of the  */
/* form (Y^2 - X^2)*Z^2 = Z^4 + d*X^2*Y^2.                                   */
/*****************************************************************************/

int ted_validate(const PROPOINT *p, const ECDPARAM *m)
{
  int r, len = m->len; UINT16 c = m->c;
  UINT16 tmp[3*_len]; // temporary space for three gfp elements
  UINT16 *t1 = tmp, *t2 = &tmp[len], *t3 = &tmp[2*len];
  UINT16 *t4 = (UINT16 *) p->slack, *prod = (UINT16 *) &(p->slack[len]);
  UINT16 *x = p->x, *y = p->y, *z = p->z;
  
  // compute t1 = (Y^2 - X^2)*Z^2 and t2 = Z^4 + d*X^2*Y^2
  gfp_sqr(t1, x, c, len);               // t1 := X^2;
  gfp_sqr(t2, y, c, len);               // t2 := Y^2;
  gfp_mul(t3, t1, t2, c, len);          // t3 := t1*t2;
  gfp_sub(t2, t2, t1, c, len);          // t2 := t2-t1;
  gfp_mul(t4, t3, m->d, c, len);        // t4 := t3*d;
  gfp_sqr(t3, z, c, len);               // t3 := Z^2;
  gfp_mul(t1, t3, t2, c, len);          // t1 := t3*t2
  gfp_sqr(t2, t3, c, len);              // t2 := t3^2
  gfp_add(t2, t2, t4, c, len);          // t2 := t2 + t4
  
  // compare t1 and t2 in constant time
  r = gfp_cmp(t1, t2, c, len);
  if (r != 0) return MSPECC_ERR_INVALID_POINT;
  
  return MSPECC_NO_ERROR;
}


/*****************************************************************************/
/* Scalar multiplication R = k*P on a twisted Edwards curve according to the */
/* binary method (also referred to as "double and add" method). The result R */
/* is given in extended projective coordinates of the form (X,Y,Z,E,H) with  */
/* E*H = T = X*Y/Z.                                                          */
/*****************************************************************************/

void ted_mul_binary(PROPOINT *r, const UINT16 *k, const AFFPOINT *p,
                    const ECDPARAM *m)
{
  int ki, len = m->len, i = (len<<4)-1;
  UINT16 tmp[3*_len]; // temporary space for three gfp elements
  PROPOINT q = { tmp, &tmp[len], &tmp[2*len], NULL, r->slack };
  
  // find position of first non-zero bit in k
  while ((get_bit(k, i) == 0) && (i >= 0)) i--;
  if (i < 0) { ted_set0_pro(r, len); return; }
  
  // initialize Q with P
  ted_affine_extaff(&q, p, m);
  // initialize R with Q
  ted_extaff_extpro(r, &q, m);
  
  // binary method for scalar multiplication
  for(i = i-1; i >= 0; i--)
  {
    ted_double(r, m);
    ki = get_bit(k, i);
    if (ki) ted_add(r, &q, m);
  }
}


/*****************************************************************************/
/* Conversion of a point P given in standard projective coordinates of the   */
/* form (X,Y,Z) into a point R in standard affine coordinates (x,y). This    */
/* conversion requires the inversion of Z, which can leak information about  */
/* the secret scalar when implemented in a straightforward way using e.g.    */
/* the extended Euclidean algorithm; see the paper "Projective Coordinates   */
/* Leak" (Proceedings of EUROCRYPT 2004). To prevent this leakage, a simple  */
/* masking technique is applied, which means that, instead of inverting Z    */
/* "directly", it is first multiplied by a field element 'u' that is unknown */
/* the attacker, and then the product Z*u is inverted. Finally, the inverse  */
/* (Z*u)^-1 is multiplied by 'u' to obtain Z^-1. This approach for secure    */
/* inversion is described in Section 6 of the paper "SPA Vulnerabilities of  */
/* the Binary Extended Euclidean Algorithm" (Journal of Cryptographic        */
/* Engineering, 2016).                                                       */
/*****************************************************************************/

int ted_proj_affine(PROPOINT *r, const PROPOINT *p, const ECDPARAM *m)
{
  int err, len = m->len; UINT16 c = m->c;
  UINT16 *xp = p->x, *yp = p->y, *zp = p->z;
  UINT16 *xr = r->x, *yr = r->y, *zr = r->z;
  UINT16 *t1 = r->slack, *prod = &(r->slack[len]);
  
  // "masked" inversion of Z to thwart timing attacks
  gfp_mul(t1, zp, INV_MASK, c, len);
  err = gfp_inv(t1, t1, c, len);
  if (err != MSPECC_NO_ERROR) { ted_set0_pro(r, len); return err; }
  gfp_mul(zr, t1, INV_MASK, c, len);
  
  // get least non-negative residue of x = X*(1/Z)
  gfp_mul(t1, xp, zr, c, len);
  gfp_lnr(xr, t1, c, len);
  
  // get least non-negative residue of y = Y*(1/Z)
  gfp_mul(t1, yp, zr, c, len);
  gfp_lnr(yr, t1, c, len);
  
  // overwrite Z with 1
  int_set(zr, 1, len);
  
  return MSPECC_NO_ERROR;
}


/*****************************************************************************/
/* Variable-base scalar multiplication R = k*P on a twisted Edwards curve,   */
/* including some tests to ensure the validity of inputs and outputs. The    */
/* base point P is expected to be given in standard affine coordinates       */
/* (x,y). The result R is also given in standard affine coordinates.         */
/*****************************************************************************/

int ted_mul_varbase(AFFPOINT *r, const UINT16 *k, const AFFPOINT *p,
                    const ECDPARAM *m)
{
  int err, len = m->len;
  UINT16 tmp[8*_len]; // temporary space for eight gfp elements
  PROPOINT q = { tmp, &tmp[len], &tmp[2*len], &tmp[3*len], &tmp[5*len] };
  
  // validate point P (does P satisfy curve equation?)
  ted_aff_to_pro(&q, p, m);
  err = ted_validate(&q, m);
  if (err != MSPECC_NO_ERROR) { ted_set0_aff(r, len); return err; }
  
  // scalar multiplication according to the binary method
  ted_mul_binary(&q, k, p, m);
  
  // convert result from projective to affine coordinates
  err = ted_proj_affine(&q, &q, m);
  if (err != MSPECC_NO_ERROR) { ted_set0_aff(r, len); return err; }
  
  // validate point Q (does Q satisfy curve equation?)
  err = ted_validate(&q, m);
  if (err != MSPECC_NO_ERROR) { ted_set0_aff(r, len); return err; }
  
  // assign x and y-coordinate of point Q to output R
  int_copy(r->x, q.x, len); int_copy(r->y, q.y, len);
  
  return MSPECC_NO_ERROR;
}


/*****************************************************************************/
/* Extract the 4-bit digit di = 8*k[3*maxd+i] + 4*k[2*maxd+i] + 2*k[maxd+i]  */
/* + k[i] from a scalar 'k', where 'maxd' is the number of 4-bit digits the  */
/* scalar consists of, i.e. maxd = 2*len.                                    */
/*****************************************************************************/

int get_digit(const UINT16 *k, int i, int len)
{
  int maxd = (len<<2);
  int d = get_bit(k, i);
  
  i += maxd;
  d += (get_bit(k, i) << 1);
  i += maxd;
  d += (get_bit(k, i) << 2);
  i += maxd;
  d += (get_bit(k, i) << 3);
  
  return d;
}


/*****************************************************************************/
/* Load the point with index 'i' from the pre-computed comb table.           */
/*****************************************************************************/

void ted_load_point(PROPOINT *r, int i, const ECDPARAM *m)
{
  int len = m->len;
  FIXPOINT *f = (FIXPOINT *) &(m->ctab[(i&0xF)]);
  
  int_copy(r->x, f->u, len);
  int_copy(r->y, f->v, len);
  int_copy(r->z, f->w, len);
}


/*****************************************************************************/
/* Scalar multiplication R = k*P on a twisted Edwards curve according to the */
/* fixed-base comb method. Four bits of the scalar 'k' are processed at a    */
/* time, which means the loop is iterated 2*'len'-1 times. The base point P  */
/* is not passed as a parameter but "hard-coded" in the pre-computed comb    */
/* table, which can be accessed via the ECDPARAM structure 'm'. The result R */
/* is given in extended projective coordinates of the form (X,Y,Z,E,H) with  */
/* E*H = T = X*Y/Z.                                                          */
/*****************************************************************************/

void ted_mul_comb4b(PROPOINT *r, const UINT16 *k, const ECDPARAM *m)
{
  int di, len = m->len, i = (len<<2)-1;
  UINT16 tmp[3*_len]; // temporary space for three gfp elements
  PROPOINT q = { tmp, &tmp[len], &tmp[2*len], NULL, r->slack };
  
  // initialize R with di-th element from the pre-computed comb table
  di = get_digit(k, i, len);
  ted_load_point(&q, di, m);
  ted_extaff_extpro(r, &q, m);
  
  // fixed-base comb method (4 bits of k are processed per iteration)
  for(i = i-1; i >= 0; i--)
  {
    ted_double(r, m);
    di = get_digit(k, i, len);
    ted_load_point(&q, di, m);
    ted_add(r, &q, m);
  }
}


/*****************************************************************************/
/* Fixed-base scalar multiplication R = k*P on a twisted Edwards curve,      */
/* including some tests to ensure the validity of inputs and outputs. The    */
/* base point P is not passed as a parameter since it is fixed (normally     */
/* specified in the domain parameters). The result R is given in standard    */
/* affine coordinates.                                                       */
/*****************************************************************************/

int ted_mul_fixbase(AFFPOINT *r, const UINT16 *k, const ECDPARAM *m)
{
  int err, len = m->len;
  UINT16 tmp[8*_len]; // temporary space for eight gfp elements
  PROPOINT q = { tmp, &tmp[len], &tmp[2*len], &tmp[3*len], &tmp[5*len] };
  
  // perform scalar multiplication via fixed-base comb method
  ted_mul_comb4b(&q, k, m);
  
  // convert result from projective to affine coordinates
  err = ted_proj_affine(&q, &q, m);
  if (err != MSPECC_NO_ERROR) { ted_set0_aff(r, len); return err; }
  
  // validate point Q (does Q satisfy curve equation?)
  err = ted_validate(&q, m);
  if (err != MSPECC_NO_ERROR) { ted_set0_aff(r, len); return err; }
  
  // assign x, y-coordinate of affine point to output
  int_copy(r->x, q.x, len); int_copy(r->y, q.y, len);
  
  return MSPECC_NO_ERROR;
}


/*****************************************************************************/
/* Convert a point P on a twisted Edwards curve to the corresponding point R */
/* on the birationally equivalent Montgomery curve. The point P is expected  */
/* to be given in standard projective coordinates and the result R is also   */
/* given in standard projective coordinates. This conversion requires a      */
/* pre-computed constant c = SquareRoot(-(A+2)/B) where A, B are the curve   */
/* parameters of the Montgomery curve.                                       */
/*****************************************************************************/

void ted_to_mon(PROPOINT *r, const PROPOINT *p, const ECDPARAM *m)
{
  int len = m->len; UINT16 c = m->c;
  UINT16 tmp[2*_len]; // temporary space for two gfp elements
  UINT16 *t1 = tmp, *t2 = &tmp[len], *t3 = r->slack, *prod = &(r->slack[len]);
  UINT16 *xt = p->x, *yt = p->y, *zt = p->z;
  UINT16 *xm = r->x, *ym = r->y, *zm = r->z;
  
  gfp_add(t1, zt, yt, c, len);
  gfp_sub(t2, zt, yt, c, len);
  gfp_mul(t3, zt, m->cpc, c, len);
  gfp_mul(ym, t3, t1, c, len);          // ym := c*(zt + yt)*zt;
  gfp_mul(zm, t2, xt, c, len);          // zm := (zt - yt)*xt;
  gfp_mul(xm, t1, xt, c, len);          // xm := (zt + yt)*xt;
}


/*****************************************************************************/
/* Double-base scalar multiplication R = k*P + l*G on a twisted Edwards      */
/* curve, including some tests to ensure the validity of inputs and outputs. */
/* This operation is typically performed in the verification of an ECDSA or  */
/* EdDSA signature and consists of a variable-base scalar multiplication k*P */
/* and a fixed-base scalar multiplication l*G. The base point P for the      */
/* variable-base scalar multiplication is expected to be given in standard   */
/* affine coordinates (x,y). On the other hand, the base point G for the     */
/* fixed-base scalar multiplication is not passed as a parameter since it is */
/* fixed (normally specified as part of the domain parameters). The result R */
/* is given in standard affine coordinates.                                  */
/*****************************************************************************/

int ted_mul_dblbase(AFFPOINT *r, const DBLSCALAR *k, const AFFPOINT *p,
                    const ECDPARAM *m)
{
  int err, len = m->len;
  UINT16 tmp[11*_len]; // temporary space for eleven gfp elements
  PROPOINT q = { tmp, &tmp[len], &tmp[2*len], &tmp[3*len], &tmp[5*len] };
  PROPOINT s = { &tmp[8*len], &tmp[9*len], &tmp[10*len], NULL, q.slack };
  
  // validate point P (does P satisfy curve equation?)
  ted_aff_to_pro(&q, p, m);  
  err = ted_validate(&q, m);
  if (err != MSPECC_NO_ERROR) { ted_set0_aff(r, len); return err; }
  
  // perform variable-base scalar multiplication with k->var
  // ted_mul_binary(&q, k->var, p, m); ted_copy(&s, &q, len, 3);
  ted_to_mon(&s, &q, m);
  mon_proj_affine(&s, &s, m);
  mon_mul_ladder(&q, k->var, s.x, m);
  mon_recover_y(&q, &q, &s, m);
  mon_to_ted(&s, &q, m);
  
  // perform fixed-base scalar multiplication with k->fix
  ted_mul_comb4b(&q, k->fix, m);
  
  // add result of first scalar multiplication to point Q
  ted_add_pro(&q, &s, m);
  
  // convert result from projective to affine coordinates
  err = ted_proj_affine(&q, &q, m);
  if (err != MSPECC_NO_ERROR) { ted_set0_aff(r, len); return err; }
  
  // validate point S (does S satisfy curve equation?)
  err = ted_validate(&q, m);
  if (err != MSPECC_NO_ERROR) { ted_set0_aff(r, len); return err; }
  
  // assign x, y-coordinate of affine point to output
  int_copy(r->x, q.x, len); int_copy(r->y, q.y, len);
  
  return MSPECC_NO_ERROR;
}


/*****************************************************************************/
/* Exponentiation of a field element 'a' by the exponent e = (p-5)/8 where p */
/* is expected to be congruent to 5 mod 8.                                   */
/*****************************************************************************/

void gfp_exp_p58(UINT16 *r, const UINT16 *a, UINT16 c, int len)
{
  int i;
  UINT16 e0 = ((~c)-4)>>3;  // least significant word of e = (p-5)/8
  UINT16 tmp[4*_len];       // temporary space for four gfp elements
  UINT16 *b = tmp, *t = &tmp[len], *prod = &tmp[2*len];
  
  // compute b = a^15 mod p
  gfp_sqr(b, a, c, len);
  gfp_mul(b, b, a, c, len);
  gfp_sqr(b, b, c, len);
  gfp_mul(b, b, a, c, len);
  gfp_sqr(b, b, c, len);
  gfp_mul(b, b, a, c, len);
  
  // the exponent e = (p-5)/8 consists of 'len' words, whereby e[len-1] is
  // always 0x0FFF, e[len-2] to e[1] are all 0xFFFF, and the most significant 4
  // bits of e[0] are 1. Hence, e contains 4*(len-1) 4-bit blocks that are 0xF.
  // the following loop processes these 4-bit blocks; we initialize t with b,
  // which allows us to skip the first 4-bit block
  int_copy(t, b, len);
  for (i = (((len-1)<<2) - 2); i >= 0; i --)
  {
    gfp_sqr(t, t, c, len);
    gfp_sqr(t, t, c, len);
    gfp_sqr(t, t, c, len);
    gfp_sqr(t, t, c, len);
    gfp_mul(t, t, b, c, len);
  }
  
  // process the 12 LSBs of exponent e = (p-5)/8
  for (i = 11; i >= 0; i --)
  {
    gfp_sqr(t, t, c, len);
    if ((e0>>i)&1) gfp_mul(t, t, a, c, len);
  }
  
  int_copy(r, t, len);
}


void ted_test_modexp(void)  // test function for modular exponentiation
{
  UINT16 a[16] = { 0xBA98, 0xFEDC, 0x3210, 0x7654, 0xBA98, 0xFEDC, 0x3210, 0x7654, 0xBA98, 0xFEDC, 0x3210, 0x7654, 0xBA98, 0xFEDC, 0x3210, 0x7654 };
  UINT16 r[16];
  
  #ifdef MSPECC_DEBUG_PRINT
  int_print("a = ", a, 16);
  #endif
  
  gfp_exp_p58(r, a, 19, 16);
  
  #ifdef MSPECC_DEBUG_PRINT
  int_print("r = ", r, 16);
  #endif
}


/*****************************************************************************/
/* Decompression of a point given by its y-coordinate and the LSB of the     */
/* x-coordinate. This function is roughly based on RFC 8032, Section 5.1.3.  */
/*****************************************************************************/

int ted_decompress(AFFPOINT *r, const UINT16 *a, const ECDPARAM *m)
{
  int len = m->len, x0 = (a[len-1]>>15); // x0 is LSB of x-coordinate
  UINT16 c = m->c, tmp[6*_len]; // temporary space for six gfp elements
  UINT16 *x = r->x, *y = r->y, *u = tmp, *v = &tmp[len];
  UINT16 *t1 = &tmp[2*len], *t2 = &tmp[3*len], *prod = &tmp[4*len];
  
  // initialize x, y and t1
  int_set(x, 0, len);
  int_copy(y, a, len);
  y[len-1] &= 0x7FFF;
  int_set(t1, 1, len);
  
  // compute u = y^2 - 1 and v = d*y^2 + 1
  gfp_sqr(u, y, c, len);
  gfp_mul(v, u, m->d, c, len);
  gfp_sub(u, u, t1, c, len);
  gfp_add(v, v, t1, c, len);
  
  // compute t1 = u*v^3 and t2 = t1*v^4 = u*v^7
  gfp_sqr(t1, v, c, len);               // t1 = v^2
  gfp_sqr(t2, t1, c, len);              // t2 = v^4
  gfp_mul(t1, t1, v, c, len);           // t1 = v^3
  gfp_mul(t1, t1, u, c, len);           // t1 = u*v^3
  gfp_mul(t2, t2, t1, c, len);          // t2 = t2*t1 = v^4*u*v^3 = u*v^7
  
  // compute t1 = t1*t2^((p-5)/8) = u*v^3*(u*v^7)^((p-5)/8)
  gfp_exp_p58(t2, t2, c, len);
  gfp_mul(t1, t1, t2, c, len);
  
  // check whether v*(t1)^2 = u; if this is not the case then multiply t1 by
  // root(-1) = 2^((p-1)/4) and check again
  gfp_sqr(t2, t1, c, len);
  gfp_mul(t2, t2, v, c, len);
  if (gfp_cmp(t2, u, c, len) != 0)
  {
    gfp_mul(t1, t1, m->rm1, c, len);
    gfp_sqr(t2, t1, c, len);
    gfp_mul(t2, t2, v, c, len);
    if (gfp_cmp(t2, u, c, len) != 0) return MSPECC_ERR_INVALID_POINT;
  }
  
  // determine whether t1 or -t1 = p-t1 is the correct x-coordinate
  gfp_lnr(x, t1, c, len);  // get least non-negative residue of t1
  if (int_is0(x, len) && (x0 != 0)) return MSPECC_ERR_INVALID_POINT;
  if ((x[0]&1) != x0) { gfp_set(t1, c, len); int_sub(x, t1, x, len); }

  return MSPECC_NO_ERROR;
}
