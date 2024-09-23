#ifndef _TYPEDEFS_H
#define _TYPEDEFS_H

#include "config.h"

/*******************************************************/
/* definition of data type for an elliptic curve point */
/*******************************************************/

typedef struct pro_point // curve point in (extended) projective coordinates
{
  Word *x;     // x coordinate
  Word *y;     // y coordinate (or slack space if not used)
  Word *z;     // z coordinate
  Word *extra; // extra space (for storing additional coordinates)
  Word *slack; // slack space (can be shared among several points)
} PROPOINT;

/*************************************************************/
/* definition of data type for a point in affine coordinates */
/*************************************************************/

typedef struct aff_point // curve point in affine coordinates
{
  Word *x;
  Word *y;
} AFFPOINT;

/*************************************************************/
/* definition of data type for a fixed point, e.g. generator */
/*************************************************************/

typedef struct fix_point // curve point in extended affine coordinates
{
  const Word *u;
  const Word *v;
  const Word *w;
} FIXPOINT;

/************************************************************/
/* definition of data type for double scalar multiplication */
/************************************************************/

typedef struct dbl_scalar // two scalars, each consisting of len words
{
  Word *fix;   // scalar for fixed-base scalar multiplication
  Word *var;   // scalar for variable-base scalar multiplication
} DBLSCALAR;


/*****************************************************/
/* definition of data type for ECC domain parameters */
/*****************************************************/

typedef struct ecparam {  // elliptic curve domain parameters
  const int len;    // number of w-bit words the prime p consists of
  const Word c;     // constant c defining the prime p = 2^(w*len-1) - c
  const Word *a24;  // constant (A+2)/4 of Montgomery curve(is usually small)
  const Word *dte;  // parameter d of corresponding TED curve with a = -1
  const Word *rma;  // root of -a = -(A+2)/B (for point-conversion MON <-> TED)
  const Word *rm1;  // root of -1 (i.e., 2^((p-1)/4) mod p for decompression)
  const Word *tbl;  // table of pre-computed points for fixed-base comb method
} ECDPARAM;

#endif
