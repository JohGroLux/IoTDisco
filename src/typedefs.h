#ifndef _TYPEDEFS_H
#define _TYPEDEFS_H

/******************************************/
/* definitions for standard integer types */
/******************************************/

#ifdef _MSC_VER
typedef __int8 INT8;              // 8-bit integer (-128~127)
typedef unsigned __int8 UINT8;    // 8-bit unsigned integer (0~255)
typedef __int16 INT16;            // 16-bit integer
typedef unsigned __int16 UINT16;  // 16-bit unsigned integer
typedef __int32 INT32;            // 32-bit integer
typedef unsigned __int32 UINT32;  // 32-bit unsigned integer
typedef __int64 INT64;            // 64-bit integer
typedef unsigned __int64 UINT64;  // 64-bit unsigned integer
#else
#include <stdint.h>
typedef int8_t INT8;
typedef uint8_t UINT8;
typedef int16_t INT16;
typedef uint16_t UINT16;
typedef int32_t INT32;
typedef uint32_t UINT32;
typedef int64_t INT64;
typedef uint64_t UINT64;
#endif /* _MSC_VER */

/*******************************************************/
/* definition of data type for an elliptic curve point */
/*******************************************************/

typedef struct pro_point // curve point in (extended) projective coordinates
{
  UINT16 *x;     // x coordinate
  UINT16 *y;     // y coordinate (or slack space if not used)
  UINT16 *z;     // z coordinate
  UINT16 *extra; // extra space (for storing additional coordinates)
  UINT16 *slack; // slack space (can be shared among several points)
} PROPOINT;

/*************************************************************/
/* definition of data type for a point in affine coordinates */
/*************************************************************/

typedef struct aff_point // curve point in affine coordinates
{
  UINT16 *x;
  UINT16 *y;
} AFFPOINT;

/*************************************************************/
/* definition of data type for a fixed point, e.g. generator */
/*************************************************************/

typedef struct fix_point // curve point in extended affine coordinates
{
  const UINT16 *u;
  const UINT16 *v;
  const UINT16 *w;
} FIXPOINT;

/************************************************************/
/* definition of data type for double scalar multiplication */
/************************************************************/

typedef struct dbl_scalar // two scalars, each consisting of len words
{
  UINT16 *fix;   // scalar for fixed-base scalar multiplication
  UINT16 *var;   // scalar for variable-base scalar multiplication
} DBLSCALAR;


/*****************************************************/
/* definition of data type for ECC domain parameters */
/*****************************************************/

typedef struct ecd_param // elliptic curve domain parameters
{
  const int len;         // number of 16-words the prime p consists of
  const UINT16 c;        // constant c defining the prime p = 2^(16*len-1) - c
  const UINT16 *a24;     // constant (A+2)/4 of Montgomery curve
  const UINT16 *d;       // parameter d of twisted Edwards curve
  const UINT16 *rm1;     // root of -1 mod p, i.e. 2^((p-1)/4) mod p
  const UINT16 *cpc;     // constant for point conversion between MON and TED
  const FIXPOINT *ctab;  // pre-computed table for fixed-base comb method
} ECDPARAM;

#endif
