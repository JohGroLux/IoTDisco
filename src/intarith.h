#ifndef _INTARITH_H
#define _INTARITH_H

#include "typedefs.h"

/* prototypes of functions for which only C implementations exist (no ASM) */
int  int_is0(const UINT16 *a, int len);
int  int_is1(const UINT16 *a, int len);
void int_set(UINT16 *r, UINT16 a, int len);
int  int_cmp(const UINT16 *a, const UINT16 *b, int len);  // comparison
void int_copy(UINT16 *r, const UINT16 *a, int len);
void int_print(const char *c, const UINT16 *a, int len);  // print in hex

/* prototypes of functions for which both C and ASM implementations exist  */
int  int_shr_c99(UINT16 *r, const UINT16 *a, int len);
int  int_add_c99(UINT16 *r, const UINT16 *a, const UINT16 *b, int len);
int  int_sub_c99(UINT16 *r, const UINT16 *a, const UINT16 *b, int len);
void int_mul_c99(UINT16 *r, const UINT16 *a, const UINT16 *b, int len);
void int_mul32_c99(UINT16 *r, const UINT16 *a, const UINT16 *b, int len);
void int_sqr_c99(UINT16 *r, const UINT16 *a, int len);

#endif
