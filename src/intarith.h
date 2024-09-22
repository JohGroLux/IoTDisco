#ifndef _INTARITH_H
#define _INTARITH_H

#include "config.h"

/* prototypes of functions for which only C implementations exist (no ASM) */
int  int_is0(const Word *a, int len);
int  int_is1(const Word *a, int len);
void int_set(Word *r, Word a, int len);
int  int_cmp(const Word *a, const Word *b, int len);  // comparison
void int_copy(Word *r, const Word *a, int len);
void int_print(const char *c, const Word *a, int len);  // print in hex

/* prototypes of functions for which both C and ASM implementations exist  */
int  int_shr_c99(Word *r, const Word *a, int len);
int  int_add_c99(Word *r, const Word *a, const Word *b, int len);
int  int_sub_c99(Word *r, const Word *a, const Word *b, int len);
void int_mul_c99(Word *r, const Word *a, const Word *b, int len);
void int_mul32_c99(Word *r, const Word *a, const Word *b, int len);
void int_sqr_c99(Word *r, const Word *a, int len);

#endif
