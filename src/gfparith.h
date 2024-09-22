#ifndef _GFPARITH_H
#define _GFPARITH_H

#include "config.h"

/* prototypes of functions for which only C implementations exist (no ASM) */
void gfp_set(Word *r, Word c, int len);
int  gfp_isp(const Word *a, Word c, int len);

/* prototypes of functions for which both C and ASM implementations exist  */
void gfp_add_c99(Word *r, const Word *a, const Word *b, Word c, int len);
void gfp_sub_c99(Word *r, const Word *a, const Word *b, Word c, int len);
void gfp_cneg_c99(Word *r, const Word *a, Word c, int neg, int len);
void gfp_hlv_c99(Word *r, const Word *a, Word c, int len);
void gfp_red_c99(Word *r, const Word *a, Word c, int len);
void gfp_mul_c99(Word *r, const Word *a, const Word *b, Word c, int len);
void gfp_sqr_c99(Word *r, const Word *a, Word c, int len);
void gfp_red32_c99(Word *r, const Word *a, Word c, int len);
void gfp_mul32_c99(Word *r, const Word *a, const Word *b, Word c, int len);

/* prototypes of functions for which only C implementations exist, but they */
/* contain sub-functions with C and ASM implementations                     */
void gfp_lnr(Word *r, const Word *a, Word c, int len);
int  gfp_cmp(Word *a, Word *b, Word c, int len);
int  gfp_inv(Word *r, const Word *a, Word c, int len);

#endif
