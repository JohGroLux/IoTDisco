#ifndef _GFPARITH_H
#define _GFPARITH_H

#include "typedefs.h"

/* prototypes of functions for which only C implementations exist (no ASM) */
void gfp_set(UINT16 *r, UINT16 c, int len);
int  gfp_isp(const UINT16 *a, UINT16 c, int len);

/* prototypes of functions for which both C and ASM implementations exist  */
void gfp_add_c99(UINT16 *r, const UINT16 *a, const UINT16 *b, UINT16 c, int len);
void gfp_sub_c99(UINT16 *r, const UINT16 *a, const UINT16 *b, UINT16 c, int len);
void gfp_cneg_c99(UINT16 *r, const UINT16 *a, UINT16 c, int neg, int len);
void gfp_hlv_c99(UINT16 *r, const UINT16 *a, UINT16 c, int len);
void gfp_red_c99(UINT16 *r, const UINT16 *a, UINT16 c, int len);
void gfp_red32_c99(UINT16 *r, const UINT16 *a, UINT16 c, int len);

/* prototypes of functions for which only C implementations exist, but they */
/* contain sub-functions with C and ASM implementations                     */
void gfp_lnr(UINT16 *r, const UINT16 *a, UINT16 c, int len);
int  gfp_cmp(UINT16 *a, UINT16 *b, UINT16 c, int len);
int  gfp_inv(UINT16 *r, const UINT16 *a, UINT16 c, int len);

#endif
