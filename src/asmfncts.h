#ifndef _ASMFNCTS_H
#define _ASMFNCTS_H

#include "typedefs.h"

/*------Assembler function prototypes------*/
extern int  int_add_msp(UINT16 *r, const UINT16 *a, const UINT16 *b, int len);
extern void int_mul_msp(UINT16 *r, const UINT16 *a, const UINT16 *b, int len);
extern void int_mul32_msp(UINT16 *r, const UINT16 *a, const UINT16 *b, int len);
extern int  int_shr_msp(UINT16 *r, const UINT16 *a, int len);
extern void int_sqr_msp(UINT16 *r, const UINT16 *a, int len);
extern int  int_sub_msp(UINT16 *r, const UINT16 *a, const UINT16 *b, int len);

extern void gfp_add_msp(UINT16 *r, const UINT16 *a, const UINT16 *b, UINT16 c, int len);
extern void gfp_cneg_msp(UINT16 *r, const UINT16 *a, int neg, UINT16 c, int len);
extern void gfp_hlv_msp(UINT16 *r, const UINT16 *a, UINT16 c, int len);
extern void gfp_red_msp(UINT16 *r, const UINT16 *a, UINT16 c, int len);
extern void gfp_red32_msp(UINT16 *r, const UINT16 *a, UINT16 c, int len);
extern void gfp_sub_msp(UINT16 *r, const UINT16 *a, const UINT16 *b, UINT16 c, int len);

#endif
