#ifndef _CONFIG_H
#define _CONFIG_H

#define MSPECC_MIN_LEN 10 

#define MSPECC_MAX_LEN 16

// 
#define MSPECC_USE_ASM

// define MSPECC_USE_VLA to use Variable-Length Arrays (VLA)
// undefine it to use static arrays of length MSPECC_MAX_LEN
#define MSPECC_USE_VLA

#ifndef NDEBUG
#define MSPECC_DEBUG_PRINT
#endif

#define MSPECC_NO_ERROR           0
#define MSPECC_ERR_INVERSION_ZERO 1
#define MSPECC_ERR_INVALID_POINT  2
#define MSPECC_ERR_INVALID_SCALAR 4

#ifdef MSPECC_USE_ASM
#include "asmfncts.h"
#define int_add(r, a, b, len) int_add_msp((r), (a), (b), (len))
#define int_mul(r, a, b, len) int_mul_msp((r), (a), (b), (len))
#define int_mul32(r, a, b, len) int_mul32_msp((r), (a), (b), (len))
#define int_shr(r, a, len) int_shr_msp((r), (a), (len))
#define int_sqr(r, a, len) int_sqr_msp((r), (a), (len))
#define int_sub(r, a, b, len) int_sub_msp((r), (a), (b), (len))
#define gfp_add(r, a, b, c, len) gfp_add_msp((r), (a), (b), (c), (len))
#define gfp_cneg(r, a, neg, c, len) gfp_cneg_msp((r), (a), (neg), (c), (len))
#define gfp_hlv(r, a, c, len) gfp_hlv_msp((r), (a), (c), (len))
#define gfp_red(r, a, c, len) gfp_red_msp((r), (a), (c), (len))
#define gfp_red32(r, a, c, len) gfp_red32_msp((r), (a), (c), (len))
#define gfp_sub(r, a, b, c, len) gfp_sub_msp((r), (a), (b), (c), (len))
#else
#define int_add(r, a, b, len) int_add_c99((r), (a), (b), (len))
#define int_mul(r, a, b, len) int_mul_c99((r), (a), (b), (len))
#define int_mul32(r, a, b, len) int_mul32_c99((r), (a), (b), (len))
#define int_shr(r, a, len) int_shr_c99((r), (a), (len))
#define int_sqr(r, a, len) int_sqr_c99((r), (a), (len))
#define int_sub(r, a, b, len) int_sub_c99((r), (a), (b), (len))
#define gfp_add(r, a, b, c, len) gfp_add_c99((r), (a), (b), (c), (len))
#define gfp_cneg(r, a, neg, c, len) gfp_cneg_c99((r), (a), (neg), (c), (len))
#define gfp_hlv(r, a, c, len) gfp_hlv_c99((r), (a), (c), (len))
#define gfp_red(r, a, c, len) gfp_red_c99((r), (a), (c), (len))
#define gfp_red32(r, a, c, len) gfp_red32_c99((r), (a), (c), (len))
#define gfp_sub(r, a, b, c, len) gfp_sub_c99((r), (a), (b), (c), (len))
#endif  // MSPECC_USE_ASM

#endif  // _CONFIG_H
