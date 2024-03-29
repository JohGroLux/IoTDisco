///////////////////////////////////////////////////////////////////////////////
// mspcfg.h: Macros to Push/Pop Callee-Saved Registers on MSP430(X) Devices. //
// This file is part of SECC430, a Scalable ECC implementation for MSP430.   //
// Version 1.1.0 (2023-11-13), see <http://www.cryptolux.org/> for updates.  //
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


#ifndef _MSPCFG_H
#define _MSPCFG_H

// choose the target device
#include <msp430f1611.h>
// #include <msp430fr5969.h>
// include <msp430f2617.h>
// #include <msp430f5438a.h>

/*
// The Assembler functions of SECC430 are supposed to work correctly with any
// combination of code model and data model. This requires __CODE_MODEL__ and
// __DATA_MODEL__ to be defined, which is always the case for IAR Assembler
// version >= 6.40.1. However, when using an earlier version of IAR Assembler,
// these definitions must be made manually as shown below for large code model
// and medium data model. Of course, the definitions must match the settings
// for code model and data model in Project -> Options -> General Options.
#define __CODE_MODEL_SMALL__ 0
#define __CODE_MODEL_LARGE__ 1
#define __DATA_MODEL_SMALL__ 0
#define __DATA_MODEL_MEDIUM__ 1
#define __DATA_MODEL_LARGE__ 2
#define __CODE_MODEL__ __CODE_MODEL_LARGE__
#define __DATA_MODEL__ __DATA_MODEL_MEDIUM__
*/

#if defined(__MSP430_HAS_MSP430X_CPU__)||defined(__MSP430_HAS_MSP430XV2_CPU__)

// When using the medium and large data model, all 20 bits of the callee-saved 
// registers need to be preserved, which means the registers have to be pushed
// via PUSHX.A or PUSHM.A and popped via POPX.A or POPM.A. Each pushed 20-bit
// register occupies four bytes on the stack.
#if defined(__DATA_MODEL__)
#if (__DATA_MODEL__ == __DATA_MODEL_SMALL__)
#define PUSHX_ PUSHX.W
#define POPX_ POPX.W
#define PUSHM_ PUSHM.W
#define POPM_ POPM.W
#define BYTES_PER_PUSH 2
#else  // Data Model is Medium or Large
#define PUSHX_ PUSHX.A
#define POPX_ POPX.A
#define PUSHM_ PUSHM.A
#define POPM_ POPM.A
#define BYTES_PER_PUSH 4
#endif
#else  // __DATA_MODEL__ not defined
#error __DATA_MODEL__ not defined (IAR Assembler version >= 6.40.1 required!)
#endif

#if defined(__CODE_MODEL__)
#if (__CODE_MODEL__ == __CODE_MODEL_SMALL__)
#define OFFSET(n) (BYTES_PER_PUSH*(n) + 2)
#define RET_ RET
#else  // Code Model is Large
#define OFFSET(n) (BYTES_PER_PUSH*(n) + 4)
#define RET_ RETA
#endif
#else  // __CODE_MODEL__ not defined
#error __CODE_MODEL__ not defined (IAR Assembler version >= 6.40.1 required!)
#endif

PUSH_R11 MACRO
    PUSHX_ R11
    ENDM

POP_R11 MACRO
    POPX_  R11
    ENDM

PUSH_R10_TO_R11 MACRO
    PUSHM_ #2, R11
    ENDM

POP_R11_TO_R10 MACRO
    POPM_  #2, R11
    ENDM

PUSH_R9_TO_R11 MACRO
    PUSHM_ #3, R11
    ENDM

POP_R11_TO_R9 MACRO
    POPM_  #3, R11
    ENDM

PUSH_R8_TO_R11 MACRO
    PUSHM_ #4, R11
    ENDM

POP_R11_TO_R8 MACRO
    POPM_  #4, R11
    ENDM

PUSH_R7_TO_R11 MACRO
    PUSHM_ #5, R11
    ENDM

POP_R11_TO_R7 MACRO
    POPM_  #5, R11
    ENDM

PUSH_R6_TO_R11 MACRO
    PUSHM_ #6, R11
    ENDM

POP_R11_TO_R6 MACRO
    POPM_  #6, R11
    ENDM

PUSH_R5_TO_R11 MACRO
    PUSHM_ #7, R11
    ENDM

POP_R11_TO_R5 MACRO
    POPM_  #7, R11
    ENDM

PUSH_R4_TO_R11 MACRO
    PUSHM_ #8, R11
    ENDM

POP_R11_TO_R4 MACRO
    POPM_  #8, R11
    ENDM

EINT_ MACRO
    NOP                     ; required to get rid of a warning in IAR ASM 6.50
    EINT                    ; enable interrupts
    NOP                     ; required to get rid of a warning in IAR ASM 6.50
    ENDM

#else  /* we have a conventional MSP430 device, no MSP430X or MSP430XV2 !!! */

#define OFFSET(n) (2*(n) + 2)

#define RET_ RET

PUSH_R11 MACRO
    PUSH R11
    ENDM

POP_R11 MACRO
    POP  R11
    ENDM

PUSH_R10_TO_R11 MACRO
    PUSH R10
    PUSH R11
    ENDM

POP_R11_TO_R10 MACRO
    POP  R11
    POP  R10
    ENDM

PUSH_R9_TO_R11 MACRO
    PUSH R9
    PUSH R10
    PUSH R11
    ENDM

POP_R11_TO_R9 MACRO
    POP  R11
    POP  R10
    POP  R9
    ENDM

PUSH_R8_TO_R11 MACRO
    PUSH R8
    PUSH R9
    PUSH R10
    PUSH R11
    ENDM

POP_R11_TO_R8 MACRO
    POP  R11
    POP  R10
    POP  R9
    POP  R8
    ENDM

PUSH_R7_TO_R11 MACRO
    PUSH R7
    PUSH R8
    PUSH R9
    PUSH R10
    PUSH R11
    ENDM

POP_R11_TO_R7 MACRO
    POP  R11
    POP  R10
    POP  R9
    POP  R8
    POP  R7
    ENDM

PUSH_R6_TO_R11 MACRO
    PUSH R6
    PUSH R7
    PUSH R8
    PUSH R9
    PUSH R10
    PUSH R11
    ENDM

POP_R11_TO_R6 MACRO
    POP  R11
    POP  R10
    POP  R9
    POP  R8
    POP  R7
    POP  R6
    ENDM

PUSH_R5_TO_R11 MACRO
    PUSH R5
    PUSH R6
    PUSH R7
    PUSH R8
    PUSH R9
    PUSH R10
    PUSH R11
    ENDM

POP_R11_TO_R5 MACRO
    POP  R11
    POP  R10
    POP  R9
    POP  R8
    POP  R7
    POP  R6
    POP  R5
    ENDM

PUSH_R4_TO_R11 MACRO
    PUSH R4
    PUSH R5
    PUSH R6
    PUSH R7
    PUSH R8
    PUSH R9
    PUSH R10
    PUSH R11
    ENDM

POP_R11_TO_R4 MACRO
    POP  R11
    POP  R10
    POP  R9
    POP  R8
    POP  R7
    POP  R6
    POP  R5
    POP  R4
    ENDM

EINT_ MACRO
    EINT
    ENDM

#endif

#endif  /* _MSPCFG_H */
