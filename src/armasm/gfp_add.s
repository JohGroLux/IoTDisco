///////////////////////////////////////////////////////////////////////////////
// gfp_add.s: Multiple-Precision Addition Modulo a Pseudo-Mersenne Prime.    //
// This file is part of SECCCM3, a Scalable ECC implementation for Cortex-M3 //
// Version 0.9.0 (16-08-24), see <http://github.com/johgrolux/> for updates. //
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


    AREA gfparith, CODE, READONLY ; Name this block of code ARMex
    
    EXPORT gfp_add_asm
    ALIGN 2
    
    
;;;;;;;;;;;;;;;;;;;;
;; Register Names ;;
;;;;;;;;;;;;;;;;;;;;
    
rPtr     RN r0
aPtr     RN r1
bPtr     RN r2
cWord    RN r3
Len      RN r4
opWord   RN r5
sumLo    RN r6
sumHi    RN r7
msWord   RN r8
    
    
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Scalable Modular Addition ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    
gfp_add_asm FUNCTION
    
    push  {r4-r8}
    ldr   Len, [sp, #20]  
    sub   Len, Len, #1
    
    ldr   sumLo, [aPtr, Len, LSL #2]
    ldr   opWord, [bPtr, Len, LSL #2]
    adds  sumLo, sumLo, opWord
    and   msWord, sumLo, #0x7fffffff
    lsr   sumLo, sumLo, #31 ; right-shift sumLo 31 bits (C flag not affected!)
    adc   sumLo, sumLo, #0  ; add C flag generated by adds sumLo, sumLo, opWord
    adc   sumLo, sumLo, #0  ; add C flag generated by adds sumLo, sumLo, opWord
    umull sumLo, sumHi, cWord, sumLo  ; sumLo contains lower part of product
    
loop
    ldr   opWord, [aPtr], #4
    adds  sumLo, sumLo, opWord
    adc   sumHi, sumHi, #0
    ldr   opWord, [bPtr], #4
    adds  sumLo, sumLo, opWord
    str   sumLo, [rPtr], #4
    adc   sumLo, sumHi, #0
    mov   sumHi, #0
    subs  Len, Len, #1
    bne   loop
    
    add   sumLo, sumLo, msWord
    str   sumLo, [rPtr]
    
    pop   {r4-r8}
    bx    lr
    
    ENDFUNC
    
    
    END
