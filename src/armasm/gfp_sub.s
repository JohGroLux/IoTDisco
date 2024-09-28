;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; gfp_sub.s: Multiple-Precision Subtraction Modulo a Pseudo-Mersenne Prime. ;;
;; This file is part of SECCCM3, a Scalable ECC implementation for Cortex-M3 ;;
;; Version 0.9.0 (16-08-24), see <http:;;github.com/johgrolux/> for updates. ;;
;; License: GPLv3 (see LICENSE file), other licenses available upon request. ;;
;; ------------------------------------------------------------------------- ;;
;; This program is free software: you can redistribute it and/or modify it   ;;
;; under the terms of the GNU General Public License as published by the     ;;
;; Free Software Foundation, either version 3 of the License, or (at your    ;;
;; option) any later version. This program is distributed in the hope that   ;;
;; it will be useful, but WITHOUT ANY WARRANTY; without even the implied     ;;
;; warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the  ;;
;; GNU General Public License for more details. You should have received a   ;;
;; copy of the GNU General Public License along with this program. If not,   ;;
;; see <http:;;www.gnu.org/licenses/>.                                       ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


    AREA gfparith, CODE, READONLY ; Name this block of code ARMex
    
    EXPORT gfp_sub_asm
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
    
    
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Generic Modular Subtraction ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    
gfp_sub_asm PROC
    
    push  {r4-r8}
    ldr   Len, [sp, #20]  
    sub   Len, Len, #1
    
    mov   sumLo, #0xfffffffc
    mov   sumHi, #1
    ldr   opWord, [aPtr, Len, LSL #2]
    adds  sumLo, sumLo, opWord
    adc   sumHi, sumHi, #0
    ldr   opWord, [bPtr, Len, LSL #2]
    subs  sumLo, sumLo, opWord
    sbc   sumHi, sumHi, #0
    and   msWord, sumLo, #0x7fffffff
    lsr   sumLo, sumLo, #31
    orr   sumLo, sumLo, sumHi, LSL #1
    umull sumLo, sumHi, cWord, sumLo  ; sumLo contains lower part of product
    subs  sumLo, sumLo, cWord, LSL #2
    sbc   sumHi, sumHi, #0
    
loop
    ldr   opWord, [aPtr], #4
    adds  sumLo, sumLo, opWord
    adc   sumHi, sumHi, #0
    ldr   opWord, [bPtr], #4
    subs  sumLo, sumLo, opWord
    str   sumLo, [rPtr], #4
    sbc   sumLo, sumHi, #0
    asr   sumHi, sumHi, #32
    subs  Len, Len, #1
    bne   loop
    
    add   sumLo, sumLo, msWord
    add   sumLo, sumLo, #4
    str   sumLo, [rPtr]
    
    pop   {r4-r8}
    bx    lr
    
    ENDP
    
    
    END
