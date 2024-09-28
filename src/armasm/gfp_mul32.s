;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; gfp_mul32.s: Multiplication by 32-bit Int Modulo a Pseudo-Mersenne Prime. ;;
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


    AREA gfparith, CODE, READONLY
    
    EXPORT gfp_mul32_asm
    ALIGN 2
    
    
;;;;;;;;;;;;;;;;;;;;
;; Register Names ;;
;;;;;;;;;;;;;;;;;;;;
    
rPtr     RN r0
aPtr     RN r1
bPtr     RN r2
bWord    RN r2
cWord    RN r3
Len      RN r4
mulOp1   RN r5
prodLo   RN r5
prodHi   RN r6
hiWord   RN r7
lStop    RN r8
    
mask     RN r2
    
    
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Generic Modular Multiplication by a 32-bit Integer ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    
gfp_mul32_asm PROC
    
    push  {r4-r8}       ; push callee-saved registers
    ldr   Len, [sp, #20]
    
    ;; Initialize pointers and variables
    
    ldr   bWord, [bPtr]  ; we only use the very first word b[0] of array b
    add   lStop, aPtr, Len, LSL #2
    sub   lStop, lStop, #4
    mov   hiWord, #0
    
    ;; START OF FIRST OUTER LOOP
    
OutLoop1
    
    ldr   mulOp1, [aPtr], #4
    umull prodLo, prodHi, mulOp1, bWord
    adds  prodLo, prodLo, hiWord
    adc   prodHi, prodHi, #0
    str   prodLo, [rPtr], #4
    mov   hiWord, prodHi
    cmp   aPtr, lStop
    blt   OutLoop1   
    
    ;; END OF FIRST OUTER LOOP
    
    ldr   mulOp1, [aPtr], #4
    umull prodLo, prodHi, mulOp1, bWord
    adds  prodLo, prodLo, hiWord
    adc   prodHi, prodHi, #0
    
    ;; start the reduction
    
    adds  prodLo, prodLo, prodLo
    adcs  prodHi, prodHi, prodHi
    sbc   mask, mask, mask   ; mask is 0 is C was 1 and 0xFFFFFFFF otherwise
    bic   mask, cWord, mask  ; mask is now either cWord (if C was 1) or 0
    lsr   prodLo, prodLo, #1
    str   prodLo, [rPtr], #4
    
    mov   lStop, rPtr
    sub   rPtr, rPtr, Len, LSL #2
    
    umull prodLo, prodHi, prodHi, cWord
    add   prodHi, prodHi, mask
    ldr   hiWord, [rPtr]
    adds  prodLo, prodLo, hiWord
    str   prodLo, [rPtr], #4
    ldr   hiWord, [rPtr]
    adcs  prodHi, prodHi, hiWord
    str   prodHi, [rPtr], #4
    
    ;; START OF SECOND OUTER LOOP
    
OutLoop2
    ldr   hiWord, [rPtr]
    adcs  hiWord, hiWord, #0
    str   hiWord, [rPtr], #4
    teq   rPtr, lStop
    bne   OutLoop2
    
    ;; END OF SECOND OUTER LOOP
    
    pop   {r4-r8}
    bx    lr
    
    ENDP
    
    
    END
