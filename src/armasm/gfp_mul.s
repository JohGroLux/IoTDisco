;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; gfp_mul.s: Multi-Precision Multiplication Modulo a Pseudo-Mersenne Prime. ;;
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
    
    EXPORT gfp_mul_asm
    ALIGN 2
    
    
;;;;;;;;;;;;;;;;;;;;
;; Register Names ;;
;;;;;;;;;;;;;;;;;;;;
    
rPtr     RN r0
aPtr     RN r1
bPtr     RN r2
Len      RN r3
mulOp1   RN r4
mulOp2   RN r5
prodLo   RN r4
prodHi   RN r5
accuLo   RN r6
accuHi   RN r7
accuEx   RN r8
innStop  RN r9
outStop  RN r10
aStart   RN r11
bStart   RN r12
    
cWord    RN r9
loWord   RN r11
hiWord   RN r12
    
    
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Generic Modular Multiplication ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    
gfp_mul_asm PROC
    
    push  {r0, r3-r12}       ; push callee-saved registers
    
    ldr   Len, [sp, #44]
    sub   sp, sp, Len, LSL #3
    mov   rPtr, sp
    
    ;; Initialize pointers and variables
    
    mov   aStart, aPtr       ;
    mov   bStart, bPtr       ;
    mov   innStop, bPtr      ;
    add   outStop, rPtr, Len, LSL #2 ;
    
    ;; Do the very first multiplication a[0]*b[0] and store LSW in r[0] and the
    ;; MSW in the accu registers
    
    ldr   mulOp1, [aPtr]
    ldr   mulOp2, [bPtr]
    umull mulOp2, accuLo, mulOp1, mulOp2
    str   mulOp2, [rPtr], #4
    mov   accuHi, #0
    mov   accuEx, #0
    
    ;; START OF FIRST OUTER LOOP
    
OutLoop1
    
    add   bStart, bStart, #4
    mov   aPtr, aStart
    mov   bPtr, bStart
    
InnLoop1
    ldr   mulOp1, [aPtr], #4
    ldr   mulOp2, [bPtr], #-4
    umull prodLo, prodHi, mulOp1, mulOp2
    adds  accuLo, accuLo, prodLo
    adcs  accuHi, accuHi, prodHi
    adc   accuEx, accuEx, #0
    cmp   bPtr, innStop
    bge   InnLoop1
    
    str   accuLo, [rPtr], #4
    mov   accuLo, accuHi
    mov   accuHi, accuEx
    mov   accuEx, #0
    
    cmp   rPtr, outStop
    blt   OutLoop1
    
    ;; END OF FIRST OUTER LOOP
    
    add   innStop, aStart, Len, LSL #2
    add   outStop, outStop, Len, LSL #2
    sub   outStop, outStop, #8
    
    ;; START OF SECOND OUTER LOOP
    
OutLoop2
    
    add   aStart, aStart, #4
    mov   aPtr, aStart
    mov   bPtr, bStart
    
InnLoop2
    ldr   mulOp1, [aPtr], #4
    ldr   mulOp2, [bPtr], #-4
    umull prodLo, prodHi, mulOp1, mulOp2
    adds  accuLo, accuLo, prodLo
    adcs  accuHi, accuHi, prodHi
    adc   accuEx, accuEx, #0
    cmp   aPtr, innStop
    blt   InnLoop2
    
    str   accuLo, [rPtr], #4
    mov   accuLo, accuHi
    mov   accuHi, accuEx
    mov   accuEx, #0
    
    cmp   rPtr, outStop
    blt   OutLoop2
    
    ;; END OF SECOND OUTER LOOP
    
    ;; Do the very last multiplication a[len-1]*b[len-1]
    
    ldr   mulOp1, [aPtr, #-4]
    ldr   mulOp2, [bPtr, #8]
    umull prodLo, prodHi, mulOp1, mulOp2
    adds  accuLo, accuLo, prodLo
    adc   accuHi, accuHi, prodHi
    stmia rPtr!, {accuLo, accuHi}
    
    ;; Start of the reduction operation
    
    sub   aPtr, rPtr, Len, LSL #3
    sub   bPtr, rPtr, Len, LSL #2
    sub   outStop, rPtr, #4
    add   rPtr, sp, Len, LSL #3
    ldr   cWord, [rPtr, #4]
    ldr   rPtr, [rPtr]
    lsl   cWord, cWord, #1
    mov   hiWord, #0
    
    ;; START OF THIRD OUTER LOOP
    
OutLoop3
    ldr   mulOp1, [bPtr], #4
    umull prodLo, prodHi, mulOp1, cWord
    ldr   loWord, [aPtr], #4
    adds  prodLo, prodLo, loWord
    adc   prodHi, prodHi, #0
    adds  prodLo, prodLo, hiWord
    adc   prodHi, prodHi, #0
    str   prodLo, [rPtr], #4
    mov   hiWord, prodHi
    cmp   bPtr, outStop
    blt   OutLoop3
    
    ;; END OF THIRD OUTER LOOP
    
    ldr   mulOp1, [bPtr], #4
    umull prodLo, prodHi, mulOp1, cWord
    ldr   loWord, [aPtr], #4
    adds  prodLo, prodLo, loWord
    adc   prodHi, prodHi, #0
    adds  prodLo, prodLo, hiWord
    adc   prodHi, prodHi, #0
    lsls  prodLo, prodLo, #1
    adc   prodHi, prodHi, prodHi
    lsr   prodLo, prodLo, #1
    str   prodLo, [rPtr], #4
    
    mov   outStop, rPtr
    sub   rPtr, rPtr, Len, LSL #2
    lsr   cWord, cWord, #1
    umull prodLo, prodHi, prodHi, cWord
    ldmia rPtr, {loWord, hiWord}
    adds  prodLo, prodLo, loWord
    adcs  prodHi, prodHi, hiWord
    stmia rPtr!, {prodLo, prodHi}
    
    ;; START OF FOURTH OUTER LOOP
    
OutLoop4
    ldr   loWord, [rPtr]
    adcs  loWord, loWord, #0
    str   loWord, [rPtr], #4
    teq   rPtr, outStop
    bne   OutLoop4
    
    ;; END OF FOURTH OUTER LOOP
    
    add   sp, sp, Len, LSL #3
    pop   {r0, r3-r12}
    bx    lr
    
    ENDP
    
    
    END
