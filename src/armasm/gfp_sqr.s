///////////////////////////////////////////////////////////////////////////////
// gfp_sqr.s: Multiple-Precision Squaring Modulo a Pseudo-Mersenne Prime.    //
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


    AREA gfparith, CODE, READONLY
    
    EXPORT gfp_sqr_asm
    ALIGN 2
    
    
;;;;;;;;;;;;;;;;;;;;
;; Register Names ;;
;;;;;;;;;;;;;;;;;;;;
    
rPtr     RN r0
aPtr     RN r1
cWord    RN r2
Len      RN r3
mulOp1   RN r4
mulOp2   RN r5
prodLo   RN r4
prodHi   RN r5
accuLo   RN r6
accuHi   RN r7
accuEx   RN r8
bPtr     RN r9
outStop  RN r10
aStart   RN r11
bStart   RN r12
    
loWord   RN r11
hiWord   RN r12
    
    
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; generic modular squaring ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    
gfp_sqr_asm FUNCTION
    
    push  {r0, r4-r12}       ; push callee-saved registers
    
    sub   sp, sp, Len, LSL #3
    mov   rPtr, sp
    
    ;; Initialize pointers and variables
    mov   aStart, aPtr       ;
    mov   bStart, aPtr       ;
    add   outStop, rPtr, Len, LSL #2 ;
    mov   accuLo, #0
    mov   accuHi, #0
    mov   accuEx, #0
    str   accuLo, [rPtr], #4
    
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
    cmp   aPtr, bPtr        ; check whether aPtr < bPtr
    blo   InnLoop1          ; if yes then branch back to InnLoop1
    
    str   accuLo, [rPtr], #4
    mov   accuLo, accuHi
    mov   accuHi, accuEx
    mov   accuEx, #0
    
    cmp   rPtr, outStop
    blt   OutLoop1
    
    ;; END OF FIRST OUTER LOOP
    
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
    cmp   aPtr, bPtr        ; check whether aPtr < bPtr
    blt   InnLoop2          ; if yes then branch back to InnLoop2
    
    str   accuLo, [rPtr], #4
    mov   accuLo, accuHi
    mov   accuHi, accuEx
    mov   accuEx, #0
    
    cmp   rPtr, outStop
    blt   OutLoop2
    
    ;; END OF SECOND OUTER LOOP
    
    stmia rPtr!, {accuLo, accuHi}
    
    ;; Double the result obtained so far and add the squares a[i]^2 from the
    ;; main diagonal
    
    add   outStop, bStart, #4   ; outStop contains now address of a[len]
    sub   aPtr, outStop, Len, LSL #2
    sub   rPtr, rPtr, Len, LSL #3
    
    ;; START OF THIRD OUTER LOOP
    
OutLoop3
    
    ldr   mulOp1, [aPtr], #4
    umull prodLo, prodHi, mulOp1, mulOp1
    ldmia rPtr, {accuLo, accuHi}
    adds  prodLo, prodLo, accuEx
    adc   prodHi, prodHi, #0
    mov   accuEx, #0
    adds  accuLo, accuLo, accuLo
    adcs  accuHi, accuHi, accuHi
    adc   accuEx, accuEx, #0
    adds  accuLo, accuLo, prodLo
    adcs  accuHi, accuHi, prodHi
    adc   accuEx, accuEx, #0
    stmia rPtr!, {accuLo, accuHi}
    cmp   aPtr, outStop
    blt   OutLoop3
    
    ;; END OF THIRD OUTER LOOP
    
    ;; Start of the reduction operation
    
    sub   aPtr, rPtr, Len, LSL #3
    sub   bPtr, rPtr, Len, LSL #2
    sub   outStop, rPtr, #4
    ldr   rPtr, [sp, +Len, LSL #3]
    lsl   cWord, cWord, #1
    mov   hiWord, #0
    
    ;; START OF FOURTH OUTER LOOP
    
OutLoop4
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
    blt   OutLoop4
    
    ;; END OF FOURTH OUTER LOOP
    
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
    
    ;; START OF FIFTH OUTER LOOP
    
OutLoop5
    ldr   loWord, [rPtr]
    adcs  loWord, loWord, #0
    str   loWord, [rPtr], #4
    teq   rPtr, outStop
    bne   OutLoop5
    
    ;; END OF FIFTH OUTER LOOP
    
    add   sp, sp, Len, LSL #3
    pop   {r0, r4-r12}
    bx    lr
    
    ENDFUNC
    
    
    END
