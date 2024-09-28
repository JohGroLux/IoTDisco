;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; int_sqr.s: Multiple-Precision Squaring (Product-Scanning Method).         ;;
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


    AREA intarith, CODE, READONLY ; Name this block of code ARMex
    
    EXPORT int_sqr_asm
    ALIGN 2
    
    
;;;;;;;;;;;;;;;;;;;;
;; Register Names ;;
;;;;;;;;;;;;;;;;;;;;
    
rPtr     RN r0
aPtr     RN r1
Len      RN r2
bPtr     RN r3
mulOp1   RN r4
mulOp2   RN r5
prodLo   RN r4
prodHi   RN r5
accuLo   RN r6
accuHi   RN r7
accuEx   RN r8
aStart   RN r9
bStart   RN r10
outStop  RN r11
    
    
;;;;;;;;;;;;;;;;;;;;;;
;; Generic Squaring ;;
;;;;;;;;;;;;;;;;;;;;;;
    
int_sqr_asm PROC
    
    push  {r4-r12}          ; push callee-saved registers
    
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
    
    pop   {r4-r12}
    bx    lr
    
    ENDP
    
    
    END
