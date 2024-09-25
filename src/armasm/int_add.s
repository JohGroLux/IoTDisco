///////////////////////////////////////////////////////////////////////////////
// int_add.s: Addition of two Multiple-Precision Integers (returns Carry).   //
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


    AREA intarith, CODE, READONLY ; Name this block of code ARMex
    
    EXPORT int_add_asm
    ALIGN 2
    
    
int_add_asm FUNCTION
    push {r4, r5}
    add  r3, r0, r3, LSL #2 ; r3 contains now the address of r[len]
    adds r4, r4, #0         ; clear C flag
loop
    ldr  r4, [r1], #4
    ldr  r5, [r2], #4
    adcs r4, r4, r5
    str  r4, [r0], #4
    teq  r0, r3             ; check whether r0 exceeds address of r[len-1]
    bne  loop               ; if not then branch back to beginning of loop
    mov  r0, #0
    adc  r0, r0, r0
    pop  {r4, r5}
    bx   lr
    ENDFUNC
    
    
    END
