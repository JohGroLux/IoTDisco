///////////////////////////////////////////////////////////////////////////////
// moncurve.h: Function prototypes for arithmetic on twisted Edwards curves. //
// This file is part of SECC430, a Scalable ECC implementation for MSP430.   //
// Version 1.0.1 (2023-06-24), see <http://www.cryptolux.org/> for updates.  //
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


#ifndef _TEDCURVE_H
#define _TEDCURVE_H

#include "typedefs.h"

/***********************/
/* function prototypes */
/***********************/

void ted_copy(PROPOINT *r, const PROPOINT *p, int len, int num);
void ted_add(PROPOINT *p, const PROPOINT *q, const ECDPARAM *m);
void ted_double(PROPOINT *p, const ECDPARAM *m);
int  ted_validate(const PROPOINT *p, const ECDPARAM *m);
void ted_mul_binary(PROPOINT *r, const UINT16 *k, const AFFPOINT *p, const ECDPARAM *m);
void ted_mul_comb4b(PROPOINT *r, const UINT16 *k, const ECDPARAM *m);
void ted_to_mon(PROPOINT *r, const PROPOINT *p, const ECDPARAM *m);
int  ted_proj_affine(PROPOINT *r, const PROPOINT *p, const ECDPARAM *m);
int  ted_mul_varbase(AFFPOINT *q, const UINT16 *k, const AFFPOINT *p, const ECDPARAM *m);
int  ted_mul_fixbase(AFFPOINT *r, const UINT16 *k, const ECDPARAM *m);
int  ted_mul_dblbase(AFFPOINT *r, const DBLSCALAR *k, const AFFPOINT *p, const ECDPARAM *m);

#endif
