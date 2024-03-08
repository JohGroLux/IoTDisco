///////////////////////////////////////////////////////////////////////////////
// moncurve.h: Function prototypes for arithmetic on Montgomery curves.      //
// This file is part of SECC430, a Scalable ECC implementation for MSP430.   //
// Version 1.0.1 (2016-06-24), see <http://www.cryptolux.org/> for updates.  //
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


#ifndef _MONCURVE_H
#define _MONCURVE_H

#include "typedefs.h"

/***********************/
/* function prototypes */
/***********************/

void mon_add(PROPOINT *p, const PROPOINT *q, const UINT16 *x, const ECDPARAM *m);
void mon_double(PROPOINT *p, const ECDPARAM *m);
int  mon_check_order(PROPOINT *r, const UINT16 *xp, const ECDPARAM *m);
void mon_mul_ladder(PROPOINT *r, const UINT16 *k, const UINT16 *x, const ECDPARAM *m);
void mon_to_ted(PROPOINT *r, const PROPOINT *p, const ECDPARAM *m);
int  mon_proj_affine(PROPOINT *r, const PROPOINT *p, const ECDPARAM *m);
void mon_recover_y(PROPOINT *r, const PROPOINT *q, const PROPOINT *p, const ECDPARAM *m);
int  mon_mul_varbase(UINT16 *r, const UINT16 *k, const UINT16 *p, const ECDPARAM *m);
int mon_mul_fixbase(UINT16 *r, const UINT16 *k, const ECDPARAM *m);

void mon_test25519(void);
#endif
