/*


!
!  Dalton, a molecular electronic structure program
!  Copyright (C) 2018 by the authors of Dalton.
!
!  This program is free software; you can redistribute it and/or
!  modify it under the terms of the GNU Lesser General Public
!  License version 2.1 as published by the Free Software Foundation.
!
!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!  Lesser General Public License for more details.
!
!  If a copy of the GNU LGPL v2.1 was not distributed with this
!  code, you can obtain one at https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html.
!

!

*/
/* fun-test.c:
   implementation of a test GGA-class functional.
   This is a second example functional.
   (c) Pawel Salek, pawsa@theochem.kth.se, aug 2001
   NOTE:
   this file may seem unnecessarily complex but the structure really pays off
   when implementing multiple functionals depending on different parameters.
*/

#include <stddef.h>
#include "general.h"
#define __CVERSION__

#include "functionals.h"

/* INTERFACE PART */
static integer example9_isgga(void) { return 1; }
static integer example9_read(const char* conf_line);
static real example9_energy(const FunDensProp* dp);
static void example9_first(FunFirstFuncDrv *ds,   real factor, const FunDensProp* dp);
static void example9_second(FunSecondFuncDrv *ds, real factor, const FunDensProp* dp);
static void example9_third(FunThirdFuncDrv *ds,   real factor, const FunDensProp* dp);

Functional Example9Functional = {
  "Example9",         /* name */
  example9_isgga,     /* gga-corrected */
   1,
  example9_read, 
  NULL,              /* reporter */
  example9_energy, 
  example9_first,
  example9_second,
  example9_third
};

/* IMPLEMENTATION PART */
static integer
example9_read(const char* conf_line)
{
  fun_set_hf_weight(0.0);
  return 1;
}

static const real EPREF= -5e-2;

static real
example9_energy(const FunDensProp* dp)
{
  return EPREF*(dp->gradab);
}
/* example_first:
   derivatives with respect to rho_alpha, and grad_alpha
 */
static void
example9_first(FunFirstFuncDrv *ds, real factor, const FunDensProp* dp)
{
  ds->df00001 +=  EPREF*factor;
}

static void
example9_second(FunSecondFuncDrv *ds, real factor, const FunDensProp* dp)
{
  ds->df00001 +=  EPREF*factor;
}

/* example_third:
   Test functional derivatives.
*/
static void
example9_third(FunThirdFuncDrv *ds, real factor, const FunDensProp* dp)
{
  ds->df00001 +=  EPREF*factor;
}
