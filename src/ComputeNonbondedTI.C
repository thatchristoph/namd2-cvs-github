/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
   Common operations for ComputeNonbonded classes
*/

#include "ComputeNonbondedInl.h"

// 3 inline functions to handle explicit calculation of separation-shifted
// vdW for FEP and TI, and a shifted electrostatics potential for decoupling

/* ********************************** */
/* vdW energy, force and dU/dl for TI */
/* ********************************** */
inline void ti_vdw_force_energy_dUdl (BigReal A, BigReal B, BigReal r2, 
  BigReal myVdwShift, BigReal switchdist2, BigReal cutoff2, 
  BigReal myVdwLambda, BigReal alchVdwShiftCoeff, BigReal switchfactor, 
  BigReal* alch_vdw_energy, BigReal* alch_vdw_force, BigReal* alch_vdw_dUdl) {
  const BigReal r2_1 = r2 + myVdwShift;  //myVdwShift already multplied by relevant (1-vdwLambda)
  const BigReal r6_1 = r2_1*r2_1*r2_1;
    
  // switching function (this is correct whether switching is active or not)
  const BigReal switchmul = r2 > switchdist2? \
           switchfactor*(cutoff2 - r2)*(cutoff2 - r2)*(cutoff2 - 3.*switchdist2 + 2.*r2) \
           : 1.;
  const BigReal switchmul2 = (r2 > switchdist2)? \
           12.*switchfactor*(cutoff2 - r2)*(r2 - switchdist2) : 0.;
  // separation-shifted vdW force and energy
  *alch_vdw_energy = A/(r6_1*r6_1) - B/r6_1;
  *alch_vdw_force =  myVdwLambda * (  \
                    + (12.*(*alch_vdw_energy) + 6.*B/r6_1)/r2_1 * switchmul \
                    + (*alch_vdw_energy) * switchmul2);
  *alch_vdw_dUdl = (myVdwLambda * alchVdwShiftCoeff * (
            (6*A/(r6_1*r6_1) - 3*B/r6_1) / r2_1 ) + (*alch_vdw_energy))
            * switchmul;  // dU/dlambda
  *alch_vdw_energy *= myVdwLambda*switchmul;
}


/*************THERMODYNAMIC INTEGRATION*************/
#define TIFLAG
#define CALCENERGY

#define NBTYPE NBPAIR
#include "ComputeNonbondedBase.h"
#define FULLELECT
#include "ComputeNonbondedBase.h"
#define MERGEELECT
#include "ComputeNonbondedBase.h"
#undef MERGEELECT
#define SLOWONLY
#include "ComputeNonbondedBase.h"
#undef SLOWONLY
#undef FULLELECT
#undef  NBTYPE

#define NBTYPE NBSELF
#include "ComputeNonbondedBase.h"
#define FULLELECT
#include "ComputeNonbondedBase.h"
#define MERGEELECT
#include "ComputeNonbondedBase.h"
#undef MERGEELECT
#define SLOWONLY
#include "ComputeNonbondedBase.h"
#undef SLOWONLY
#undef FULLELECT
#undef  NBTYPE

#undef CALCENERGY

#define NBTYPE NBPAIR
#include "ComputeNonbondedBase.h"
#define FULLELECT
#include "ComputeNonbondedBase.h"
#define MERGEELECT
#include "ComputeNonbondedBase.h"
#undef MERGEELECT
#define SLOWONLY
#include "ComputeNonbondedBase.h"
#undef SLOWONLY
#undef FULLELECT
#undef  NBTYPE

#define NBTYPE NBSELF
#include "ComputeNonbondedBase.h"
#define FULLELECT
#include "ComputeNonbondedBase.h"
#define MERGEELECT
#include "ComputeNonbondedBase.h"
#undef MERGEELECT
#define SLOWONLY
#include "ComputeNonbondedBase.h"
#undef SLOWONLY
#undef FULLELECT
#undef  NBTYPE

#undef TIFLAG

