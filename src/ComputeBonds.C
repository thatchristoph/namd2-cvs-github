/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include "InfoStream.h"
#include "ComputeBonds.h"
#include "Molecule.h"
#include "Parameters.h"
#include "Node.h"
#include "ReductionMgr.h"
#include "Lattice.h"
#include "PressureProfile.h"
#include "Debug.h"


// static initialization
int BondElem::pressureProfileSlabs = 0;
int BondElem::pressureProfileAtomTypes = 1;
BigReal BondElem::pressureProfileThickness = 0;
BigReal BondElem::pressureProfileMin = 0;

void BondElem::getMoleculePointers
    (Molecule* mol, int* count, int32*** byatom, Bond** structarray)
{
#ifdef MEM_OPT_VERSION
  NAMD_die("Should not be called in BondElem::getMoleculePointers in memory optimized version!");
#else
  *count = mol->numBonds;
  *byatom = mol->bondsByAtom;
  *structarray = mol->bonds;
#endif
}

void BondElem::getParameterPointers(Parameters *p, const BondValue **v) {
  *v = p->bond_array;
}

void BondElem::computeForce(BigReal *reduction, 
                            BigReal *pressureProfileData)
{
  DebugM(1, "::computeForce() localIndex = " << localIndex[0] << " "
               << localIndex[1] << std::endl);

  // skip Lonepair bonds (other k=0. bonds have been filtered out)
  if (0. == value->k) return;

  BigReal scal;         // force scaling
  BigReal energy;	// energy from the bond

  //DebugM(3, "::computeForce() -- starting with bond type " << bondType << std::endl);

  // get the bond information
  Real k = value->k * scale;
  Real x0 = value->x0;

  // compute vectors between atoms and their distances
  const Lattice & lattice = p[0]->p->lattice;
  const Vector r12 = lattice.delta(p[0]->x[localIndex[0]].position,
					p[1]->x[localIndex[1]].position);

  if (0. == x0) {  // for Drude bonds
    SimParameters *simParams = Node::Object()->simParameters;
    BigReal drudeBondLen = simParams->drudeBondLen;
    BigReal drudeBondConst = simParams->drudeBondConst;

    BigReal r2 = r12.length2();

    scal = -2.0*k;    // bond interaction for equilibrium length 0
    energy = k * r2;

    if ( (drudeBondConst > 0) && ( ! simParams->drudeHardWallOn ) &&
        (r2 > drudeBondLen*drudeBondLen) ) {
      // add a quartic restraining potential to keep Drude bond short
      BigReal r = sqrt(r2);
      BigReal diff = r - drudeBondLen;
      BigReal diff2 = diff*diff;

      scal += -4*drudeBondConst * diff2 * diff / r;
      energy += drudeBondConst * diff2 * diff2;
    }
  }
  else {
    BigReal r = r12.length();  // Distance between atoms
    BigReal diff = r - x0;     // Compare it to the rest bond

    //  Add the energy from this bond to the total energy
    energy = k*diff*diff;

    //  Determine the magnitude of the force
    diff *= -2.0*k;

    //  Scale the force vector accordingly
    scal = (diff/r);
  }
  const Force f12 = scal * r12;


  //  Now add the forces to each force vector
  p[0]->f[localIndex[0]] += f12;
  p[1]->f[localIndex[1]] -= f12;

  DebugM(3, "::computeForce() -- ending with delta energy " << energy << std::endl);
  reduction[bondEnergyIndex] += energy;
  reduction[virialIndex_XX] += f12.x * r12.x;
  reduction[virialIndex_XY] += f12.x * r12.y;
  reduction[virialIndex_XZ] += f12.x * r12.z;
  reduction[virialIndex_YX] += f12.y * r12.x;
  reduction[virialIndex_YY] += f12.y * r12.y;
  reduction[virialIndex_YZ] += f12.y * r12.z;
  reduction[virialIndex_ZX] += f12.z * r12.x;
  reduction[virialIndex_ZY] += f12.z * r12.y;
  reduction[virialIndex_ZZ] += f12.z * r12.z;

  if (pressureProfileData) {
    BigReal z1 = p[0]->x[localIndex[0]].position.z;
    BigReal z2 = p[1]->x[localIndex[1]].position.z;
    int n1 = (int)floor((z1-pressureProfileMin)/pressureProfileThickness);
    int n2 = (int)floor((z2-pressureProfileMin)/pressureProfileThickness);
    pp_clamp(n1, pressureProfileSlabs);
    pp_clamp(n2, pressureProfileSlabs);
    int p1 = p[0]->x[localIndex[0]].partition;
    int p2 = p[1]->x[localIndex[1]].partition;
    int pn = pressureProfileAtomTypes;
    pp_reduction(pressureProfileSlabs,
                n1, n2, 
                p1, p2, pn, 
                f12.x * r12.x, f12.y * r12.y, f12.z * r12.z,
                pressureProfileData);
  } 
}

void BondElem::submitReductionData(BigReal *data, SubmitReduction *reduction)
{
  reduction->item(REDUCTION_BOND_ENERGY) += data[bondEnergyIndex];
  ADD_TENSOR(reduction,REDUCTION_VIRIAL_NORMAL,data,virialIndex);
}

