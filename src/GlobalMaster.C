/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include "Node.h"
#include "Molecule.h"
#include "NamdTypes.h"
#include "GlobalMaster.h"

//#define DEBUGM
#define MIN_DEBUG_LEVEL 1
#include "Debug.h"

void GlobalMaster::processData(AtomIDList::iterator a_i,
			       AtomIDList::iterator a_e,
			       PositionList::iterator p_i,
			       PositionList::iterator g_i,
			       PositionList::iterator g_e,
			       AtomIDList::iterator last_atoms_forced_i,
			       AtomIDList::iterator last_atoms_forced_e,
			       ForceList::iterator last_forces_i,
			       AtomIDList::iterator forceid_i,
			       AtomIDList::iterator forceid_e,
			       ForceList::iterator totalforce_i) {
  atomIdBegin = a_i;
  atomIdEnd = a_e;
  atomPositionBegin = p_i;
  groupPositionBegin = g_i;
  groupPositionEnd = g_e;
  lastAtomsForcedBegin = last_atoms_forced_i;
  lastAtomsForcedEnd = last_atoms_forced_e;
  lastForcesBegin = last_forces_i;
  forceIdBegin = forceid_i;
  forceIdEnd = forceid_e;
  totalForceBegin = totalforce_i;

  calculate();

  /* check to make sure the force arrays still match */
  if(appForcesChanged) {
    if(fAtoms.size() != appForces.size())
      NAMD_die("# of atoms forced != # of forces given");
  }
  if(appForcesChanged) {
    if(reqGroups.size() != grpForces.size())
      NAMD_die("# of groups forced != # of groups requested");
  }

  /* check if the groups have changed, and update the masses */
  if(reqGroupsChanged) {
    groupMasses.resize(0);

    // is there a way to do this non-globally?
    Molecule *mol = Node::Object()->molecule;

    // update each group mass
    int g;
    for(g=0;g<reqGroups.size();g++) {
      AtomIDList &atom_list = reqGroups[g];
      BigReal mass = 0; // the total mass of the group
      int a;
      for(a=0;a<atom_list.size();a++) { // get the total
	mass += mol->atommass(a);
      }
      groupMasses.add(mass); // add the mass to the group
    }
  }
}

void GlobalMaster::check() const {
  /* check to make sure the force arrays still match */
  if(fAtoms.size() != appForces.size())
    NAMD_die("# of atoms forced != # of forces given");
  if(reqGroups.size() != grpForces.size())
    NAMD_die("# of groups forced != # of groups requested");
}

void GlobalMaster::clearChanged() {
  reqAtomsChanged = false;
  appForcesChanged = false;
  reqGroupsChanged = false;
}

void GlobalMaster::calculate() {
  NAMD_die("Internal error: pure virtual function called");
}

GlobalMaster::GlobalMaster() {
  step = -1;
  clearChanged();
}

bool GlobalMaster::changedAtoms() {
  return reqAtomsChanged;
}

bool GlobalMaster::changedForces() {
  return appForcesChanged;
}

bool GlobalMaster::changedGroups() {
  return reqGroupsChanged;
}

const AtomIDList &GlobalMaster::requestedAtoms() {
  return reqAtoms;
}

AtomIDList &GlobalMaster::modifyRequestedAtoms() {
  reqAtomsChanged = true;
  return reqAtoms;
}

const AtomIDList &GlobalMaster::forcedAtoms() {
  return fAtoms;
}

const ForceList &GlobalMaster::appliedForces() {
  return appForces;
}

const ForceList &GlobalMaster::groupForces() {
  return grpForces;
}

const ResizeArray<AtomIDList> &GlobalMaster::requestedGroups() {
  return reqGroups;
}

AtomIDList &GlobalMaster::modifyForcedAtoms() {
  appForcesChanged = true;
  return fAtoms;
}

ForceList &GlobalMaster::modifyAppliedForces() {
  appForcesChanged = true;
  return appForces;
}

ForceList &GlobalMaster::modifyGroupForces() {
  // XXX should we mark something else here?
  appForcesChanged = true;
  return grpForces;
}

ResizeArray<AtomIDList> &GlobalMaster::modifyRequestedGroups() {
  reqGroupsChanged = true;
  DebugM(1,"Groups have changed.\n");
  return reqGroups;
}

AtomIDList::const_iterator GlobalMaster::getAtomIdBegin() {
  return atomIdBegin;
}

AtomIDList::const_iterator GlobalMaster::getAtomIdEnd() {
  return atomIdEnd;
}

PositionList::const_iterator GlobalMaster::getAtomPositionBegin() {
  return atomPositionBegin;
}

PositionList::const_iterator GlobalMaster::getGroupPositionBegin() {
  return groupPositionBegin;
}

PositionList::const_iterator GlobalMaster::getGroupPositionEnd() {
  return groupPositionEnd;
}

ResizeArray<BigReal>::const_iterator GlobalMaster::getGroupMassBegin()
{
  return groupMasses.begin();
}

ResizeArray<BigReal>::const_iterator GlobalMaster::getGroupMassEnd() {
  return groupMasses.end();
}

AtomIDList::const_iterator GlobalMaster::getLastAtomsForcedBegin() {
  return lastAtomsForcedBegin;
}

AtomIDList::const_iterator GlobalMaster::getLastAtomsForcedEnd() {
  return lastAtomsForcedEnd;
}

ForceList::const_iterator GlobalMaster::getLastForcesBegin() {
  return lastForcesBegin;
}

AtomIDList::const_iterator GlobalMaster::getForceIdBegin()
{
  return forceIdBegin;
}

AtomIDList::const_iterator GlobalMaster::getForceIdEnd()
{
  return forceIdEnd;
}

ForceList::const_iterator GlobalMaster::getTotalForce()
{
  return totalForceBegin;
}
