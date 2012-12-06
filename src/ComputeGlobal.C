/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
   Forwards atoms to master node for force evaluation.
*/

#include "InfoStream.h"
#include "Node.h"
#include "PatchMap.h"
#include "PatchMap.inl"
#include "AtomMap.h"
#include "ComputeGlobal.h"
#include "ComputeGlobalMsgs.h"
#include "PatchMgr.h"
#include "Molecule.h"
#include "ReductionMgr.h"
#include "ComputeMgr.h"
#include "ComputeMgr.decl.h"
#include "SimParameters.h"
#include <stdio.h>

//#define DEBUGM
#define MIN_DEBUG_LEVEL 1
#include "Debug.h"

// CLIENTS

ComputeGlobal::ComputeGlobal(ComputeID c, ComputeMgr *m)
	: ComputeHomePatches(c)
{
  DebugM(3,"Constructing client\n");
  aid.resize(0);
  gdef.resize(0);
  comm = m;
  firsttime = 1;
  int i, numAtoms=Node::Object()->molecule->numAtoms;
  isRequested = new int[numAtoms];
  for (i=0; i<numAtoms; ++i)
    isRequested[i] = 0;
  SimParameters *sp = Node::Object()->simParameters;
  dofull = (sp->GBISserOn || sp->GBISOn || sp->fullDirectOn || sp->FMAOn || sp->PMEOn);
  fid = new AtomIDList;
  totalForce = new ForceList;
  reduction = ReductionMgr::Object()->willSubmit(REDUCTIONS_BASIC);
}

ComputeGlobal::~ComputeGlobal()
{
  delete[] isRequested;
  delete fid;
  delete totalForce;
  delete reduction;
}

void ComputeGlobal::configure(AtomIDList newaid, AtomIDList newgdef) {
  DebugM(4,"Receiving configuration (" << newaid.size() <<
	" atoms and " << newgdef.size() << " atoms/groups) on client\n");

  AtomIDList::iterator a, a_e;
  
  for (a=aid.begin(),a_e=aid.end(); a!=a_e; ++a)
    isRequested[*a] = 0;

  // store data
  aid = newaid;
  gdef = newgdef;
  
  for (a=aid.begin(),a_e=aid.end(); a!=a_e; ++a)
    isRequested[*a] = 1;

  // calculate group masses
  Molecule *mol = Node::Object()->molecule;
  gmass.resize(0);
  AtomIDList::iterator g_i, g_e;
  g_i = gdef.begin(); g_e = gdef.end();
  for ( ; g_i != g_e; ++g_i ) {
    BigReal mass = 0;
    for ( ; *g_i != -1; ++g_i ) {
      mass += mol->atommass(*g_i);
    }
    gmass.add(mass);
  }
}

#if 0
void ComputeGlobal::recvConfig(ComputeGlobalConfigMsg *msg) {
  DebugM(3,"Receiving configure on client\n");
  configure(msg->aid,msg->gdef);
  delete msg;
  sendData();
}
#endif

void ComputeGlobal::recvResults(ComputeGlobalResultsMsg *msg) {
  DebugM(3,"Receiving results (" << msg->aid.size() << " forces, "
	 << msg->newgdef.size() << " new group atoms) on client\n");

  // set the forces only if we aren't going to resend the data
  int setForces = !msg->resendCoordinates;

  if(setForces) { // we are requested to 
    // Store forces to patches
    PatchMap *patchMap = PatchMap::Object();
    int numPatches = patchMap->numPatches();
    AtomMap *atomMap = AtomMap::Object();
    const Lattice & lattice = patchList[0].p->lattice;
    ResizeArrayIter<PatchElem> ap(patchList);
    Force **f = new Force*[numPatches];
    FullAtom **t = new FullAtom*[numPatches];
    for ( int i = 0; i < numPatches; ++i ) { f[i] = 0; t[i] = 0; }
    Force extForce = 0.;
    Tensor extVirial;

    for (ap = ap.begin(); ap != ap.end(); ap++) {
      (*ap).r = (*ap).forceBox->open();
      f[(*ap).patchID] = (*ap).r->f[Results::normal];
      t[(*ap).patchID] = (*ap).p->getAtomList().begin();
    }

    AtomIDList::iterator a = msg->aid.begin();
    AtomIDList::iterator a_e = msg->aid.end();
    ForceList::iterator f2 = msg->f.begin();
    for ( ; a != a_e; ++a, ++f2 ) {
      DebugM(1,"processing atom "<<(*a)<<", F="<<(*f2)<<"...\n");
      /* XXX if (*a) is out of bounds here we get a segfault */
      LocalID localID = atomMap->localID(*a);
      if ( localID.pid == notUsed || ! f[localID.pid] ) continue;
      Force f_atom = (*f2);
      f[localID.pid][localID.index] += f_atom;
      Position x_orig = t[localID.pid][localID.index].position;
      Transform trans = t[localID.pid][localID.index].transform;
      Position x_atom = lattice.reverse_transform(x_orig,trans);
      extForce += f_atom;
      extVirial += outer(f_atom,x_atom);
    }
    DebugM(1,"done with the loop\n");

  // calculate forces for atoms in groups
    Molecule *mol = Node::Object()->molecule;
    AtomIDList::iterator g_i, g_e;
    g_i = gdef.begin(); g_e = gdef.end();
    ResizeArray<BigReal>::iterator gm_i = gmass.begin();
    ForceList::iterator gf_i = msg->gforce.begin();
    //iout << iDEBUG << "recvResults\n" << endi;
    for ( ; g_i != g_e; ++g_i, ++gm_i, ++gf_i ) {
      //iout << iDEBUG << *gf_i << '\n' << endi;
      Vector accel = (*gf_i) / (*gm_i);
      for ( ; *g_i != -1; ++g_i ) {
	//iout << iDEBUG << *g_i << '\n' << endi;
	LocalID localID = atomMap->localID(*g_i);
	if ( localID.pid == notUsed || ! f[localID.pid] ) continue;
	Force f_atom = accel * mol->atommass(*g_i);
	f[localID.pid][localID.index] += f_atom;
        Position x_orig = t[localID.pid][localID.index].position;
        Transform trans = t[localID.pid][localID.index].transform;
        Position x_atom = lattice.reverse_transform(x_orig,trans);
        extForce += f_atom;
        extVirial += outer(f_atom,x_atom);
      }
    }
    DebugM(1,"done with the groups\n");

    for (ap = ap.begin(); ap != ap.end(); ap++) {
      (*ap).forceBox->close(&((*ap).r));
    }

    delete [] f;
    delete [] t;

    ADD_VECTOR_OBJECT(reduction,REDUCTION_EXT_FORCE_NORMAL,extForce);
    ADD_TENSOR_OBJECT(reduction,REDUCTION_VIRIAL_NORMAL,extVirial);
    reduction->submit();
  }
  // done setting the forces

  // Get reconfiguration if present
  if ( msg->reconfig ) configure(msg->newaid, msg->newgdef);

  // send another round of data if requested

  if(msg->resendCoordinates) {
    DebugM(3,"Sending requested data right away\n");
    sendData();
  }

  delete msg;
  DebugM(3,"Done processing results\n");
}

void ComputeGlobal::doWork()
{
  DebugM(2,"doWork\n");
  if(!firsttime) sendData();
  else {
    ComputeGlobalDataMsg *msg = new ComputeGlobalDataMsg;
    comm->sendComputeGlobalData(msg);
    firsttime = 0;
  }
  DebugM(2,"done with doWork\n");
}

void ComputeGlobal::sendData()
{
  DebugM(2,"sendData\n");
  // Get positions from patches
  PatchMap *patchMap = PatchMap::Object();
  int numPatches = patchMap->numPatches();
  AtomMap *atomMap = AtomMap::Object();
  const Lattice & lattice = patchList[0].p->lattice;
  ResizeArrayIter<PatchElem> ap(patchList);
  CompAtom **x = new CompAtom*[numPatches];
  FullAtom **t = new FullAtom*[numPatches];
  for ( int i = 0; i < numPatches; ++i ) { x[i] = 0; t[i] = 0; }

  int step = -1;
  for (ap = ap.begin(); ap != ap.end(); ap++) {
    x[(*ap).patchID] = (*ap).positionBox->open();
    t[(*ap).patchID] = (*ap).p->getAtomList().begin();
    step = (*ap).p->flags.step;
  }

  ComputeGlobalDataMsg *msg = new  ComputeGlobalDataMsg;

  msg->step = step;

  AtomIDList::iterator a = aid.begin();
  AtomIDList::iterator a_e = aid.end();
  for ( ; a != a_e; ++a ) {
    LocalID localID = atomMap->localID(*a);
    if ( localID.pid == notUsed || ! x[localID.pid] ) continue;
    msg->aid.add(*a);
    Position x_orig = x[localID.pid][localID.index].position;
    Transform trans = t[localID.pid][localID.index].transform;
    msg->p.add(lattice.reverse_transform(x_orig,trans));
  }

  // calculate group centers of mass
  Molecule *mol = Node::Object()->molecule;
  AtomIDList::iterator g_i, g_e;
  g_i = gdef.begin(); g_e = gdef.end();
  ResizeArray<BigReal>::iterator gm_i = gmass.begin();
  for ( ; g_i != g_e; ++g_i, ++gm_i ) {
    Vector com(0,0,0);
    for ( ; *g_i != -1; ++g_i ) {
      LocalID localID = atomMap->localID(*g_i);
      if ( localID.pid == notUsed || ! x[localID.pid] ) continue;
      Position x_orig = x[localID.pid][localID.index].position;
      Transform trans = t[localID.pid][localID.index].transform;
      com += lattice.reverse_transform(x_orig,trans) * mol->atommass(*g_i);
    }
    com /= *gm_i;
    DebugM(1,"Adding center of mass "<<com<<"\n");
    msg->gcom.add(com);
  }

  for (ap = ap.begin(); ap != ap.end(); ap++) {
    (*ap).positionBox->close(&(x[(*ap).patchID]));
  }
  
  msg->fid = *fid;
  msg->tf = *totalForce;
  delete fid;
  delete totalForce;
  fid = new AtomIDList;
  totalForce = new ForceList;

  delete [] x;
  delete [] t;

  DebugM(3,"Sending data (" << msg->aid.size() << " positions) on client\n");
  comm->sendComputeGlobalData(msg);
}


// This function is called by each HomePatch after force
// evaluation. It stores the indices and forces of the requested
// atoms here, to be sent to GlobalMasterServer during the next
// time step. The total force is the sum of three components:
// "normal", "nbond" and "slow", the latter two may be calculated
// less frequently, so their most recent values are stored in
// "f_saved" and used here. If we don't do full electrostatics,
// there's no "slow" part.
void ComputeGlobal::saveTotalForces(HomePatch *homePatch)
{
  if ( Node::Object()->simParameters->accelMDOn && Node::Object()->simParameters->accelMDDebugOn && Node::Object()->simParameters->accelMDdihe ) {
    int i, index, num=homePatch->numAtoms;
    FullAtomList atoms = homePatch->atom;
    ForceList af=homePatch->f[Results::amdf];

    for (i=0; i<num; ++i)
      if (isRequested[index=atoms[i].id]) {
        fid->add(index);
        totalForce->add(af[i]);
      }
    return;
  }

  int fixedAtomsOn = Node::Object()->simParameters->fixedAtomsOn;
  int i, index, num=homePatch->numAtoms;
  FullAtomList atoms = homePatch->atom;
  ForceList f1=homePatch->f[Results::normal], f2=homePatch->f_saved[Results::nbond],
            f3=homePatch->f_saved[Results::slow];
  Force f_sum;
  
  for (i=0; i<num; ++i)
    if (isRequested[index=atoms[i].id])
    { f_sum = f1[i]+f2[i];
      if (dofull)
        f_sum += f3[i];
      if ( fixedAtomsOn && atoms[i].atomFixed ) f_sum = 0.;
      fid->add(index);
      totalForce->add(f_sum);
    }
}
