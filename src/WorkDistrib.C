/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*****************************************************************************
 * $Source: /Projects/namd2/cvsroot/namd2/src/WorkDistrib.C,v $
 * $Author: jlai7 $
 * $Date: 2012/11/27 21:13:20 $
 * $Revision: 1.1251 $
 *****************************************************************************/

/** \file WorkDistrib.C
 *  Currently, WorkDistrib generates the layout of the Patches,
 *  directs the construction and distribution of Computes and
 *  associates Computes with Patches.
 */

#include <stdio.h>

#include "InfoStream.h"
#include "ProcessorPrivate.h"
#include "BOCgroup.h"
#include "WorkDistrib.decl.h"
#include "WorkDistrib.h"
#include "Lattice.h"
#include "ComputeMsmMsa.h"  // needed for MsmMsaData definition
#include "main.decl.h"
#include "main.h"
#include "Node.h"
#include "PatchMgr.h"
#include "PatchMap.inl"
#include "NamdTypes.h"
#include "PDB.h"
#include "SimParameters.h"
#include "Molecule.h"
#include "NamdOneTools.h"
#include "Compute.h"
#include "ComputeMap.h"
#include "RecBisection.h"
#include "Random.h"
#include "varsizemsg.h"
#include "ProxyMgr.h"
#include "Priorities.h"
#include "SortAtoms.h"

//#define DEBUGM
#define MIN_DEBUG_LEVEL 2
#include "Debug.h"
#ifdef MEM_OPT_VERSION
extern int isOutputProcessor(int); 
#endif
class ComputeMapChangeMsg : public CMessage_ComputeMapChangeMsg
{
public:

  int numNewNodes;
  int numNewNumPartitions;
  int *newNodes;
  char *newNumPartitions;

//  VARSIZE_DECL(ComputeMapChangeMsg);
};

/*
VARSIZE_MSG(ComputeMapChangeMsg,
  VARSIZE_ARRAY(newNodes);
)
*/

static int eventMachineProgress;

//======================================================================
// Public functions
//----------------------------------------------------------------------
WorkDistrib::WorkDistrib()
{
  CkpvAccess(BOCclass_group).workDistrib = thisgroup;
  patchMapArrived = false;
  computeMapArrived = false;
  eventMachineProgress = traceRegisterUserEvent("CmiMachineProgressImpl",233);
}

//----------------------------------------------------------------------
WorkDistrib::~WorkDistrib(void)
{ }


//----------------------------------------------------------------------
void WorkDistrib::saveComputeMapChanges(int ep, CkGroupID chareID)
{
  saveComputeMapReturnEP = ep;
  saveComputeMapReturnChareID = chareID;

  ComputeMap *computeMap = ComputeMap::Object();

  int i;
  int nc = computeMap->numComputes();
  
  ComputeMapChangeMsg *mapMsg = new (nc, nc, 0) ComputeMapChangeMsg ;

  mapMsg->numNewNodes = nc;
  for(i=0; i<nc; i++)
    mapMsg->newNodes[i] = computeMap->newNode(i);
  mapMsg->numNewNumPartitions = nc;
  for(i=0; i<nc; i++)
    mapMsg->newNumPartitions[i] = computeMap->newNumPartitions(i);

  CProxy_WorkDistrib(thisgroup).recvComputeMapChanges(mapMsg);

/*
    // store the latest compute map
  SimParameters *simParams = Node::Object()->simParameters;
  if (simParams->storeComputeMap) {
    computeMap->saveComputeMap(simParams->computeMapFilename);
    CkPrintf("ComputeMap has been stored in %s.\n", simParams->computeMapFilename);
  }
*/
}

void WorkDistrib::recvComputeMapChanges(ComputeMapChangeMsg *msg) {
  
  if ( ! CkMyRank() ) {
    ComputeMap *computeMap = ComputeMap::Object();
    int nc = computeMap->numComputes();
    if ( nc != msg->numNewNodes ) NAMD_bug("recvComputeMapChanges 1");
    int i;
    for(i=0; i<nc; i++)
      computeMap->setNewNode(i,msg->newNodes[i]);
    if ( msg->numNewNumPartitions ) {
      if ( nc != msg->numNewNumPartitions ) NAMD_bug("recvComputeMapChanges 2");
      for(i=0; i<nc; i++)
        computeMap->setNewNumPartitions(i,msg->newNumPartitions[i]);
    }
  }

  delete msg;

  CkCallback cb(CkIndex_WorkDistrib::doneSaveComputeMap(NULL), 0, thisgroup);
  contribute(0, NULL, CkReduction::random, cb);
}

void WorkDistrib::doneSaveComputeMap(CkReductionMsg *msg) {
  delete msg;

  CkSendMsgBranch(saveComputeMapReturnEP, CkAllocMsg(0,0,0), 0, saveComputeMapReturnChareID);
}

#ifdef MEM_OPT_VERSION
//All basic info already exists for each atom inside the FullAtomList because
//it is loaded when reading the binary per-atom file. This function will fill
//the info regarding transform, nonbondedGroupSize etc. Refer to
//WorkDistrib::createAtomLists
void WorkDistrib::fillAtomListForOnePatch(int pid, FullAtomList &alist){
  PatchMap *patchMap = PatchMap::Object();

  ScaledPosition center(0.5*(patchMap->min_a(pid)+patchMap->max_a(pid)),
                          0.5*(patchMap->min_b(pid)+patchMap->max_b(pid)),
                          0.5*(patchMap->min_c(pid)+patchMap->max_c(pid)));

    int n = alist.size();
    FullAtom *a = alist.begin();
/* 
    //Those options are not supported in MEM_OPT_VERSIOn -Chao Mei 
//Modifications for alchemical fep
    Bool alchFepOn = params->alchFepOn;
    Bool alchThermIntOn = params->alchThermIntOn;
//fepe
    Bool lesOn = params->lesOn;
    Bool pairInteractionOn = params->pairInteractionOn;

    Bool pressureProfileTypes = (params->pressureProfileAtomTypes > 1);
*/
    SimParameters *params = Node::Object()->simParameters;
    const Lattice lattice = params->lattice;
    Transform mother_transform;
    for(int j=0; j < n; j++)
    {
      int aid = a[j].id;
      a[j].nonbondedGroupSize = 0;  // must be set based on coordinates
      
      // a[j].fixedPosition = a[j].position;  ParallelIOMgr stores ref coord here.

      if ( a[j].migrationGroupSize ) {
       if ( a[j].migrationGroupSize != a[j].hydrogenGroupSize ) {
            Position pos = a[j].position;
            int mgs = a[j].migrationGroupSize;
            int c = 1;

            for ( int k=a[j].hydrogenGroupSize; k<mgs;
                                k+=a[j+k].hydrogenGroupSize ) {
              pos += a[j+k].position;
              ++c;
            }

            pos *= 1./c;
            mother_transform = a[j].transform;  // should be 0,0,0
            pos = lattice.nearest(pos,center,&mother_transform);
            a[j].position = lattice.apply_transform(a[j].position,mother_transform);
            a[j].transform = mother_transform;        
       } else {
        a[j].position = lattice.nearest(a[j].position, center, &(a[j].transform));
        mother_transform = a[j].transform;
       }
      } else {
        a[j].position = lattice.apply_transform(a[j].position,mother_transform);
        a[j].transform = mother_transform;
      }

/* 
    //Those options are not supported in MEM_OPT_VERSIOn -Chao Mei 
//Modifications for alchemical fep
      if ( alchFepOn || alchThermIntOn || lesOn || pairInteractionOn || pressureProfileTypes) {
        a[j].partition = molecule->get_fep_type(aid);
      } 
      else {
        a[j].partition = 0;
      }
//fepe
*/
      a[j].partition = 0;

      //set langevinParams based on atom status
      if(params->langevinOn) {
        BigReal bval = params->langevinDamping;
        if(!params->langevinHydrogen && 
           ((a[j].status & HydrogenAtom)!=0)) {
          bval = 0;
        }else if ((a[j].status & LonepairAtom)!=0) {
          bval = 0;
        }else if ((a[j].status & DrudeAtom)!=0) {
          bval = params->langevinDamping;
        }
        a[j].langevinParam = bval;
      }

    }

    int size, allfixed;
    for(int j=0; j < n; j+=size) {
      size = a[j].hydrogenGroupSize;
      if ( ! size ) {
        NAMD_bug("Mother atom with hydrogenGroupSize of 0!");
      }
      allfixed = 1;
      for (int k = 0; k < size; ++k ) {
        allfixed = ( allfixed && (a[j+k].atomFixed) );
      }
      for (int k = 0; k < size; ++k ) {
        a[j+k].groupFixed = allfixed ? 1 : 0;
      }
    }

    if ( params->outputPatchDetails ) {
      int patchId = pid;
      int numAtomsInPatch = n;
      int numFixedAtomsInPatch = 0;
      int numAtomsInFixedGroupsInPatch = 0;
      for(int j=0; j < n; j++) {
        numFixedAtomsInPatch += ( a[j].atomFixed ? 1 : 0 );
        numAtomsInFixedGroupsInPatch += ( a[j].groupFixed ? 1 : 0 );
      }
      iout << "PATCH_DETAILS:"
           << " on proc " << CkMyPe()
           << " patch " << patchId
           << " atoms " << numAtomsInPatch
           << " fixed_atoms " << numFixedAtomsInPatch
           << " fixed_groups " << numAtomsInFixedGroupsInPatch
           << "\n" << endi;
    }

}

void WorkDistrib::random_velocities_parallel(BigReal Temp,InputAtomList &inAtoms)
{
  int i, j;             //  Loop counter
  BigReal kbT;          //  Boltzman constant * Temp
  BigReal randnum;      //  Random number from -6.0 to 6.0
  BigReal kbToverM;     //  sqrt(Kb*Temp/Mass)
  SimParameters *simParams = Node::Object()->simParameters;
  Bool lesOn = simParams->lesOn;
  Random vel_random(simParams->randomSeed);
  int lesReduceTemp = lesOn && simParams->lesReduceTemp;
  BigReal tempFactor = lesReduceTemp ? 1.0 / simParams->lesFactor : 1.0;

  kbT = Temp*BOLTZMANN;
  int count=0;
  int totalAtoms = inAtoms.size();
  for(i=0;i<totalAtoms;i++)
  {
    Real atomMs=inAtoms[i].mass;

    if (atomMs <= 0.) {
      kbToverM = 0.;
    } else {
      /*
       * lesOn is not supported in MEM_OPT_VERSION, so the original assignment
       * is simplified. --Chao Mei
       */
      //kbToverM = sqrt(kbT *
        //( lesOn && structure->get_fep_type(aid) ? tempFactor : 1.0 ) /
        //                  atomMs );
      kbToverM = sqrt(kbT * 1.0 / atomMs);
    }
    for (randnum=0.0, j=0; j<12; j++)
    {
      randnum += vel_random.uniform();
    }

    randnum -= 6.0;

    inAtoms[i].velocity.x = randnum*kbToverM;

    for (randnum=0.0, j=0; j<12; j++)
    {
      randnum += vel_random.uniform();
    }

    randnum -= 6.0;

    inAtoms[i].velocity.y = randnum*kbToverM;

    for (randnum=0.0, j=0; j<12; j++)
    {
      randnum += vel_random.uniform();
    }

    randnum -= 6.0;
    
    inAtoms[i].velocity.z = randnum*kbToverM;
  }
}
#endif

//----------------------------------------------------------------------
// This should only be called on node 0.
//----------------------------------------------------------------------
FullAtomList *WorkDistrib::createAtomLists(void)
{
  int i;
  StringList *current;	//  Pointer used to retrieve configuration items
  CProxy_Node nd(CkpvAccess(BOCclass_group).node);
  Node *node = nd.ckLocalBranch();
  PatchMap *patchMap = PatchMap::Object();
  CProxy_PatchMgr pm(CkpvAccess(BOCclass_group).patchMgr);
  PatchMgr *patchMgr = pm.ckLocalBranch();
  SimParameters *params = node->simParameters;
  Molecule *molecule = node->molecule;
  PDB *pdb = node->pdb;

  int numPatches = patchMap->numPatches();
  int numAtoms = pdb->num_atoms();

  Vector *positions = new Position[numAtoms];
  pdb->get_all_positions(positions);

  Vector *velocities = new Velocity[numAtoms];

  if ( params->initialTemp < 0.0 ) {
    Bool binvels=FALSE;

    //  Reading the veolcities from a PDB
    current = node->configList->find("velocities");

    if (current == NULL) {
      current = node->configList->find("binvelocities");
      binvels = TRUE;
    }

    if (!binvels) {
      velocities_from_PDB(current->data, velocities, numAtoms);
    }
    else {
      velocities_from_binfile(current->data, velocities, numAtoms);
    }
  }
  else {
    // Random velocities for a given temperature
    random_velocities(params->initialTemp, molecule, velocities, numAtoms);
  }

  //  If COMMotion == no, remove center of mass motion
  if (!(params->comMove)) {
    remove_com_motion(velocities, molecule, numAtoms);
  }

  FullAtomList *atoms = new FullAtomList[numPatches];

  const Lattice lattice = params->lattice;

    if ( params->staticAtomAssignment ) {
      FullAtomList sortAtoms;
      for ( i=0; i < numAtoms; i++ ) {
        HydrogenGroupID &h = molecule->hydrogenGroup[i];
        if ( ! h.isMP ) continue;
        FullAtom a;
        a.id = i;
        a.migrationGroupSize = h.isMP ? h.atomsInMigrationGroup : 0;
        a.position = positions[h.atomID];
        sortAtoms.add(a);
      } 
      int *order = new int[sortAtoms.size()];
      for ( i=0; i < sortAtoms.size(); i++ ) {
        order[i] = i;
      } 
      int *breaks = new int[numPatches];
      sortAtomsForPatches(order,breaks,sortAtoms.begin(),
                        sortAtoms.size(),numAtoms,
                        patchMap->gridsize_c(),
                        patchMap->gridsize_b(),
                        patchMap->gridsize_a());

      i = 0;
      for ( int pid = 0; pid < numPatches; ++pid ) {
        int iend = breaks[pid];
        for ( ; i<iend; ++i ) {
          FullAtom &sa = sortAtoms[order[i]];
          int mgs = sa.migrationGroupSize;
/*
CkPrintf("patch %d (%d %d %d) has group %d atom %d size %d at %.2f %.2f %.2f\n",
          pid, patchMap->index_a(pid), patchMap->index_b(pid), 
          patchMap->index_c(pid), order[i], sa.id, mgs,
          sa.position.x, sa.position.y, sa.position.z);
*/
          for ( int k=0; k<mgs; ++k ) {
            HydrogenGroupID &h = molecule->hydrogenGroup[sa.id + k];
            int aid = h.atomID;
            FullAtom a;
            a.id = aid;
            a.position = positions[aid];
            a.velocity = velocities[aid];
            a.vdwType = molecule->atomvdwtype(aid);
            a.status = molecule->getAtoms()[aid].status;
            a.langevinParam = molecule->langevin_param(aid);
            a.hydrogenGroupSize = h.isGP ? h.atomsInGroup : 0;
            a.migrationGroupSize = h.isMP ? h.atomsInMigrationGroup : 0;
            if(params->rigidBonds != RIGID_NONE) {
              a.rigidBondLength = molecule->rigid_bond_length(aid);
            }else{
              a.rigidBondLength = 0.0;
            }
            atoms[pid].add(a);
          }
        }
CkPrintf("patch %d (%d %d %d) has %d atoms\n",
          pid, patchMap->index_a(pid), patchMap->index_b(pid), 
          patchMap->index_c(pid), atoms[pid].size());
      }
      delete [] order;
      delete [] breaks;
    } else
    {
    // split atoms into patches based on migration group and position
    int aid, pid=0;
    for(i=0; i < numAtoms; i++)
      {
      // Assign atoms to patches without splitting hydrogen groups.
      // We know that the hydrogenGroup array is sorted with group parents
      // listed first.  Thus, only change the pid if an atom is a group parent.
      HydrogenGroupID &h = molecule->hydrogenGroup[i];
      aid = h.atomID;
      FullAtom a;
      a.id = aid;
      a.position = positions[aid];
      a.velocity = velocities[aid];
      a.vdwType = molecule->atomvdwtype(aid);
      a.status = molecule->getAtoms()[aid].status;
      a.langevinParam = molecule->langevin_param(aid);
      a.hydrogenGroupSize = h.isGP ? h.atomsInGroup : 0;
      a.migrationGroupSize = h.isMP ? h.atomsInMigrationGroup : 0;
      if(params->rigidBonds != RIGID_NONE) {
        a.rigidBondLength = molecule->rigid_bond_length(aid);
      }else{
        a.rigidBondLength = 0.0;
      }
      if (h.isMP) {
	pid = patchMap->assignToPatch(positions[aid],lattice);
      } // else: don't change pid
      atoms[pid].add(a);
      }
    }

  delete [] positions;
  delete [] velocities;

  for(i=0; i < numPatches; i++)
  {
    ScaledPosition center(0.5*(patchMap->min_a(i)+patchMap->max_a(i)),
			  0.5*(patchMap->min_b(i)+patchMap->max_b(i)),
			  0.5*(patchMap->min_c(i)+patchMap->max_c(i)));

    int n = atoms[i].size();
    FullAtom *a = atoms[i].begin();
    int j;
//Modifications for alchemical fep
    Bool alchFepOn = params->alchFepOn;
    Bool alchThermIntOn = params->alchThermIntOn;
//fepe
    Bool lesOn = params->lesOn;
  
    Bool pairInteractionOn = params->pairInteractionOn;

    Bool pressureProfileTypes = (params->pressureProfileAtomTypes > 1);

    Transform mother_transform;
    for(j=0; j < n; j++)
    {
      int aid = a[j].id;

      a[j].nonbondedGroupSize = 0;  // must be set based on coordinates

      a[j].atomFixed = molecule->is_atom_fixed(aid) ? 1 : 0;
      a[j].fixedPosition = a[j].position;

      if ( a[j].migrationGroupSize ) {
       if ( a[j].migrationGroupSize != a[j].hydrogenGroupSize ) {
            Position pos = a[j].position;
            int mgs = a[j].migrationGroupSize;
            int c = 1;
            for ( int k=a[j].hydrogenGroupSize; k<mgs;
                                k+=a[j+k].hydrogenGroupSize ) {
              pos += a[j+k].position;
              ++c;
            }
            pos *= 1./c;
            mother_transform = a[j].transform;  // should be 0,0,0
            pos = lattice.nearest(pos,center,&mother_transform);
            a[j].position = lattice.apply_transform(a[j].position,mother_transform);
            a[j].transform = mother_transform;
       } else {
        a[j].position = lattice.nearest(
		a[j].position, center, &(a[j].transform));
        mother_transform = a[j].transform;
       }
      } else {
        a[j].position = lattice.apply_transform(a[j].position,mother_transform);
        a[j].transform = mother_transform;
      }

      a[j].mass = molecule->atommass(aid);
      a[j].charge = molecule->atomcharge(aid);

//Modifications for alchemical fep
      if ( alchFepOn || alchThermIntOn || lesOn || pairInteractionOn || pressureProfileTypes) {
        a[j].partition = molecule->get_fep_type(aid);
      } 
      else {
        a[j].partition = 0;
      }
//fepe

    }

    int size, allfixed, k;
    for(j=0; j < n; j+=size) {
      size = a[j].hydrogenGroupSize;
      if ( ! size ) {
        NAMD_bug("Mother atom with hydrogenGroupSize of 0!");
      }
      allfixed = 1;
      for ( k = 0; k < size; ++k ) {
        allfixed = ( allfixed && (a[j+k].atomFixed) );
      }
      for ( k = 0; k < size; ++k ) {
        a[j+k].groupFixed = allfixed ? 1 : 0;
      }
    }

    if ( params->outputPatchDetails ) {
      int patchId = i;
      int numAtomsInPatch = n;
      int numFixedAtomsInPatch = 0;
      int numAtomsInFixedGroupsInPatch = 0;
      for(j=0; j < n; j++) {
        numFixedAtomsInPatch += ( a[j].atomFixed ? 1 : 0 );
        numAtomsInFixedGroupsInPatch += ( a[j].groupFixed ? 1 : 0 );
      }
      iout << "PATCH_DETAILS:"
           << " patch " << patchId
           << " atoms " << numAtomsInPatch
           << " fixed_atoms " << numFixedAtomsInPatch
           << " fixed_groups " << numAtomsInFixedGroupsInPatch
           << "\n" << endi;
    }
  }

  return atoms;

}

//----------------------------------------------------------------------
// This should only be called on node 0.
//----------------------------------------------------------------------
void WorkDistrib::createHomePatches(void)
{
  int i;
  PatchMap *patchMap = PatchMap::Object();
  CProxy_PatchMgr pm(CkpvAccess(BOCclass_group).patchMgr);
  PatchMgr *patchMgr = pm.ckLocalBranch();

  int numPatches = patchMap->numPatches();

  FullAtomList *atoms = createAtomLists();
    
#ifdef MEM_OPT_VERSION
/*  CProxy_Node nd(CkpvAccess(BOCclass_group).node);
  Node *node = nd.ckLocalBranch();
  node->molecule->delEachAtomSigs();
  node->molecule->delMassChargeSpace();
*/
#endif

  int maxAtoms = -1;
  int maxPatch = -1;
  for(i=0; i < numPatches; i++) {
    int numAtoms = atoms[i].size();
    if ( numAtoms > maxAtoms ) { maxAtoms = numAtoms; maxPatch = i; }
  }
  iout << iINFO << "LARGEST PATCH (" << maxPatch <<
	") HAS " << maxAtoms << " ATOMS\n" << endi;

  for(i=0; i < numPatches; i++)
  {
    if ( ! ( i % 100 ) )
    {
      DebugM(3,"Created " << i << " patches so far.\n");
    }

    patchMgr->createHomePatch(i,atoms[i]);
  }

  delete [] atoms;
}

void WorkDistrib::distributeHomePatches() {
  // ref BOC
  CProxy_Node nd(CkpvAccess(BOCclass_group).node);
  Node *node = nd.ckLocalBranch();
  CProxy_PatchMgr pm(CkpvAccess(BOCclass_group).patchMgr);
  PatchMgr *patchMgr = pm.ckLocalBranch();
  // ref singleton
  PatchMap *patchMap = PatchMap::Object();

  // Move patches to the proper node
  for(int i=0;i < patchMap->numPatches(); i++)
  {
    if (patchMap->node(i) != node->myid() )
    {
      DebugM(3,"patchMgr->movePatch("
	<< i << "," << patchMap->node(i) << ")\n");
      patchMgr->movePatch(i,patchMap->node(i));
    }
  }
  patchMgr->sendMovePatches();
}

void WorkDistrib::reinitAtoms() {

  PatchMap *patchMap = PatchMap::Object();
  CProxy_PatchMgr pm(CkpvAccess(BOCclass_group).patchMgr);
  PatchMgr *patchMgr = pm.ckLocalBranch();

  int numPatches = patchMap->numPatches();

  FullAtomList *atoms = createAtomLists();

  for(int i=0; i < numPatches; i++) {
    patchMgr->sendAtoms(i,atoms[i]);
  }

  delete [] atoms;

}


//----------------------------------------------------------------------

class PatchMapMsg : public CMessage_PatchMapMsg {
  public:
    char *patchMapData;
};

void WorkDistrib::sendPatchMap(void)
{
  if ( CkNumPes() == 1 ) {
    patchMapArrived = true;
    return;
  }

  //Automatically enable spanning tree
  CProxy_Node nd(CkpvAccess(BOCclass_group).node);
  Node *node = nd.ckLocalBranch();
  SimParameters *params = node->simParameters;
  if( ( PatchMap::Object()->numPatches() <= CkNumPes()/4
#ifdef NODEAWARE_PROXY_SPANNINGTREE 
      || CkNumPes() > CkNumNodes()
      ) && ( CkNumNodes() > 1
#endif
    ) && params->isSendSpanningTreeUnset() )
    ProxyMgr::Object()->setSendSpanning();

  int size = PatchMap::Object()->packSize();

  PatchMapMsg *mapMsg = new (size, 0) PatchMapMsg;

  PatchMap::Object()->pack(mapMsg->patchMapData);

  CProxy_WorkDistrib workProxy(thisgroup);
  workProxy[0].savePatchMap(mapMsg);
}

// saveMaps() is called when the map message is received
void WorkDistrib::savePatchMap(PatchMapMsg *msg)
{
  // Use a resend to forward messages before processing.  Otherwise the
  // map distribution is slow on many CPUs.  We need to use a tree
  // rather than a broadcast because some implementations of broadcast
  // generate a copy of the message on the sender for each recipient.
  // This is because MPI doesn't allow re-use of an outstanding buffer.

  if ( CkMyRank() ) patchMapArrived = true;

  if ( patchMapArrived && CkMyPe() ) {
    PatchMap::Object()->unpack(msg->patchMapData);

    //Automatically enable spanning tree
    CProxy_Node nd(CkpvAccess(BOCclass_group).node);
    Node *node = nd.ckLocalBranch();
    SimParameters *params = node->simParameters;
    if( ( PatchMap::Object()->numPatches() <= CkNumPes()/4
#ifdef NODEAWARE_PROXY_SPANNINGTREE 
        || CkNumPes() > CkNumNodes()
        ) && ( CkNumNodes() > 1
#endif
      ) && params->isSendSpanningTreeUnset() )
      ProxyMgr::Object()->setSendSpanning();
  }

  if ( patchMapArrived ) {
    if ( CkMyRank() + 1 < CkNodeSize(CkMyNode()) ) {
      ((CProxy_WorkDistrib(thisgroup))[CkMyPe()+1]).savePatchMap(msg);
    } else {
      delete msg;
    }
    return;
  }

  patchMapArrived = true;

  int self = CkMyNode();
  int range_begin = 0;
  int range_end = CkNumNodes();
  while ( self != range_begin ) {
    ++range_begin;
    int split = range_begin + ( range_end - range_begin ) / 2;
    if ( self < split ) { range_end = split; }
    else { range_begin = split; }
  }
  int send_near = self + 1;
  int send_far = send_near + ( range_end - send_near ) / 2;

  int pids[3];
  int npid = 0;
  if ( send_far < range_end ) pids[npid++] = CkNodeFirst(send_far);
  if ( send_near < send_far ) pids[npid++] = CkNodeFirst(send_near);
  pids[npid++] = CkMyPe();  // always send the message to ourselves
  CProxy_WorkDistrib(thisgroup).savePatchMap(msg,npid,pids);
}


class ComputeMapMsg : public CMessage_ComputeMapMsg {
  public:
    int nComputes;
    ComputeMap::ComputeData *computeMapData;
};

void WorkDistrib::sendComputeMap(void)
{
  if ( CkNumNodes() == 1 ) {
    computeMapArrived = true;
    ComputeMap::Object()->initPtrs();
    return;
  }

  int size = ComputeMap::Object()->numComputes();

  ComputeMapMsg *mapMsg = new (size, 0) ComputeMapMsg;

  mapMsg->nComputes = size;
  ComputeMap::Object()->pack(mapMsg->computeMapData);

  CProxy_WorkDistrib workProxy(thisgroup);
  workProxy[0].saveComputeMap(mapMsg);
}

// saveMaps() is called when the map message is received
void WorkDistrib::saveComputeMap(ComputeMapMsg *msg)
{
  // Use a resend to forward messages before processing.  Otherwise the
  // map distribution is slow on many CPUs.  We need to use a tree
  // rather than a broadcast because some implementations of broadcast
  // generate a copy of the message on the sender for each recipient.
  // This is because MPI doesn't allow re-use of an outstanding buffer.

  if ( CkMyRank() ) {
    NAMD_bug("WorkDistrib::saveComputeMap called on non-rank-zero pe");
  }

  if ( computeMapArrived && CkMyPe() ) {
    ComputeMap::Object()->unpack(msg->nComputes, msg->computeMapData);
  }

  if ( computeMapArrived ) {
    delete msg;
    ComputeMap::Object()->initPtrs();
    return;
  }

  computeMapArrived = true;

  int self = CkMyNode();
  int range_begin = 0;
  int range_end = CkNumNodes();
  while ( self != range_begin ) {
    ++range_begin;
    int split = range_begin + ( range_end - range_begin ) / 2;
    if ( self < split ) { range_end = split; }
    else { range_begin = split; }
  }
  int send_near = self + 1;
  int send_far = send_near + ( range_end - send_near ) / 2;

  int pids[3];
  int npid = 0;
  if ( send_far < range_end ) pids[npid++] = CkNodeFirst(send_far);
  if ( send_near < send_far ) pids[npid++] = CkNodeFirst(send_near);
  pids[npid++] = CkMyPe();  // always send the message to ourselves
  CProxy_WorkDistrib(thisgroup).saveComputeMap(msg,npid,pids);
}


//----------------------------------------------------------------------
void WorkDistrib::patchMapInit(void)
{
  PatchMap *patchMap = PatchMap::Object();
  CProxy_Node nd(CkpvAccess(BOCclass_group).node);
  Node *node = nd.ckLocalBranch();
  SimParameters *params = node->simParameters;
  Lattice lattice = params->lattice;

  BigReal patchSize = params->patchDimension;

  ScaledPosition xmin, xmax;

  double maxNumPatches = 1.e9;  // need to adjust fractional values
  if ( params->minAtomsPerPatch > 0 )
#ifndef MEM_OPT_VERSION
    maxNumPatches = node->pdb->num_atoms() / params->minAtomsPerPatch;
#else
    maxNumPatches = node->molecule->numAtoms / params->minAtomsPerPatch;
#endif

  DebugM(3,"Mapping patches\n");
  if ( lattice.a_p() && lattice.b_p() && lattice.c_p() ) {
    xmin = 0.;  xmax = 0.;
  }
  else if ( params->FMAOn || params->MSMOn ) {
  // Need to use full box for FMA to match NAMD 1.X results.
#if 0
    node->pdb->find_extremes(&(xmin.x),&(xmax.x),lattice.a_r());
    node->pdb->find_extremes(&(xmin.y),&(xmax.y),lattice.b_r());
    node->pdb->find_extremes(&(xmin.z),&(xmax.z),lattice.c_r());
#endif
    node->pdb->find_extremes(lattice);
    node->pdb->get_extremes(xmin, xmax);
#if 0
    printf("+++ center=%.4f %.4f %.4f\n",
        lattice.origin().x, lattice.origin().y, lattice.origin().z);
    printf("+++ xmin=%.4f  xmax=%.4f\n", xmin.x, xmax.x);
    printf("+++ ymin=%.4f  ymax=%.4f\n", xmin.y, xmax.y);
    printf("+++ zmin=%.4f  zmax=%.4f\n", xmin.z, xmax.z);
#endif
  // Otherwise, this allows a small number of stray atoms.
  }
  else {
#if 0
    node->pdb->find_extremes(&(xmin.x),&(xmax.x),lattice.a_r(),0.9);
    node->pdb->find_extremes(&(xmin.y),&(xmax.y),lattice.b_r(),0.9);
    node->pdb->find_extremes(&(xmin.z),&(xmax.z),lattice.c_r(),0.9);
#endif
    node->pdb->find_extremes(lattice, 0.9);
    node->pdb->get_extremes(xmin, xmax);
  }

#if 0
  BigReal origin_shift;
  origin_shift = lattice.a_r() * lattice.origin();
  xmin.x -= origin_shift;
  xmax.x -= origin_shift;
  origin_shift = lattice.b_r() * lattice.origin();
  xmin.y -= origin_shift;
  xmax.y -= origin_shift;
  origin_shift = lattice.c_r() * lattice.origin();
  xmin.z -= origin_shift;
  xmax.z -= origin_shift;
#endif

  // SimParameters default is -1 for unset
  int twoAwayX = params->twoAwayX;
  int twoAwayY = params->twoAwayY;
  int twoAwayZ = params->twoAwayZ;

  // SASA implementation is not compatible with twoAway patches
  if (params->LCPOOn) {
    if ( twoAwayX > 0 || twoAwayY > 0 || twoAwayZ > 0 ) {
      iout << iWARN << "Ignoring twoAway[XYZ] due to LCPO SASA implementation.\n" << endi;
    }
    twoAwayX = twoAwayY = twoAwayZ = 0;
  }

  // if you think you know what you're doing go right ahead
  if ( twoAwayX > 0 ) maxNumPatches = 1.e9;
  if ( twoAwayY > 0 ) maxNumPatches = 1.e9;
  if ( twoAwayZ > 0 ) maxNumPatches = 1.e9;
  if ( params->maxPatches > 0 ) {
      maxNumPatches = params->maxPatches;
      iout << iINFO << "LIMITING NUMBER OF PATCHES TO " <<
                                maxNumPatches << "\n" << endi;
  }

  int numpes = CkNumPes();
  SimParameters *simparam = Node::Object()->simParameters;
  if(simparam->simulateInitialMapping) {
    numpes = simparam->simulatedPEs;
    delete [] patchMap->nPatchesOnNode;
    patchMap->nPatchesOnNode = new int[numpes];
    memset(patchMap->nPatchesOnNode, 0, numpes*sizeof(int));	
  }

#ifdef NAMD_CUDA
  // for CUDA be sure there are more patches than pes

  int numPatches = patchMap->sizeGrid(
	xmin,xmax,lattice,patchSize,maxNumPatches,params->staticAtomAssignment,
	twoAwayX>0 ? 2 : 1, twoAwayY>0 ? 2 : 1, twoAwayZ>0 ? 2 : 1);
  if ( numPatches < numpes && twoAwayX < 0 ) {
    twoAwayX = 1;
    numPatches = patchMap->sizeGrid(
	xmin,xmax,lattice,patchSize,maxNumPatches,params->staticAtomAssignment,
	twoAwayX>0 ? 2 : 1, twoAwayY>0 ? 2 : 1, twoAwayZ>0 ? 2 : 1);
  }
  if ( numPatches < numpes && twoAwayY < 0 ) {
    twoAwayY = 1;
    numPatches = patchMap->sizeGrid(
	xmin,xmax,lattice,patchSize,maxNumPatches,params->staticAtomAssignment,
	twoAwayX>0 ? 2 : 1, twoAwayY>0 ? 2 : 1, twoAwayZ>0 ? 2 : 1);
  }
  if ( numPatches < numpes && twoAwayZ < 0 ) {
    twoAwayZ = 1;
    numPatches = patchMap->sizeGrid(
	xmin,xmax,lattice,patchSize,maxNumPatches,params->staticAtomAssignment,
	twoAwayX>0 ? 2 : 1, twoAwayY>0 ? 2 : 1, twoAwayZ>0 ? 2 : 1);
  }
  if ( numPatches < numpes ) {
    NAMD_die("CUDA-enabled NAMD requires at least one patch per process.");
  }
  if ( numPatches % numpes && numPatches <= 1.4 * numpes ) {
    int exactFit = numPatches - numPatches % numpes;
    int newNumPatches = patchMap->sizeGrid(
	xmin,xmax,lattice,patchSize,exactFit,params->staticAtomAssignment,
	twoAwayX>0 ? 2 : 1, twoAwayY>0 ? 2 : 1, twoAwayZ>0 ? 2 : 1);
    if ( newNumPatches == exactFit ) {
      iout << iINFO << "REDUCING NUMBER OF PATCHES TO IMPROVE LOAD BALANCE\n" << endi;
      maxNumPatches = exactFit;
    }
  }

  patchMap->makePatches(xmin,xmax,lattice,patchSize,maxNumPatches,params->staticAtomAssignment,
	twoAwayX>0 ? 2 : 1, twoAwayY>0 ? 2 : 1, twoAwayZ>0 ? 2 : 1);

#else

  int availPes = numpes;
  if ( params->noPatchesOnZero && numpes > 1 ) {
      availPes -= 1;
      if(params->noPatchesOnOne && numpes > 2)
        availPes -= 1;
  }
#ifdef MEM_OPT_VERSION
  if(params->noPatchesOnOutputPEs && numpes - params->numoutputprocs >2)
    {
      availPes -= params->numoutputprocs;
      if ( params->noPatchesOnZero && numpes > 1 && isOutputProcessor(0)){
	availPes++;
      }
      if ( params->noPatchesOnOne && numpes > 2 && isOutputProcessor(1)){
	availPes++;
      }
    }
#endif

  int numPatches = patchMap->sizeGrid(
	xmin,xmax,lattice,patchSize,1.e9,params->staticAtomAssignment,
	twoAwayX ? 2 : 1, twoAwayY ? 2 : 1, twoAwayZ ? 2 : 1);
  if ( ( numPatches > (0.3*availPes) || numPatches > maxNumPatches
       ) && twoAwayZ < 0 ) {
    twoAwayZ = 0;
    numPatches = patchMap->sizeGrid(
	xmin,xmax,lattice,patchSize,1.e9,params->staticAtomAssignment,
	twoAwayX ? 2 : 1, twoAwayY ? 2 : 1, twoAwayZ ? 2 : 1);
  }
  if ( ( numPatches > (0.6*availPes) || numPatches > maxNumPatches
       ) && twoAwayY < 0 ) {
    twoAwayY = 0;
    numPatches = patchMap->sizeGrid(
	xmin,xmax,lattice,patchSize,1.e9,params->staticAtomAssignment,
	twoAwayX ? 2 : 1, twoAwayY ? 2 : 1, twoAwayZ ? 2 : 1);
  }
  if ( ( numPatches > availPes || numPatches > maxNumPatches
       ) && twoAwayX < 0 ) {
    twoAwayX = 0;
    numPatches = patchMap->sizeGrid(
	xmin,xmax,lattice,patchSize,1.e9,params->staticAtomAssignment,
	twoAwayX ? 2 : 1, twoAwayY ? 2 : 1, twoAwayZ ? 2 : 1);
  }
  if ( numPatches > availPes && numPatches <= (1.4*availPes) && availPes <= maxNumPatches ) {
    int newNumPatches = patchMap->sizeGrid(
	xmin,xmax,lattice,patchSize,availPes,params->staticAtomAssignment,
	twoAwayX ? 2 : 1, twoAwayY ? 2 : 1, twoAwayZ ? 2 : 1);
    if ( newNumPatches <= availPes && numPatches <= (1.4*newNumPatches) ) {
      iout << iINFO << "REDUCING NUMBER OF PATCHES TO IMPROVE LOAD BALANCE\n" << endi;
      maxNumPatches = availPes;
    }
  }

  patchMap->makePatches(xmin,xmax,lattice,patchSize,maxNumPatches,params->staticAtomAssignment,
	twoAwayX ? 2 : 1, twoAwayY ? 2 : 1, twoAwayZ ? 2 : 1);

#endif

}


//----------------------------------------------------------------------
void WorkDistrib::assignNodeToPatch()
{
  PatchMap *patchMap = PatchMap::Object();
  int nNodes = Node::Object()->numNodes();
  SimParameters *simparam = Node::Object()->simParameters;
  if(simparam->simulateInitialMapping) {
	  nNodes = simparam->simulatedPEs;
  }

#if USE_TOPOMAP 
  TopoManager tmgr;
  int numPes = tmgr.getDimNX() * tmgr.getDimNY() * tmgr.getDimNZ();
  if (numPes > patchMap->numPatches() && (assignPatchesTopoGridRecBisection() > 0)) {
    CkPrintf ("Blue Gene/L topology partitioner finished successfully \n");
  }
  else
#endif
    if (nNodes > patchMap->numPatches())
      assignPatchesBitReversal();
    else
#if ! (CMK_BLUEGENEP | CMK_BLUEGENEL)
    if ( simparam->PMEOn )
      assignPatchesSpaceFillingCurve();	  
    else
#endif
      assignPatchesRecursiveBisection();
      // assignPatchesRoundRobin();
      // assignPatchesToLowestLoadNode();
  
  int *nAtoms = new int[nNodes];
  int numAtoms=0;
  int i;
  for(i=0; i < nNodes; i++)
    nAtoms[i] = 0;

  for(i=0; i < patchMap->numPatches(); i++)
  {
    //    iout << iINFO << "Patch " << i << " has " 
    //	 << patchMap->patch(i)->getNumAtoms() << " atoms and "
    //	 << patchMap->patch(i)->getNumAtoms() * 
    //            patchMap->patch(i)->getNumAtoms() 
    //	 << " pairs.\n" << endi;	 
#ifdef MEM_OPT_VERSION
      numAtoms += patchMap->numAtoms(i);
      nAtoms[patchMap->node(i)] += patchMap->numAtoms(i);	  
#else
    if (patchMap->patch(i)) {
      numAtoms += patchMap->patch(i)->getNumAtoms();
      nAtoms[patchMap->node(i)] += patchMap->patch(i)->getNumAtoms();	  
    }
#endif
  }

  if ( numAtoms != Node::Object()->molecule->numAtoms ) {
    for(i=0; i < nNodes; i++)
      iout << iINFO << nAtoms[i] << " atoms assigned to node " << i << "\n" << endi;
    iout << iINFO << "Assigned " << numAtoms << " atoms but expected " << Node::Object()->molecule->numAtoms << "\n" << endi;
    NAMD_die("Incorrect atom count in WorkDistrib::assignNodeToPatch\n");
  }

  delete [] nAtoms;
 
  //  PatchMap::Object()->printPatchMap();
}

//----------------------------------------------------------------------
// void WorkDistrib::assignPatchesSlices() 
// {
//   int pid; 
//   int assignedNode = 0;
//   PatchMap *patchMap = PatchMap::Object();
//   Node *node = CLocalBranch(Node, CkpvAccess(BOCclass_group).node);

//   int *numAtoms = new int[node->numNodes()];
//   for (int i=0; i<node->numNodes(); i++) {
//     numAtoms[i] = 0;
//   }

//   // Assign patch to node with least atoms assigned.
//   for(pid=0; pid < patchMap->numPatches(); pid++) {
//     assignedNode = 0;
//     for (i=1; i < node->numNodes(); i++) {
//       if (numAtoms[i] < numAtoms[assignedNode]) assignedNode = i;
//     }
//     patchMap->assignNode(pid, assignedNode);
//     numAtoms[assignedNode] += patchMap->patch(pid)->getNumAtoms();

//     /*
//     iout << iINFO << "Patch (" << pid << ") has " 
//       << patchMap->patch(pid)->getNumAtoms() 
//       << " atoms:  Assigned to Node(" << assignedNode << ")\n" 
//       << endi;
//     */
//   }

//   delete[] numAtoms;
// }

//----------------------------------------------------------------------
void WorkDistrib::assignPatchesToLowestLoadNode() 
{
  int pid; 
  int assignedNode = 0;
  PatchMap *patchMap = PatchMap::Object();
  CProxy_Node nd(CkpvAccess(BOCclass_group).node);
  Node *node = nd.ckLocalBranch();
  SimParameters *simParams = node->simParameters;
  int ncpus = node->numNodes();
  if(simParams->simulateInitialMapping) {
	  ncpus = simParams->simulatedPEs;
  }

  int *load = new int[ncpus];
  int *assignedNodes = new int[patchMap->numPatches()];
  for (int i=0; i<ncpus; i++) {
    load[i] = 0;
  }
  CkPrintf("assignPatchesToLowestLoadNode\n");
  int defaultNode = 0;
  if ( simParams->noPatchesOnZero && ncpus > 1 ){
    defaultNode = 1;
    if( simParams->noPatchesOnOne && ncpus > 2)
      defaultNode = 2;
  }
  // Assign patch to node with least atoms assigned.
  for(pid=0; pid < patchMap->numPatches(); pid++) {
    assignedNode = defaultNode;
    for (int i=assignedNode + 1; i < ncpus; i++) {
      if (load[i] < load[assignedNode]) assignedNode = i;
    }
    assignedNodes[pid] = assignedNode;
#ifdef MEM_OPT_VERSION
    load[assignedNode] += patchMap->numAtoms(pid) + 1;
#else
    load[assignedNode] += patchMap->patch(pid)->getNumAtoms() + 1;
#endif
  }

  delete[] load;
  sortNodesAndAssign(assignedNodes);
  delete[] assignedNodes;
}

//----------------------------------------------------------------------
void WorkDistrib::assignPatchesBitReversal() 
{
  int pid; 
  PatchMap *patchMap = PatchMap::Object();
  CProxy_Node nd(CkpvAccess(BOCclass_group).node);
  Node *node = nd.ckLocalBranch();
  SimParameters *simparam = node->simParameters;

  int ncpus = node->numNodes();
  if(simparam->simulateInitialMapping) {
	  ncpus = simparam->simulatedPEs;
  }
  int npatches = patchMap->numPatches();
  if ( ncpus <= npatches )
    NAMD_bug("WorkDistrib::assignPatchesBitReversal called improperly");

  // find next highest power of two
  int npow2 = 1;  int nbits = 0;
  while ( npow2 < ncpus ) { npow2 *= 2; nbits += 1; }

  // build bit reversal sequence
  SortableResizeArray<int> seq(ncpus);
  // avoid using node 0 (reverse of 0 is 0 so start at 1)
  int i = 1;
  for ( int icpu=0; icpu<(ncpus-1); ++icpu ) {
    int ri;
    for ( ri = ncpus; ri >= ncpus; ++i ) {
      ri = 0;
      int pow2 = 1;
      int rpow2 = npow2 / 2;
      for ( int j=0; j<nbits; ++j ) {
        ri += rpow2 * ( ( i / pow2 ) % 2 );
        pow2 *= 2;  rpow2 /= 2;
      }
    }
    seq[icpu] = ri;
  }

  // extract and sort patch locations
  sortNodesAndAssign(seq.begin());
  if ( ncpus > 2*npatches ) sortNodesAndAssign(seq.begin()+npatches, 1);
}

//----------------------------------------------------------------------
struct nodesort {
  int node;
  int a_total;
  int b_total;
  int c_total;
  int npatches;
  nodesort() : node(-1),a_total(0),b_total(0),c_total(0),npatches(0) { ; }
  int operator==(const nodesort &o) const {
    float a1 = ((float)a_total)/((float)npatches);
    float a2 = ((float)o.a_total)/((float)o.npatches);
    float b1 = ((float)b_total)/((float)npatches);
    float b2 = ((float)o.b_total)/((float)o.npatches);
    float c1 = ((float)c_total)/((float)npatches);
    float c2 = ((float)o.c_total)/((float)o.npatches);
    return ((a1 == a2) && (b1 == b2) && (c1 == c2));
  }
  int operator<(const nodesort &o) const {
    float a1 = ((float)a_total)/((float)npatches);
    float a2 = ((float)o.a_total)/((float)o.npatches);
    float b1 = ((float)b_total)/((float)npatches);
    float b2 = ((float)o.b_total)/((float)o.npatches);
    float c1 = ((float)c_total)/((float)npatches);
    float c2 = ((float)o.c_total)/((float)o.npatches);
    return ( (a1 < a2) || ((a1 == a2) && (b1 < b2)) ||
		((a1 == a2) && (b1 == b2) && (c1 < c2)) );
  }
};

void WorkDistrib::sortNodesAndAssign(int *assignedNode, int baseNodes) {
  // if baseNodes is zero (default) then set both nodes and basenodes
  // if baseNodes is nonzero then this is a second call to set basenodes only
  int i, pid; 
  PatchMap *patchMap = PatchMap::Object();
  int npatches = patchMap->numPatches();
  CProxy_Node nd(CkpvAccess(BOCclass_group).node);
  Node *node = nd.ckLocalBranch();
  int nnodes = node->numNodes();
  SimParameters *simparam = node->simParameters;
  if(simparam->simulateInitialMapping) {
	  nnodes = simparam->simulatedPEs;
  }

  ResizeArray<nodesort> allnodes(nnodes);
  for ( i=0; i < nnodes; ++i ) {
    allnodes[i].node = i;
  }
  for ( pid=0; pid<npatches; ++pid ) {
    // iout << pid << " " << assignedNode[pid] << "\n" << endi;
    allnodes[assignedNode[pid]].npatches++;
    allnodes[assignedNode[pid]].a_total += patchMap->index_a(pid);
    allnodes[assignedNode[pid]].b_total += patchMap->index_b(pid);
    allnodes[assignedNode[pid]].c_total += patchMap->index_c(pid);
  }
  SortableResizeArray<nodesort> usednodes(nnodes);
  usednodes.resize(0);
  for ( i=0; i < nnodes; ++i ) {
    if ( allnodes[i].npatches ) usednodes.add(allnodes[i]);
  }
  usednodes.sort();
  int nused = usednodes.size();
  int i2 = nused/2;
  for ( i=0; i < nnodes; ++i ) {
    if ( allnodes[i].npatches ) allnodes[usednodes[i2++].node].node = i;
    if ( i2 == nused ) i2 = 0;
  }

  for ( pid=0; pid<npatches; ++pid ) {
    // iout << pid << " " <<  allnodes[assignedNode[pid]].node << "\n" << endi;
    if ( ! baseNodes ) {
      patchMap->assignNode(pid, allnodes[assignedNode[pid]].node);	
    }
    patchMap->assignBaseNode(pid, allnodes[assignedNode[pid]].node);	
  }
}

//----------------------------------------------------------------------
void WorkDistrib::assignPatchesRoundRobin() 
{
  int pid; 
  PatchMap *patchMap = PatchMap::Object();
  CProxy_Node nd(CkpvAccess(BOCclass_group).node);
  Node *node = nd.ckLocalBranch();
  SimParameters *simparam = node->simParameters;
  int ncpus = node->numNodes();
  if(simparam->simulateInitialMapping) {
	  ncpus = simparam->simulatedPEs;
  }
  int *assignedNode = new int[patchMap->numPatches()];

  for(pid=0; pid < patchMap->numPatches(); pid++) {
    assignedNode[pid] = pid % ncpus;
  }

  sortNodesAndAssign(assignedNode);
  delete [] assignedNode;
}

//----------------------------------------------------------------------
void WorkDistrib::assignPatchesRecursiveBisection() 
{
  PatchMap *patchMap = PatchMap::Object();
  int *assignedNode = new int[patchMap->numPatches()];
  SimParameters *simParams = Node::Object()->simParameters;
  int numNodes = Node::Object()->numNodes();
  if(simParams->simulateInitialMapping) {
	  numNodes = simParams->simulatedPEs;
  }
  
  int usedNodes = numNodes;
  int unusedNodes = 0;
  CkPrintf("assignPatchesRecursiveBisection\n");
  if ( simParams->noPatchesOnZero && numNodes > 1 ){
    usedNodes -= 1;
    if(simParams->noPatchesOnOne && numNodes > 2)
      usedNodes -= 1;
  }
  unusedNodes = numNodes - usedNodes;
  RecBisection recBisec(usedNodes,PatchMap::Object());
  if ( recBisec.partition(assignedNode) ) {
    if ( unusedNodes !=0 ) {
      for ( int i=0; i<patchMap->numPatches(); ++i ) {
        assignedNode[i] += unusedNodes;
      }
    }
    sortNodesAndAssign(assignedNode);
    delete [] assignedNode;
  } else {
    //free the array here since a same array will be allocated
    //in assignPatchesToLowestLoadNode function, thus reducting
    //temporary memory usage
    delete [] assignedNode; 
    
    iout << iWARN 
	 << "WorkDistrib: Recursive bisection fails, "
	 << "invoking space-filling curve algorithm\n";
    assignPatchesSpaceFillingCurve();
  }
}

//----------------------------------------------------------------------
void WorkDistrib::assignPatchesSpaceFillingCurve() 
{
  PatchMap *patchMap = PatchMap::Object();
  int *assignedNode = new int[patchMap->numPatches()];
  int numNodes = Node::Object()->numNodes();
  SimParameters *simParams = Node::Object()->simParameters;
  if(simParams->simulateInitialMapping) {
	  numNodes = simParams->simulatedPEs;
  }
  int usedNodes = numNodes;
  int unusedNodes = 0;
  // make exclusion list
  if ( simParams->noPatchesOnZero && numNodes > 1 ){
    usedNodes -= 1;
    if(simParams->noPatchesOnOne && numNodes > 2)
      {
	usedNodes -= 1;
      }
  }  
#ifdef MEM_OPT_VERSION
  if(simParams->noPatchesOnOutputPEs && numNodes-simParams->numoutputprocs >2)
    {
      usedNodes -= simParams->numoutputprocs;
      if ( simParams->noPatchesOnZero && numNodes > 1 && isOutputProcessor(0)){
	usedNodes++;
      }
      if ( simParams->noPatchesOnOne && numNodes > 2 && isOutputProcessor(1)){
	usedNodes++;
      }
    }
#endif
  unusedNodes = numNodes - usedNodes;

  int numPatches = patchMap->numPatches();
  if ( numPatches < usedNodes )
    NAMD_bug("WorkDistrib::assignPatchesSpaceFillingCurve() called with more nodes than patches");

  ResizeArray<double> patchLoads(numPatches);
  SortableResizeArray<double> sortedLoads(numPatches);
  for ( int i=0; i<numPatches; ++i ) {
#ifdef MEM_OPT_VERSION
    double load = patchMap->numAtoms(i) + 10;      
#else
    double load = patchMap->patch(i)->getNumAtoms() + 10;
#endif
    patchLoads[i] = load;
    sortedLoads[i] = load;
  }
  sortedLoads.sort();

  // limit maxPatchLoad to adjusted average load per node
  double sumLoad = 0;
  double maxPatchLoad = 1;
  for ( int i=0; i<numPatches; ++i ) {
    double load = sortedLoads[i];
    double total = sumLoad + (numPatches-i) * load;
    if ( usedNodes * load > total ) break;
    sumLoad += load;
    maxPatchLoad = load;
  }
  double totalLoad = 0;
  for ( int i=0; i<numPatches; ++i ) {
    if ( patchLoads[i] > maxPatchLoad ) patchLoads[i] = maxPatchLoad;
    totalLoad += patchLoads[i];
  }
  if ( usedNodes * maxPatchLoad > totalLoad )
    NAMD_bug("algorithm failure in WorkDistrib::assignPatchesSpaceFillingCurve()");

  // walk through patches in space-filling curve
  sumLoad = 0;
  int node = 0;
  int skipIndex =0;
  int adim = patchMap->gridsize_a();
  int bdim = patchMap->gridsize_b();
  int cdim = patchMap->gridsize_c();
  int b = 0;
  int c = 0;
  int binc = 1;
  int cinc = 1;
  for ( int a = 0; a < adim; ++a ) {
    while ( b >= 0 && b < bdim ) {
      while ( c >= 0 && c < cdim ) {
	if ( simParams->noPatchesOnZero && numNodes > 1 && node==0) ++node;
	if ( simParams->noPatchesOnOne && numNodes > 2 && node==1) ++node;
#ifdef MEM_OPT_VERSION
	if(simParams->noPatchesOnOutputPEs)
	  { //advance through reserved pe list and skip output PEs
	    if(isOutputProcessor(node))
	      {
		if ( node+1 < numNodes )
		  {
		    CkPrintf("Patch Map Skipping %d to avoid Output PE\n",node);
		    ++node;
		  }
	      }
	  }
#endif
        int pid = patchMap->pid(a,b,c);
        assignedNode[pid] = node;
        sumLoad += patchLoads[pid];
        double targetLoad = (double)(node+1) / (double)numNodes;
        targetLoad *= totalLoad;
	//	CkPrintf("Patch Map Using node+1 %d numNodes %d sumload %f targetLoad %f\n",node+1, numNodes, sumLoad, targetLoad);
        if ( node+1 < numNodes && sumLoad >= targetLoad ) ++node;
        c += cinc;
      }
      cinc *= -1;  c += cinc;
      b += binc;
    }
    binc *= -1;  b += binc;
  }

  sortNodesAndAssign(assignedNode);
  delete [] assignedNode; 
}

//----------------------------------------------------------------------
void WorkDistrib::mapComputes(void)
{
  PatchMap *patchMap = PatchMap::Object();
  ComputeMap *computeMap = ComputeMap::Object();
  CProxy_Node nd(CkpvAccess(BOCclass_group).node);
  Node *node = nd.ckLocalBranch();

  DebugM(3,"Mapping computes\n");

  computeMap->allocateCids();

  // Handle full electrostatics
  if ( node->simParameters->fullDirectOn )
    mapComputeHomePatches(computeFullDirectType);
  if ( node->simParameters->FMAOn )
#ifdef DPMTA
    mapComputeHomePatches(computeDPMTAType);
#else
    NAMD_die("This binary does not include DPMTA (FMA).");
#endif
  if ( node->simParameters->PMEOn ) {
#ifdef DPME
    if ( node->simParameters->useDPME )
      mapComputeHomePatches(computeDPMEType);
    else {
      if (node->simParameters->useOptPME) {
	mapComputeHomePatches(optPmeType);
	if ( node->simParameters->pressureProfileEwaldOn )
	  mapComputeHomePatches(computeEwaldType);
      }
      else {
	mapComputeHomePatches(computePmeType);
	if ( node->simParameters->pressureProfileEwaldOn )
	  mapComputeHomePatches(computeEwaldType);
      }
    }
#else
    if (node->simParameters->useOptPME) {
      mapComputeHomePatches(optPmeType);
      if ( node->simParameters->pressureProfileEwaldOn )
	mapComputeHomePatches(computeEwaldType);
    }
    else {      
      mapComputePatch(computePmeType);
      if ( node->simParameters->pressureProfileEwaldOn )
	mapComputeHomePatches(computeEwaldType);
    }
#endif
  }

  if ( node->simParameters->globalForcesOn ) {
    DebugM(2,"adding ComputeGlobal\n");
    mapComputeHomePatches(computeGlobalType);
  }

  if ( node->simParameters->extForcesOn )
    mapComputeHomePatches(computeExtType);

  if ( node->simParameters->GBISserOn )
    mapComputeHomePatches(computeGBISserType);

  if ( node->simParameters->MsmSerialOn )
    mapComputeHomePatches(computeMsmSerialType);
#ifdef CHARM_HAS_MSA
  else if ( node->simParameters->MSMOn )
    mapComputeHomePatches(computeMsmMsaType);
#else
  else if ( node->simParameters->MSMOn )
    mapComputeHomePatches(computeMsmType);
#endif

#ifdef NAMD_CUDA
  mapComputeNode(computeNonbondedCUDAType);
  mapComputeHomeTuples(computeExclsType);
  mapComputePatch(computeSelfExclsType);
#endif

  mapComputeNonbonded();

  if ( node->simParameters->LCPOOn ) {
    mapComputeLCPO();
  }

  // If we're doing true pair interactions, no need for bonded terms.
  // But if we're doing within-group interactions, we do need them.
  if ( !node->simParameters->pairInteractionOn || 
      node->simParameters->pairInteractionSelf) { 
    mapComputeHomeTuples(computeBondsType);
    mapComputeHomeTuples(computeAnglesType);
    mapComputeHomeTuples(computeDihedralsType);
    mapComputeHomeTuples(computeImpropersType);
    mapComputeHomeTuples(computeCrosstermsType);
    mapComputePatch(computeSelfBondsType);
    mapComputePatch(computeSelfAnglesType);
    mapComputePatch(computeSelfDihedralsType);
    mapComputePatch(computeSelfImpropersType);
    mapComputePatch(computeSelfCrosstermsType);
  }

  if ( node->simParameters->drudeOn ) {
    mapComputeHomeTuples(computeTholeType);
    mapComputePatch(computeSelfTholeType);
    mapComputeHomeTuples(computeAnisoType);
    mapComputePatch(computeSelfAnisoType);
  }

  if ( node->simParameters->eFieldOn )
    mapComputePatch(computeEFieldType);
  /* BEGIN gf */
  if ( node->simParameters->mgridforceOn )
    mapComputePatch(computeGridForceType);
  /* END gf */
  if ( node->simParameters->stirOn )
    mapComputePatch(computeStirType);
  if ( node->simParameters->sphericalBCOn )
    mapComputePatch(computeSphericalBCType);
  if ( node->simParameters->cylindricalBCOn )
    mapComputePatch(computeCylindricalBCType);
  if ( node->simParameters->tclBCOn ) {
    mapComputeHomePatches(computeTclBCType);
  }
  if ( node->simParameters->constraintsOn )
    mapComputePatch(computeRestraintsType);
  if ( node->simParameters->consForceOn )
    mapComputePatch(computeConsForceType);
  if ( node->simParameters->consTorqueOn )
    mapComputePatch(computeConsTorqueType);

    // store the latest compute map
  SimParameters *simParams = Node::Object()->simParameters;
  if (simParams->storeComputeMap) {
    computeMap->saveComputeMap(simParams->computeMapFilename);
  }
    // override mapping decision
  if (simParams->loadComputeMap) {
    computeMap->loadComputeMap(simParams->computeMapFilename);
    CkPrintf("ComputeMap has been loaded from %s.\n", simParams->computeMapFilename);
  }
}

//----------------------------------------------------------------------
void WorkDistrib::mapComputeHomeTuples(ComputeType type)
{
  PatchMap *patchMap = PatchMap::Object();
  ComputeMap *computeMap = ComputeMap::Object();
  CProxy_Node nd(CkpvAccess(BOCclass_group).node);
  Node *node = nd.ckLocalBranch();

  int numNodes = node->numNodes();
  SimParameters *simparam = node->simParameters;
  if(simparam->simulateInitialMapping) {
	  numNodes = simparam->simulatedPEs;
  }

  char *isBaseNode = new char[numNodes];
  memset(isBaseNode,0,numNodes*sizeof(char));

  int numPatches = patchMap->numPatches();
  for(int j=0; j<numPatches; j++) {
    isBaseNode[patchMap->basenode(j)] = 1;
  }

  for(int i=0; i<numNodes; i++) {
    if ( isBaseNode[i] ) {
      computeMap->storeCompute(i,0,type);
    }
  }

  delete [] isBaseNode;
}

//----------------------------------------------------------------------
void WorkDistrib::mapComputeHomePatches(ComputeType type)
{
  PatchMap *patchMap = PatchMap::Object();
  ComputeMap *computeMap = ComputeMap::Object();
  CProxy_Node nd(CkpvAccess(BOCclass_group).node);
  Node *node = nd.ckLocalBranch();

  int numNodes = node->numNodes();
  SimParameters *simparam = node->simParameters;
  if(simparam->simulateInitialMapping) {
	  numNodes = simparam->simulatedPEs;
  }

  for(int i=0; i<numNodes; i++) {
    if ( patchMap->numPatchesOnNode(i) ) {
      computeMap->storeCompute(i,0,type);
    }
  }
}

//----------------------------------------------------------------------
void WorkDistrib::mapComputePatch(ComputeType type)
{
  PatchMap *patchMap = PatchMap::Object();
  ComputeMap *computeMap = ComputeMap::Object();

  PatchID i;
  ComputeID cid;

  for(i=0; i<patchMap->numPatches(); i++)
  {
    cid=computeMap->storeCompute(patchMap->node(i),1,type);
    computeMap->newPid(cid,i);
    patchMap->newCid(i,cid);
  }

}


//----------------------------------------------------------------------
void WorkDistrib::mapComputeNode(ComputeType type)
{
  PatchMap *patchMap = PatchMap::Object();
  ComputeMap *computeMap = ComputeMap::Object();

  PatchID i;
  ComputeID cid;

  int ncpus = CkNumPes();
  SimParameters *simparam = Node::Object()->simParameters;
  if(simparam->simulateInitialMapping) {
	  ncpus = simparam->simulatedPEs;
  }

  for(int i=0; i<ncpus; i++) {
    computeMap->storeCompute(i,0,type);
  }

}


//----------------------------------------------------------------------
void WorkDistrib::mapComputeNonbonded(void)
{
  // For each patch, create 1 electrostatic object for self-interaction.
  // Then create 1 for each 1-away and 2-away neighbor which has a larger
  // pid.

  PatchMap *patchMap = PatchMap::Object();
  ComputeMap *computeMap = ComputeMap::Object();
  CProxy_Node nd(CkpvAccess(BOCclass_group).node);
  Node *node = nd.ckLocalBranch();
  SimParameters *simParams = Node::Object()->simParameters;
  int ncpus = CkNumPes();
  int nodesize = CkMyNodeSize();
  if(simParams->simulateInitialMapping) {
	  ncpus = simParams->simulatedPEs;
	  nodesize = simParams->simulatedNodeSize;
  }

  PatchID oneAway[PatchMap::MaxOneOrTwoAway];
  PatchID oneAwayDownstream[PatchMap::MaxOneOrTwoAway];
  int oneAwayTrans[PatchMap::MaxOneOrTwoAway];

  PatchID i;
  ComputeID cid;
  int numNeighbors;
  int j;
  double partScaling = 1.0;
  if ( ncpus < patchMap->numPatches() ) {
    partScaling = ((double)ncpus) / ((double)patchMap->numPatches());
  }

  for(i=0; i<patchMap->numPatches(); i++) // do the self 
  {

   int numPartitions = 1;
   if ( simParams->ldBalancer == LDBAL_HYBRID ) {
#ifdef  MEM_OPT_VERSION    
    int64 numFixed = patchMap->numFixedAtoms(i);
    int64 numAtoms = patchMap->numAtoms(i);
#else
    int64 numFixed = patchMap->patch(i)->getNumFixedAtoms();  // avoid overflow
    int64 numAtoms = patchMap->patch(i)->getNumAtoms();
#endif

    int divide = node->simParameters->numAtomsSelf;
    if (divide > 0) {
      numPartitions = (int) ( partScaling * ( 0.5 +
        (numAtoms*numAtoms-numFixed*numFixed) / (double)(2*divide*divide) ) );
    }
    if (numPartitions < 1) numPartitions = 1;
    if ( numPartitions > node->simParameters->maxSelfPart )
			numPartitions = node->simParameters->maxSelfPart;
    // self-interaction
    DebugM(4,"Mapping " << numPartitions << " ComputeNonbondedSelf objects for patch " << i << "\n");
//    iout <<"Self numPartitions = " <<numPartitions <<" numAtoms " <<numAtoms <<std::endl;
   }

    for(int partition=0; partition < numPartitions; partition++)
    {
      cid=computeMap->storeCompute(patchMap->node(i),1,
				   computeNonbondedSelfType,
				   partition,numPartitions);
      computeMap->newPid(cid,i);
      patchMap->newCid(i,cid);
    }
  }

  for(int p1=0; p1 <patchMap->numPatches(); p1++) // do the pairs
  {
    // this only returns half of neighbors, which is what we want
    numNeighbors=patchMap->oneOrTwoAwayNeighbors(p1,oneAway,oneAwayDownstream,oneAwayTrans);
    for(j=0;j<numNeighbors;j++)
    {
	int p2 = oneAway[j];
	int dsp = oneAwayDownstream[j];

      int numPartitions = 1;
      if ( simParams->ldBalancer == LDBAL_HYBRID ) {
#ifdef  MEM_OPT_VERSION        
        int64 numAtoms1 = patchMap->numAtoms(p1);
        int64 numAtoms2 = patchMap->numAtoms(p2);
        int64 numFixed1 = patchMap->numFixedAtoms(p1);
        int64 numFixed2 = patchMap->numFixedAtoms(p2);
#else
        int64 numAtoms1 = patchMap->patch(p1)->getNumAtoms();
        int64 numAtoms2 = patchMap->patch(p2)->getNumAtoms();
        int64 numFixed1 = patchMap->patch(p1)->getNumFixedAtoms();
        int64 numFixed2 = patchMap->patch(p2)->getNumFixedAtoms();
#endif


        const int t2 = oneAwayTrans[j];
        const int adim = patchMap->gridsize_a();
        const int bdim = patchMap->gridsize_b();
        const int cdim = patchMap->gridsize_c();
        const int nax = patchMap->numaway_a();  // 1 or 2
        const int nay = patchMap->numaway_b();  // 1 or 2
        const int naz = patchMap->numaway_c();  // 1 or 2
        const int ia1 = patchMap->index_a(p1);
        const int ia2 = patchMap->index_a(p2) + adim * Lattice::offset_a(t2);
        const int ib1 = patchMap->index_b(p1);
        const int ib2 = patchMap->index_b(p2) + bdim * Lattice::offset_b(t2);
        const int ic1 = patchMap->index_c(p1);
        const int ic2 = patchMap->index_c(p2) + cdim * Lattice::offset_c(t2);

        if ( abs(ia2-ia1) > nax ||
             abs(ib2-ib1) > nay ||
             abs(ic2-ic1) > naz )
          NAMD_bug("Bad patch distance in WorkDistrib::mapComputeNonbonded");

	int distance = 3;
 	if ( ia1 == ia2 ) --distance;
 	else if ( ia1 == ia2 + nax - 1 ) --distance;
 	else if ( ia1 + nax - 1 == ia2 ) --distance;
 	if ( ib1 == ib2 ) --distance;
 	else if ( ib1 == ib2 + nay - 1 ) --distance;
 	else if ( ib1 + nay - 1 == ib2 ) --distance;
 	if ( ic1 == ic2 ) --distance;
 	else if ( ic1 == ic2 + naz - 1 ) --distance;
 	else if ( ic1 + naz - 1 == ic2 ) --distance;
	int divide = 0;
	if ( distance == 0 ) {
	  divide = node->simParameters->numAtomsSelf2;
        } else if (distance == 1) {
	  divide = node->simParameters->numAtomsPair;
	} else {
	  divide = node->simParameters->numAtomsPair2;
	}
	if (divide > 0) {
          numPartitions = (int) ( partScaling * ( 0.5 +
	    (numAtoms1*numAtoms2-numFixed1*numFixed2)/(double)(divide*divide) ) );
	}
        if ( numPartitions < 1 ) numPartitions = 1;
        if ( numPartitions > node->simParameters->maxPairPart )
			numPartitions = node->simParameters->maxPairPart;
//	if ( numPartitions > 1 ) iout << "Mapping " << numPartitions << " ComputeNonbondedPair objects for patches " << p1 << "(" << numAtoms1 << ") and " << p2 << "(" << numAtoms2 << ")\n" << endi;
      }
		for(int partition=0; partition < numPartitions; partition++)
		{
		  cid=computeMap->storeCompute( patchMap->basenode(dsp),
			2,computeNonbondedPairType,partition,numPartitions);
		  computeMap->newPid(cid,p1);
		  computeMap->newPid(cid,p2,oneAwayTrans[j]);
		  patchMap->newCid(p1,cid);
		  patchMap->newCid(p2,cid);
		}
    }
  }
}

//----------------------------------------------------------------------
void WorkDistrib::mapComputeLCPO(void) {
  //iterate over all needed objects

  PatchMap *patchMap = PatchMap::Object();
  ComputeMap *computeMap = ComputeMap::Object();
  CProxy_Node nd(CkpvAccess(BOCclass_group).node);
  Node *node = nd.ckLocalBranch();
  SimParameters *simParams = Node::Object()->simParameters;
  int ncpus = CkNumPes();
  int nodesize = CkMyNodeSize();
  const int maxPatches = 8;

  int numPatchesInOctet;
  PatchID patchesInOctet[maxPatches];
  int oneAwayTrans[maxPatches];

  //partitioned after 1st timestep
  int numPartitions = 1;

  PatchID i;
  ComputeID cid;

  // one octet per patch
  for(i=0; i<patchMap->numPatches(); i++) {
    numPatchesInOctet =
        patchMap->getPatchesInOctet(i, patchesInOctet, oneAwayTrans);

		for(int partition=0; partition < numPartitions; partition++) {
      cid=computeMap->storeCompute(patchMap->node(i),
          numPatchesInOctet,
				  computeLCPOType,
				  partition,
          numPartitions);
      for (int p = 0; p < numPatchesInOctet; p++) {
        computeMap->newPid(cid, patchesInOctet[p], oneAwayTrans[p]);
      }
      for (int p = 0; p < numPatchesInOctet; p++) {
        patchMap->newCid(patchesInOctet[p],cid);
      }
    } // for partitions 
  } // for patches
} // mapComputeLCPO

//----------------------------------------------------------------------
void WorkDistrib::messageEnqueueWork(Compute *compute) {
  LocalWorkMsg *msg = compute->localWorkMsg;
  int seq = compute->sequence();
  int gbisPhase = compute->getGBISPhase();

  if ( seq < 0 ) {
    NAMD_bug("compute->sequence() < 0 in WorkDistrib::messageEnqueueWork");
  } else {
    SET_PRIORITY(msg,seq,compute->priority());
  }

  msg->compute = compute; // pointer is valid since send is to local Pe
  int type = compute->type();
  int cid = compute->cid;

  CProxy_WorkDistrib wdProxy(CkpvAccess(BOCclass_group).workDistrib);
  switch ( type ) {
  case computeExclsType:
  case computeSelfExclsType:
    wdProxy[CkMyPe()].enqueueExcls(msg);
    break;
  case computeBondsType:
  case computeSelfBondsType:
    wdProxy[CkMyPe()].enqueueBonds(msg);
    break;
  case computeAnglesType:
  case computeSelfAnglesType:
    wdProxy[CkMyPe()].enqueueAngles(msg);
    break;
  case computeDihedralsType:
  case computeSelfDihedralsType:
    wdProxy[CkMyPe()].enqueueDihedrals(msg);
    break;
  case computeImpropersType:
  case computeSelfImpropersType:
    wdProxy[CkMyPe()].enqueueImpropers(msg);
    break;
  case computeTholeType:
  case computeSelfTholeType:
    wdProxy[CkMyPe()].enqueueThole(msg);
    break;
  case computeAnisoType:
  case computeSelfAnisoType:
    wdProxy[CkMyPe()].enqueueAniso(msg);
    break;
  case computeCrosstermsType:
  case computeSelfCrosstermsType:
    wdProxy[CkMyPe()].enqueueCrossterms(msg);
    break;
  case computeLCPOType:
    wdProxy[CkMyPe()].enqueueLCPO(msg);
    break;
  case computeNonbondedSelfType:
    switch ( seq % 2 ) {
    case 0:
      //wdProxy[CkMyPe()].enqueueSelfA(msg);
      switch ( gbisPhase ) {
         case 1:
           wdProxy[CkMyPe()].enqueueSelfA1(msg);
           break;
         case 2:
           wdProxy[CkMyPe()].enqueueSelfA2(msg);
           break;
         case 3:
           wdProxy[CkMyPe()].enqueueSelfA3(msg);
           break;
      }
      break;
    case 1:
      //wdProxy[CkMyPe()].enqueueSelfB(msg);
      switch ( gbisPhase ) {
         case 1:
           wdProxy[CkMyPe()].enqueueSelfB1(msg);
           break;
         case 2:
           wdProxy[CkMyPe()].enqueueSelfB2(msg);
           break;
         case 3:
           wdProxy[CkMyPe()].enqueueSelfB3(msg);
           break;
      }
      break;
    default:
      NAMD_bug("WorkDistrib::messageEnqueueSelf case statement error!");
    }
    break;
  case computeNonbondedPairType:
    switch ( seq % 2 ) {
    case 0:
      //wdProxy[CkMyPe()].enqueueWorkA(msg);
      switch ( gbisPhase ) {
         case 1:
           wdProxy[CkMyPe()].enqueueWorkA1(msg);
           break;
         case 2:
           wdProxy[CkMyPe()].enqueueWorkA2(msg);
           break;
         case 3:
           wdProxy[CkMyPe()].enqueueWorkA3(msg);
           break;
      }
      break;
    case 1:
      //wdProxy[CkMyPe()].enqueueWorkB(msg);
      switch ( gbisPhase ) {
         case 1:
           wdProxy[CkMyPe()].enqueueWorkB1(msg);
           break;
         case 2:
           wdProxy[CkMyPe()].enqueueWorkB2(msg);
           break;
         case 3:
           wdProxy[CkMyPe()].enqueueWorkB3(msg);
           break;
      }
      break;
    case 2:
      wdProxy[CkMyPe()].enqueueWorkC(msg);
      break;
    default:
      NAMD_bug("WorkDistrib::messageEnqueueWork case statement error!");
    }
    break;
  case computeNonbondedCUDAType:
#ifdef NAMD_CUDA
//     CkPrintf("WorkDistrib[%d]::CUDA seq=%d phase=%d\n", CkMyPe(), seq, gbisPhase);
    //wdProxy[CkMyPe()].enqueueCUDA(msg);
    switch ( gbisPhase ) {
       case 1:
         wdProxy[CkMyPe()].enqueueCUDA(msg);
         break;
       case 2:
         wdProxy[CkMyPe()].enqueueCUDAP2(msg);
         break;
       case 3:
         wdProxy[CkMyPe()].enqueueCUDAP3(msg);
         break;
    }
#else
    msg->compute->doWork();  traceUserEvent(eventMachineProgress);  CmiMachineProgressImpl();
#endif
    break;
  case computePmeType:
    // CkPrintf("PME %d %d %x\n", CkMyPe(), seq, compute->priority());
    wdProxy[CkMyPe()].enqueuePme(msg);
    break;
  case optPmeType:
    // CkPrintf("PME %d %d %x\n", CkMyPe(), seq, compute->priority());
#ifdef NAMD_CUDA
    wdProxy[CkMyPe()].enqueuePme(msg);
#else
    msg->compute->doWork();  traceUserEvent(eventMachineProgress);  CmiMachineProgressImpl();
#endif
    break;
  default:
    wdProxy[CkMyPe()].enqueueWork(msg);
  }
}

void WorkDistrib::enqueueWork(LocalWorkMsg *msg) {
  msg->compute->doWork();  traceUserEvent(eventMachineProgress);  CmiMachineProgressImpl();
  if ( msg->compute->localWorkMsg != msg )
    NAMD_bug("WorkDistrib LocalWorkMsg recycling failed!");
}

void WorkDistrib::enqueueExcls(LocalWorkMsg *msg) {
  msg->compute->doWork();  traceUserEvent(eventMachineProgress);  CmiMachineProgressImpl();
  if ( msg->compute->localWorkMsg != msg )
    NAMD_bug("WorkDistrib LocalWorkMsg recycling failed!");
}

void WorkDistrib::enqueueBonds(LocalWorkMsg *msg) {
  msg->compute->doWork();  traceUserEvent(eventMachineProgress);  CmiMachineProgressImpl();
  if ( msg->compute->localWorkMsg != msg )
    NAMD_bug("WorkDistrib LocalWorkMsg recycling failed!");
}

void WorkDistrib::enqueueAngles(LocalWorkMsg *msg) {
  msg->compute->doWork();  traceUserEvent(eventMachineProgress);  CmiMachineProgressImpl();
  if ( msg->compute->localWorkMsg != msg )
    NAMD_bug("WorkDistrib LocalWorkMsg recycling failed!");
}

void WorkDistrib::enqueueDihedrals(LocalWorkMsg *msg) {
  msg->compute->doWork();  traceUserEvent(eventMachineProgress);  CmiMachineProgressImpl();
  if ( msg->compute->localWorkMsg != msg )
    NAMD_bug("WorkDistrib LocalWorkMsg recycling failed!");
}

void WorkDistrib::enqueueImpropers(LocalWorkMsg *msg) {
  msg->compute->doWork();  traceUserEvent(eventMachineProgress);  CmiMachineProgressImpl();
  if ( msg->compute->localWorkMsg != msg )
    NAMD_bug("WorkDistrib LocalWorkMsg recycling failed!");
}

void WorkDistrib::enqueueThole(LocalWorkMsg *msg) {
  msg->compute->doWork();  traceUserEvent(eventMachineProgress);  CmiMachineProgressImpl();
  if ( msg->compute->localWorkMsg != msg )
    NAMD_bug("WorkDistrib LocalWorkMsg recycling failed!");
}

void WorkDistrib::enqueueAniso(LocalWorkMsg *msg) {
  msg->compute->doWork();  traceUserEvent(eventMachineProgress);  CmiMachineProgressImpl();
  if ( msg->compute->localWorkMsg != msg )
    NAMD_bug("WorkDistrib LocalWorkMsg recycling failed!");
}

void WorkDistrib::enqueueCrossterms(LocalWorkMsg *msg) {
  msg->compute->doWork();  traceUserEvent(eventMachineProgress);  CmiMachineProgressImpl();
  if ( msg->compute->localWorkMsg != msg )
    NAMD_bug("WorkDistrib LocalWorkMsg recycling failed!");
}

void WorkDistrib::enqueuePme(LocalWorkMsg *msg) {
  msg->compute->doWork();  traceUserEvent(eventMachineProgress);  CmiMachineProgressImpl();
  if ( msg->compute->localWorkMsg != msg )
    NAMD_bug("WorkDistrib LocalWorkMsg recycling failed!");
}

void WorkDistrib::enqueueLCPO(LocalWorkMsg *msg) {
  msg->compute->doWork();
  if ( msg->compute->localWorkMsg != msg )
    NAMD_bug("WorkDistrib LocalWorkMsg recycling failed!");
}
void WorkDistrib::enqueueSelfA1(LocalWorkMsg *msg) {
  msg->compute->doWork();  traceUserEvent(eventMachineProgress);  CmiMachineProgressImpl();
  if ( msg->compute->localWorkMsg != msg )
    NAMD_bug("WorkDistrib LocalWorkMsg recycling failed!");
}
void WorkDistrib::enqueueSelfA2(LocalWorkMsg *msg) {
  msg->compute->doWork();  traceUserEvent(eventMachineProgress);  CmiMachineProgressImpl();
  if ( msg->compute->localWorkMsg != msg )
    NAMD_bug("WorkDistrib LocalWorkMsg recycling failed!");
}
void WorkDistrib::enqueueSelfA3(LocalWorkMsg *msg) {
  msg->compute->doWork();  traceUserEvent(eventMachineProgress);  CmiMachineProgressImpl();
  if ( msg->compute->localWorkMsg != msg )
    NAMD_bug("WorkDistrib LocalWorkMsg recycling failed!");
}

void WorkDistrib::enqueueSelfB1(LocalWorkMsg *msg) {
  msg->compute->doWork();  traceUserEvent(eventMachineProgress);  CmiMachineProgressImpl();
  if ( msg->compute->localWorkMsg != msg )
    NAMD_bug("WorkDistrib LocalWorkMsg recycling failed!");
}
void WorkDistrib::enqueueSelfB2(LocalWorkMsg *msg) {
  msg->compute->doWork();  traceUserEvent(eventMachineProgress);  CmiMachineProgressImpl();
  if ( msg->compute->localWorkMsg != msg )
    NAMD_bug("WorkDistrib LocalWorkMsg recycling failed!");
}
void WorkDistrib::enqueueSelfB3(LocalWorkMsg *msg) {
  msg->compute->doWork();  traceUserEvent(eventMachineProgress);  CmiMachineProgressImpl();
  if ( msg->compute->localWorkMsg != msg )
    NAMD_bug("WorkDistrib LocalWorkMsg recycling failed!");
}

void WorkDistrib::enqueueWorkA1(LocalWorkMsg *msg) {
  msg->compute->doWork();  traceUserEvent(eventMachineProgress);  CmiMachineProgressImpl();
  if ( msg->compute->localWorkMsg != msg )
    NAMD_bug("WorkDistrib LocalWorkMsg recycling failed!");
}
void WorkDistrib::enqueueWorkA2(LocalWorkMsg *msg) {
  msg->compute->doWork();  traceUserEvent(eventMachineProgress);  CmiMachineProgressImpl();
  if ( msg->compute->localWorkMsg != msg )
    NAMD_bug("WorkDistrib LocalWorkMsg recycling failed!");
}
void WorkDistrib::enqueueWorkA3(LocalWorkMsg *msg) {
  msg->compute->doWork();  traceUserEvent(eventMachineProgress);  CmiMachineProgressImpl();
  if ( msg->compute->localWorkMsg != msg )
    NAMD_bug("WorkDistrib LocalWorkMsg recycling failed!");
}

void WorkDistrib::enqueueWorkB1(LocalWorkMsg *msg) {
  msg->compute->doWork();  traceUserEvent(eventMachineProgress);  CmiMachineProgressImpl();
  if ( msg->compute->localWorkMsg != msg )
    NAMD_bug("WorkDistrib LocalWorkMsg recycling failed!");
}
void WorkDistrib::enqueueWorkB2(LocalWorkMsg *msg) {
  msg->compute->doWork();  traceUserEvent(eventMachineProgress);  CmiMachineProgressImpl();
  if ( msg->compute->localWorkMsg != msg )
    NAMD_bug("WorkDistrib LocalWorkMsg recycling failed!");
}
void WorkDistrib::enqueueWorkB3(LocalWorkMsg *msg) {
  msg->compute->doWork();  traceUserEvent(eventMachineProgress);  CmiMachineProgressImpl();
  if ( msg->compute->localWorkMsg != msg )
    NAMD_bug("WorkDistrib LocalWorkMsg recycling failed!");
}



void WorkDistrib::enqueueWorkC(LocalWorkMsg *msg) {
  msg->compute->doWork();  traceUserEvent(eventMachineProgress);  CmiMachineProgressImpl();
  if ( msg->compute->localWorkMsg != msg )
    NAMD_bug("WorkDistrib LocalWorkMsg recycling failed!");
}

void WorkDistrib::enqueueCUDA(LocalWorkMsg *msg) {
  msg->compute->doWork();  traceUserEvent(eventMachineProgress);  CmiMachineProgressImpl();
  // ComputeNonbondedCUDA *c = msg->compute;
  // if ( c->localWorkMsg != msg && c->localWorkMsg2 != msg )
  //   NAMD_bug("WorkDistrib LocalWorkMsg recycling failed!");
}
void WorkDistrib::enqueueCUDAP2(LocalWorkMsg *msg) {
  msg->compute->doWork();  traceUserEvent(eventMachineProgress);  CmiMachineProgressImpl();
}
void WorkDistrib::enqueueCUDAP3(LocalWorkMsg *msg) {
  msg->compute->doWork();  traceUserEvent(eventMachineProgress);  CmiMachineProgressImpl();
}


//**********************************************************************
//
//			FUNCTION velocities_from_PDB
//
//   INPUTS:
//      v - Array of vectors to populate
//	filename - name of the PDB filename to read in
//
//	This function reads in a set of initial velocities from a
//      PDB file.  It places the velocities into the array of Vectors
//      passed to it.
//
//***********************************************************************/

void WorkDistrib::velocities_from_PDB(char *filename, 
				      Vector *v, int totalAtoms)
{
  PDB *v_pdb;		//  PDB info from velocity PDB
  int i;

  //  Read the PDB
  v_pdb = new PDB(filename);
  if ( v_pdb == NULL )
  {
    NAMD_die("memory allocation failed in Node::velocities_from_PDB");
  }

  //  Make sure the number of velocities read in matches
  //  the number of atoms we have
  if (v_pdb->num_atoms() != totalAtoms)
  {
    char err_msg[129];

    sprintf(err_msg, "FOUND %d COORDINATES IN VELOCITY PDB!!",
	    v_pdb->num_atoms());

    NAMD_die(err_msg);
  }

  //  Get the entire list of atom info and loop through
  //  them assigning the velocity vector for each one
  v_pdb->get_all_positions(v);

  for (i=0; i<totalAtoms; i++)
  {
    v[i].x *= PDBVELINVFACTOR;
    v[i].y *= PDBVELINVFACTOR;
    v[i].z *= PDBVELINVFACTOR;
  }

  delete v_pdb;
}
//		END OF FUNCTION velocities_from_PDB

//**********************************************************************
//
// 			FUNCTION velocities_from_binfile
//
//    INPUTS:
// 	fname - File name to write velocities to
//	n - Number of atoms in system
//	vels - Array of velocity vectors
//					
//	This function writes out the velocities in binary format.  This is
//     done to preserve accuracy between restarts of namd.
//
//**********************************************************************

void WorkDistrib::velocities_from_binfile(char *fname, Vector *vels, int n)
{
  read_binary_file(fname,vels,n);
}
//               END OF FUNCTION velocities_from_binfile

//**********************************************************************
//
//			FUNCTION random_velocities
//
//   INPUTS:
//	v - array of vectors to populate
//	Temp - Temperature to acheive
//
//	This function assigns a random velocity distribution to a
//   simulation to achieve a desired initial temperature.  The method
//   used here was stolen from the program X-PLOR.
//
//**********************************************************************

void WorkDistrib::random_velocities(BigReal Temp,Molecule *structure,
				    Vector *v, int totalAtoms)
{
  int i, j;		//  Loop counter
  BigReal kbT;		//  Boltzman constant * Temp
  BigReal randnum;	//  Random number from -6.0 to 6.0
  BigReal kbToverM;	//  sqrt(Kb*Temp/Mass)
  SimParameters *simParams = Node::Object()->simParameters;
  Bool lesOn = simParams->lesOn;
  Random vel_random(simParams->randomSeed);

  int lesReduceTemp = lesOn && simParams->lesReduceTemp;
  BigReal tempFactor = lesReduceTemp ? 1.0 / simParams->lesFactor : 1.0;

  kbT = Temp*BOLTZMANN;

  //  Loop through all the atoms and assign velocities in
  //  the x, y and z directions for each one
  for (i=0; i<totalAtoms; i++)
  {
    if (structure->atommass(i) <= 0.) {
      kbToverM = 0.;
    } else {
      kbToverM = sqrt(kbT *
        ( lesOn && structure->get_fep_type(i) ? tempFactor : 1.0 ) /
			  structure->atommass(i) );
    }

    //  The following comment was stolen from X-PLOR where
    //  the following section of code was adapted from.
    
    //  This section generates a Gaussian random
    //  deviate of 0.0 mean and standard deviation RFD for
    //  each of the three spatial dimensions.
    //  The algorithm is a "sum of uniform deviates algorithm"
    //  which may be found in Abramowitz and Stegun,
    //  "Handbook of Mathematical Functions", pg 952.
    for (randnum=0.0, j=0; j<12; j++)
    {
      randnum += vel_random.uniform();
    }

    randnum -= 6.0;

    v[i].x = randnum*kbToverM;

    for (randnum=0.0, j=0; j<12; j++)
    {
      randnum += vel_random.uniform();
    }

    randnum -= 6.0;

    v[i].y = randnum*kbToverM;

    for (randnum=0.0, j=0; j<12; j++)
    {
      randnum += vel_random.uniform();
    }

    randnum -= 6.0;
    
    v[i].z = randnum*kbToverM;
  }
}
/*			END OF FUNCTION random_velocities		*/

//**********************************************************************
//
//			FUNCTION remove_com_motion
//
//   INPUTS:
//	vel - Array of initial velocity vectors
//
//	This function removes the center of mass motion from a molecule.
//
//**********************************************************************

void WorkDistrib::remove_com_motion(Vector *vel, Molecule *structure, int n)
{
  Vector mv(0,0,0);		//  Sum of (mv)_i
  BigReal totalMass=0; 	//  Total mass of system
  int i;			//  Loop counter

  //  Loop through and compute the net momentum
  for (i=0; i<n; i++)
  {
    BigReal mass = structure->atommass(i);
    mv += mass * vel[i];
    totalMass += mass;
  }

  mv /= totalMass;

  iout << iINFO << "REMOVING COM VELOCITY "
	<< ( PDBVELFACTOR * mv ) << "\n" << endi;

  for (i=0; i<n; i++) { vel[i] -= mv; }

}
/*			END OF FUNCTION remove_com_motion		*/

#if USE_TOPOMAP 

//Specifically designed for BGL and other 3d Tori architectures
//Partition Torus and Patch grid together using recursive bisection.
int WorkDistrib::assignPatchesTopoGridRecBisection() {
  
  PatchMap *patchMap = PatchMap::Object();
  int *assignedNode = new int[patchMap->numPatches()];
  int numNodes = Node::Object()->numNodes();
  SimParameters *simParams = Node::Object()->simParameters;
  if(simParams->simulateInitialMapping) {
	  numNodes = simParams->simulatedPEs;
  }

  int usedNodes = numNodes;
  CkPrintf("assignPatchesTopoGridRecBisection\n");
  if ( simParams->noPatchesOnZero && numNodes > 1 ) {
    usedNodes -= 1;
    if ( simParams->noPatchesOnOne && numNodes > 2 )
      usedNodes -= 1;
  }
  RecBisection recBisec(patchMap->numPatches(), PatchMap::Object());
  
  int xsize = 0, ysize = 0, zsize = 0;
  
  // Right now assumes a T*** (e.g. TXYZ) mapping
  TopoManager tmgr;
  xsize = tmgr.getDimNX();
  ysize = tmgr.getDimNY();
  zsize = tmgr.getDimNZ();
  
  //Fix to not assign patches to processor 0
  int rc = recBisec.partitionProcGrid(xsize, ysize, zsize, assignedNode);
 
  delete [] assignedNode;

  return rc;
}
#endif

#include "WorkDistrib.def.h"

