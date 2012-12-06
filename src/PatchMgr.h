/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*****************************************************************************
 * $Source: /Projects/namd2/cvsroot/namd2/src/PatchMgr.h,v $
 * $Author: jlai7 $
 * $Date: 2012/11/27 21:13:19 $
 * $Revision: 1.1028 $
 *****************************************************************************/

#ifndef PATCHMGR_H
#define PATCHMGR_H

#include "charm++.h"

#include "NamdTypes.h"
#include "SortedArray.h"
#include "HomePatch.h"
#include "HomePatchList.h"
#include "BOCgroup.h"
#include "Migration.h"
#include "MigrateAtomsMsg.h"
#include "PatchMgr.decl.h"

#if USE_TOPOMAP 
#include "TopoManager.h"
#endif

class HomePatch;

class MovePatchesMsg : public CMessage_MovePatchesMsg {
public:
    NodeID  fromNodeID;
    PatchID pid;
    FullAtomList atom;

    MovePatchesMsg(void) { ; }

    MovePatchesMsg(PatchID n, FullAtomList a) : pid(n), atom(a)
    {
      fromNodeID = CkMyPe();
    }

  // pack and unpack functions
  static void* pack(MovePatchesMsg *msg);
  static MovePatchesMsg* unpack(void *ptr);
};

class MoveAtomMsg : public CMessage_MoveAtomMsg {
public:
  int atomid;
  int moveto;
  Vector coord;
};

class MoveAllByMsg : public CMessage_MoveAllByMsg {
public:
  Vector offset;
};

class SetLatticeMsg : public CMessage_SetLatticeMsg {
public:
  Lattice lattice;
};

// PatchMgr creates and manages homepatches. There exist one instance of 
// PatchMgr on each node (derived from Charm++ Group).  // That is, when a new operator causes creation of one instance on each node. 
// In addition to creation of homepatches, it handles the atom redistribution
// at the end of each cycle (i.e., atoms can move from patch to patch at the
// cycle boundaries).
struct MovePatch 
{
    MovePatch(PatchID p=-1, NodeID n=-1) : nodeID(n), pid(p) {};
    ~MovePatch() {};

    NodeID nodeID;
    PatchID pid;

    int operator<(MovePatch m) {
      return ( nodeID < m.nodeID );
    }

    int operator==(MovePatch m) {
      return ( nodeID == m.nodeID );
    }
};

typedef SortedArray<MovePatch> MovePatchList;
typedef ResizeArrayIter<MovePatch> MovePatchListIter;

class PatchMgr : public BOCclass
{

public:
  PatchMgr();
  ~PatchMgr();

  static PatchMgr* Object() { return CkpvAccess(PatchMgr_instance); }
  
  void createHomePatch(PatchID pid, FullAtomList a);

  //atomCnt is the number of atoms patch pid has
  void preCreateHomePatch(PatchID pid, int atomCnt);

  void movePatch(PatchID, NodeID);
  void sendMovePatches();
  void recvMovePatches(MovePatchesMsg *msg);

  void sendAtoms(PatchID pid, FullAtomList &a);
  void recvAtoms(MovePatchesMsg *msg);

  // void ackMovePatches(AckMovePatchesMsg *msg);

  HomePatch *homePatch(PatchID pid) {
     return homePatches.find(HomePatchElem(pid))->patch;
  } 

  // void sendMigrationMsg(PatchID, MigrationInfo);
  void sendMigrationMsgs(PatchID, MigrationInfo*, int);
  // void recvMigrateAtoms(MigrateAtomsMsg *);
  void recvMigrateAtomsCombined(MigrateAtomsCombinedMsg *);

  void moveAtom(MoveAtomMsg *msg);
  void moveAllBy(MoveAllByMsg *msg);
  void setLattice(SetLatticeMsg *msg);


private:
  friend class PatchMap;
  PatchMap *patchMap;

  int numAllPatches;
  int numHomePatches;

  // an array of patch pointers residing on this node
  HomePatchList homePatches;

  // an array of patches to move off this node
  MovePatchList move;
  int ackMovePending;

  // data for combining migration messages
  MigrateAtomsCombinedMsg ** combineMigrationMsgs;
  ResizeArray<int> combineMigrationDestPes;
  int migrationCountdown;

public:
  void fillHomePatchAtomList(int patchId, FullAtomList *al){
      HomePatch *thisHomePatch = patchMap->homePatch(patchId);
      thisHomePatch->setAtomList(al);
  }
  void setHomePatchFixedAtomNum(int patchId, int numFixed){
      HomePatch *thisHomePatch = patchMap->homePatch(patchId);
      thisHomePatch->setNumFixedAtoms(numFixed);
  }

  void sendOneHomePatch(int patchId, int nodeId);
};



#endif /* PATCHMGR_H */

