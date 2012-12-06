/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef _WORKDISTRIB_H
#define _WORKDISTRIB_H

#include "charm++.h"

#include "main.h"

#include "NamdTypes.h"
#include "BOCgroup.h"
#include "ComputeMap.h"
#include "WorkDistrib.decl.h"

class Node;
class Compute;
class Molecule;

// For Compute objects to enqueue themselves when ready to compute
class LocalWorkMsg : public CMessage_LocalWorkMsg
{
public:
  Compute *compute;
};

enum { maxPatchDepends = 126 };

class PatchMapMsg;
class ComputeMapMsg;
class ComputeMapChangeMsg;

class WorkDistrib : public BOCclass
{
public:
  WorkDistrib();
  ~WorkDistrib(void);

  // static void messageMovePatchDone();
  // void movePatchDone();

  static void messageEnqueueWork(Compute *);
  void enqueueWork(LocalWorkMsg *msg);
  void enqueueExcls(LocalWorkMsg *msg);
  void enqueueBonds(LocalWorkMsg *msg);
  void enqueueAngles(LocalWorkMsg *msg);
  void enqueueDihedrals(LocalWorkMsg *msg);
  void enqueueImpropers(LocalWorkMsg *msg);
  void enqueueThole(LocalWorkMsg *msg);  // Drude model
  void enqueueAniso(LocalWorkMsg *msg);  // Drude model
  void enqueueCrossterms(LocalWorkMsg *msg);
  void enqueuePme(LocalWorkMsg *msg);
  void enqueueSelfA1(LocalWorkMsg *msg);
  void enqueueSelfA2(LocalWorkMsg *msg);
  void enqueueSelfA3(LocalWorkMsg *msg);
  void enqueueSelfB1(LocalWorkMsg *msg);
  void enqueueSelfB2(LocalWorkMsg *msg);
  void enqueueSelfB3(LocalWorkMsg *msg);
  void enqueueWorkA1(LocalWorkMsg *msg);
  void enqueueWorkA2(LocalWorkMsg *msg);
  void enqueueWorkA3(LocalWorkMsg *msg);
  void enqueueWorkB1(LocalWorkMsg *msg);
  void enqueueWorkB2(LocalWorkMsg *msg);
  void enqueueWorkB3(LocalWorkMsg *msg);
  void enqueueWorkC(LocalWorkMsg *msg);
  void enqueueCUDA(LocalWorkMsg *msg);
  void enqueueCUDAP2(LocalWorkMsg *msg);
  void enqueueCUDAP3(LocalWorkMsg *msg);
  void enqueueLCPO(LocalWorkMsg *msg);

  void mapComputes(void);
  void sendPatchMap(void);
  void sendComputeMap(void);
  void saveComputeMapChanges(int,CkGroupID);
  void recvComputeMapChanges(ComputeMapChangeMsg *);
  void doneSaveComputeMap(CkReductionMsg *);

  FullAtomList *createAtomLists(void);
  void createHomePatches(void);
  void distributeHomePatches(void);

  void reinitAtoms(void);
  void patchMapInit(void);
  void assignNodeToPatch(void);

  void savePatchMap(PatchMapMsg *msg);
  void saveComputeMap(ComputeMapMsg *msg);
  inline void setPatchMapArrived(bool s) {patchMapArrived=s;}

#ifdef MEM_OPT_VERSION
  void fillAtomListForOnePatch(int pid, FullAtomList &alist);
  void random_velocities_parallel(BigReal Temp,InputAtomList &inAtoms);
#endif

private:
  void mapComputeNonbonded(void);
  void mapComputeLCPO(void);
  void mapComputeNode(ComputeType);
  void mapComputeHomePatches(ComputeType);
  void mapComputeHomeTuples(ComputeType);
  void mapComputePatch(ComputeType);
  void assignPatchesToLowestLoadNode(void);
  void assignPatchesRecursiveBisection(void);
  void assignPatchesRoundRobin(void);
  void assignPatchesSpaceFillingCurve(void);
  void assignPatchesBitReversal(void);
  int  assignPatchesTopoGridRecBisection();

  void sortNodesAndAssign(int *assignedNode, int baseNodes = 0);
  void velocities_from_PDB(char *filename, 
			   Vector *v, int totalAtoms);
  void velocities_from_binfile(char *fname, Vector *vels, int n);
  void random_velocities(BigReal Temp, Molecule *structure,
			 Vector *v, int totalAtoms);
  void remove_com_motion(Vector *vel, Molecule *structure, int n);

  bool patchMapArrived;
  bool computeMapArrived;

  int saveComputeMapReturnEP;
  CkGroupID saveComputeMapReturnChareID;
};

#endif /* WORKDISTRIB_H */

