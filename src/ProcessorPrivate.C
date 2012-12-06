/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include "ProcessorPrivate.h"
#include "Debug.h"
#include "InfoStream.h"

/*
 * Variable Definitions
 */

// Instance Variables that maintain singletonness of classes

CkpvDeclare(AtomMap*, AtomMap_instance);
CkpvDeclare(BroadcastMgr*, BroadcastMgr_instance);
CkpvDeclare(CollectionMaster*, CollectionMaster_instance);
CkpvDeclare(CollectionMgr*, CollectionMgr_instance);
CkpvDeclare(LdbCoordinator*, LdbCoordinator_instance);
CkpvDeclare(Node*, Node_instance);
CkpvDeclare(PatchMap*, PatchMap_instance);
CkpvDeclare(PatchMgr*, PatchMgr_instance);
CkpvDeclare(ProxyMgr*, ProxyMgr_instance);
CkpvDeclare(ReductionMgr*, ReductionMgr_instance);

#ifdef PROCTRACE_DEBUG
CkpvDeclare(DebugFileTrace*, DebugFileTrace_instance);
#endif

// Other static variables

CkpvDeclare(PatchMgr*, PatchMap_patchMgr);
CkpvDeclare(BOCgroup, BOCclass_group);
CkpvDeclare(Communicate*, comm);
CkpvDeclare(Sync*, Sync_instance);
CkpvDeclare(infostream, iout_obj);

/*
 * Initialization Function to be called on every processor
 */

void ProcessorPrivateInit(void)
{
  CkpvInitialize(AtomMap*, AtomMap_instance);
  CkpvAccess(AtomMap_instance) = 0;
  CkpvInitialize(BroadcastMgr*, BroadcastMgr_instance);
  CkpvAccess(BroadcastMgr_instance) = 0;
  CkpvInitialize(CollectionMaster*, CollectionMaster_instance);
  CkpvAccess(CollectionMaster_instance) = 0;
  CkpvInitialize(CollectionMgr*, CollectionMgr_instance);
  CkpvAccess(CollectionMgr_instance) = 0;
  CkpvInitialize(LdbCoordinator*, LdbCoordinator_instance);
  CkpvAccess(LdbCoordinator_instance) = 0;
  CkpvInitialize(Node*, Node_instance);
  CkpvAccess(Node_instance) = 0;

  CkpvInitialize(PatchMap*, PatchMap_instance);
  CkpvAccess(PatchMap_instance) = 0;
  CkpvInitialize(PatchMgr*, PatchMgr_instance);
  CkpvAccess(PatchMgr_instance) = 0;
  CkpvInitialize(ProxyMgr*, ProxyMgr_instance);
  CkpvAccess(ProxyMgr_instance) = 0;
  CkpvInitialize(ReductionMgr*, ReductionMgr_instance);
  CkpvAccess(ReductionMgr_instance) = 0;
  CkpvInitialize(PatchMgr*, PatchMap_patchMgr);
  CkpvAccess(PatchMap_patchMgr) = 0;
  CkpvInitialize(BOCgroup, BOCclass_group);
  CkpvInitialize(Communicate*, comm);
  CkpvAccess(comm) = 0;
  CkpvInitialize(Sync*, Sync_instance);
  CkpvAccess(Sync_instance) = 0;
  CkpvInitialize(infostream, iout_obj);

#ifdef PROCTRACE_DEBUG
  CkpvInitialize(DebugFileTrace*, DebugFileTrace_instance);
  CkpvAccess(DebugFileTrace_instance) = 0;
#endif

}

