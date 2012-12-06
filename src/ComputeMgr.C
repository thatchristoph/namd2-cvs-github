/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include "InfoStream.h"
#include "ProcessorPrivate.h"

//#define DEBUGM
#define MIN_DEBUG_LEVEL 1
#include "Debug.h"

#include "BOCgroup.h"
#include "ComputeMgr.decl.h"
#include "ComputeMgr.h"
#include "ProxyMgr.decl.h"
#include "ProxyMgr.h"

#include "Node.h"
#include "ComputeMap.h"
#include "PatchMap.h"
#include "PatchMap.inl"

#include "Compute.h"
#include "ComputeNonbondedUtil.h"
#include "ComputeNonbondedSelf.h"
#include "ComputeNonbondedPair.h"
#include "ComputeNonbondedCUDA.h"
#include "ComputeAngles.h"
#include "ComputeDihedrals.h"
#include "ComputeImpropers.h"
#include "ComputeThole.h"
#include "ComputeAniso.h"
#include "ComputeCrossterms.h"
#include "ComputeBonds.h"
#include "ComputeNonbondedCUDAExcl.h"
#include "ComputeFullDirect.h"
#include "ComputeGlobal.h"
#include "ComputeGlobalMsgs.h"
#include "ComputeExt.h"
#include "ComputeGBISser.h"
#include "ComputeLCPO.h"
#include "ComputeMsmSerial.h"
#include "ComputeMsmMsa.h"
#include "ComputeMsm.h"
#include "ComputeDPMTA.h"
#include "ComputeDPME.h"
#include "ComputeDPMEMsgs.h"
#include "ComputePme.h"
#include "OptPme.h"
#include "ComputeEwald.h"
#include "ComputeEField.h"
/* BEGIN gf */
#include "ComputeGridForce.h"
/* END gf */
#include "ComputeStir.h"
#include "ComputeSphericalBC.h"
#include "ComputeCylindricalBC.h"
#include "ComputeTclBC.h"
#include "ComputeRestraints.h"
#include "ComputeConsForce.h"
#include "ComputeConsForceMsgs.h"
#include "WorkDistrib.h"

#include "LdbCoordinator.h"

/* include all of the specific masters we need here */
#include "FreeEnergyEnums.h"
#include "FreeEnergyAssert.h"
#include "FreeEnergyGroup.h"
#include "FreeEnergyVector.h"
#include "FreeEnergyRestrain.h"
#include "FreeEnergyRMgr.h"
#include "FreeEnergyLambda.h"
#include "FreeEnergyLambdMgr.h"

#include "GlobalMasterTest.h"
#include "GlobalMasterIMD.h"
#include "GlobalMasterTcl.h"
#include "GlobalMasterSMD.h"
#include "GlobalMasterTMD.h"
#include "GlobalMasterSymmetry.h"
#include "GlobalMasterEasy.h"
#include "GlobalMasterMisc.h"
#include "GlobalMasterFreeEnergy.h"
#include "GlobalMasterColvars.h"

ComputeMgr::ComputeMgr()
{
    CkpvAccess(BOCclass_group).computeMgr = thisgroup;
    computeGlobalObject = 0;
    computeDPMEObject = 0;
    computeEwaldObject = 0;
    computeNonbondedCUDAObject = 0;
    computeNonbondedWorkArrays = new ComputeNonbondedWorkArrays;
    skipSplitting = 0;
}

ComputeMgr::~ComputeMgr(void)
{
    delete computeNonbondedWorkArrays;
}

void ComputeMgr::updateComputes(int ep, CkGroupID chareID)
{
    updateComputesReturnEP = ep;
    updateComputesReturnChareID = chareID;
    updateComputesCount = CkNumPes();

    if (CkMyPe())
    {
        iout << iPE << iERRORF << "updateComputes signaled on wrong Pe!\n" << endi;
        CkExit();
        return;
    }

    CkStartQD(CkIndex_ComputeMgr::updateComputes2((CkQdMsg*)0),&thishandle);
}

void ComputeMgr::updateComputes2(CkQdMsg *msg)
{
    delete msg;

    CProxy_WorkDistrib wd(CkpvAccess(BOCclass_group).workDistrib);
    WorkDistrib  *workDistrib = wd.ckLocalBranch();
    workDistrib->saveComputeMapChanges(CkIndex_ComputeMgr::updateComputes3(),thisgroup);
}

void ComputeMgr::updateComputes3()
{
    if ( skipSplitting ) {
      CProxy_ComputeMgr(thisgroup).updateLocalComputes();
    } else {
      CProxy_ComputeMgr(thisgroup).splitComputes();
      skipSplitting = 1;
    }
}

void ComputeMgr::splitComputes()
{
  if ( ! CkMyRank() ) {
    ComputeMap *computeMap = ComputeMap::Object();
    const int nc = computeMap->numComputes();

    for (int i=0; i<nc; i++) {
      int nnp = computeMap->newNumPartitions(i);
      if ( nnp > 0 ) {
        if ( computeMap->numPartitions(i) != 1 ) {
          CkPrintf("Warning: unable to partition compute %d\n", i);
          computeMap->setNewNumPartitions(i,0);
          continue;
        }
        //CkPrintf("splitting compute %d by %d\n",i,nnp);
        computeMap->setNumPartitions(i,nnp);
        if (computeMap->newNode(i) == -1) {
          computeMap->setNewNode(i,computeMap->node(i));
        }
        for ( int j=1; j<nnp; ++j ) {
          int newcid = computeMap->cloneCompute(i,j);
          //CkPrintf("compute %d partition %d is %d\n",i,j,newcid);
        }
      }
    }
    computeMap->extendPtrs();
  }

  if (!CkMyPe())
  {
    CkStartQD(CkIndex_ComputeMgr::splitComputes2((CkQdMsg*)0), &thishandle);
  }
}

void ComputeMgr::splitComputes2(CkQdMsg *msg)
{
    delete msg;
    CProxy_ComputeMgr(thisgroup).updateLocalComputes();
}

void ComputeMgr::updateLocalComputes()
{
    ComputeMap *computeMap = ComputeMap::Object();
    CProxy_ProxyMgr pm(CkpvAccess(BOCclass_group).proxyMgr);
    ProxyMgr *proxyMgr = pm.ckLocalBranch();
    LdbCoordinator *ldbCoordinator = LdbCoordinator::Object();

     computeFlag.resize(0);

    const int nc = computeMap->numComputes();
    for (int i=0; i<nc; i++) {

        if ( computeMap->node(i) == CkMyPe() &&
             computeMap->newNumPartitions(i) > 1 ) {
           Compute *c = computeMap->compute(i);
           ldbCoordinator->Migrate(c->ldObjHandle,CkMyPe());
           delete c;
           computeMap->registerCompute(i,NULL);
           if ( computeMap->newNode(i) == CkMyPe() ) computeFlag.add(i); 
        } else
        if (computeMap->newNode(i) == CkMyPe() && computeMap->node(i) != CkMyPe())
        {
	    computeFlag.add(i);
            for (int n=0; n < computeMap->numPids(i); n++)
            {
                proxyMgr->createProxy(computeMap->pid(i,n));
            }
        }
        else if (computeMap->node(i) == CkMyPe() &&
                 (computeMap->newNode(i) != -1 && computeMap->newNode(i) != CkMyPe() ))
        {
            // CkPrintf("delete compute %d on pe %d\n",i,CkMyPe());
            delete computeMap->compute(i);
            computeMap->registerCompute(i,NULL);
        }
    }

    if (!CkMyPe())
    {
        CkStartQD(CkIndex_ComputeMgr::updateLocalComputes2((CkQdMsg*)0), &thishandle);
    }
}

void
ComputeMgr::updateLocalComputes2(CkQdMsg *msg)
{
    delete msg;
    CProxy_ComputeMgr(thisgroup).updateLocalComputes3();
}

void
ComputeMgr::updateLocalComputes3()
{
    ComputeMap *computeMap = ComputeMap::Object();
    CProxy_ProxyMgr pm(CkpvAccess(BOCclass_group).proxyMgr);
    ProxyMgr *proxyMgr = pm.ckLocalBranch();

    ProxyMgr::nodecount = 0;

    const int nc = computeMap->numComputes();

    if ( ! CkMyRank() ) {
      for (int i=0; i<nc; i++) {
        computeMap->setNewNumPartitions(i,0);
        if (computeMap->newNode(i) != -1) {
          computeMap->setNode(i,computeMap->newNode(i));
          computeMap->setNewNode(i,-1);
        }
      }
    }
 
    for(int i=0; i<computeFlag.size(); i++) createCompute(computeFlag[i], computeMap);
    computeFlag.clear();

    proxyMgr->removeUnusedProxies();

    if (!CkMyPe())
    {
        CkStartQD(CkIndex_ComputeMgr::updateLocalComputes4((CkQdMsg*)0), &thishandle);
    }
}

void
ComputeMgr::updateLocalComputes4(CkQdMsg *msg)
{
    delete msg;
    CProxy_ComputeMgr(thisgroup).updateLocalComputes5();

    // store the latest compute map
           SimParameters *simParams = Node::Object()->simParameters;
    if (simParams->storeComputeMap) {
      ComputeMap *computeMap = ComputeMap::Object();
      computeMap->saveComputeMap(simParams->computeMapFilename);
    }
}

#if 0
int firstphase = 1;
#endif

void
ComputeMgr::updateLocalComputes5()
{
    if ( ! CkMyRank() ) {
      ComputeMap::Object()->checkMap();
      PatchMap::Object()->checkMap();
    }

    // we always use the centralized building of spanning tree
    // distributed building of ST called in Node.C only
    if (proxySendSpanning || proxyRecvSpanning)
        ProxyMgr::Object()->buildProxySpanningTree2();

    // this code needs to be turned on if we want to
    // shift the creation of ST to the load balancer

#if 0
    if (proxySendSpanning || proxyRecvSpanning)
    {
        if (firstphase)
            ProxyMgr::Object()->buildProxySpanningTree2();
        else
            if (CkMyPe() == 0)
                ProxyMgr::Object()->sendSpanningTrees();

        firstphase = 0;
    }
#endif

    if (!CkMyPe())
        CkStartQD(CkIndex_ComputeMgr::doneUpdateLocalComputes(), &thishandle);
}

void ComputeMgr::doneUpdateLocalComputes()
{

//  if (!--updateComputesCount) {
    DebugM(4, "doneUpdateLocalComputes on Pe("<<CkMyPe()<<")\n");
    void *msg = CkAllocMsg(0,0,0);
    CkSendMsgBranch(updateComputesReturnEP,msg,0,updateComputesReturnChareID);
//  }
}

//
void
ComputeMgr::createCompute(ComputeID i, ComputeMap *map)
{
    Compute *c;
    PatchID pid2[2];
    PatchIDList pids;
    int trans2[2];
    SimParameters *simParams = Node::Object()->simParameters;

    PatchID pid8[8];
    int trans8[8];

    switch ( map->type(i) )
    {
    case computeNonbondedSelfType:
#ifdef NAMD_CUDA
        register_cuda_compute_self(i,map->computeData[i].pids[0].pid);
#else
        c = new ComputeNonbondedSelf(i,map->computeData[i].pids[0].pid,
                                     computeNonbondedWorkArrays,
                                     map->partition(i),map->partition(i)+1,
                                     map->numPartitions(i)); // unknown delete
        map->registerCompute(i,c);
        c->initialize();
#endif
        break;
    case computeLCPOType:
        for (int j = 0; j < 8; j++) {
          pid8[j] = map->computeData[i].pids[j].pid;
          trans8[j] = map->computeData[i].pids[j].trans;
        }
        c = new ComputeLCPO(i,pid8,trans8,
             computeNonbondedWorkArrays,
             map->partition(i),map->partition(i)+1,
             map->numPartitions(i), 8);
        map->registerCompute(i,c);
        c->initialize();
      
        break;
    case computeNonbondedPairType:
        pid2[0] = map->computeData[i].pids[0].pid;
        trans2[0] = map->computeData[i].pids[0].trans;
        pid2[1] = map->computeData[i].pids[1].pid;
        trans2[1] = map->computeData[i].pids[1].trans;
#ifdef NAMD_CUDA
        register_cuda_compute_pair(i,pid2,trans2);
#else
        c = new ComputeNonbondedPair(i,pid2,trans2,
                                     computeNonbondedWorkArrays,
                                     map->partition(i),map->partition(i)+1,
                                     map->numPartitions(i)); // unknown delete
        map->registerCompute(i,c);
        c->initialize();
#endif
        break;
#ifdef NAMD_CUDA
    case computeNonbondedCUDAType:
	c = computeNonbondedCUDAObject = new ComputeNonbondedCUDA(i,this); // unknown delete
	map->registerCompute(i,c);
	c->initialize();
	break;
#endif
    case computeExclsType:
        PatchMap::Object()->basePatchIDList(CkMyPe(),pids);
        c = new ComputeExcls(i,pids); // unknown delete
        map->registerCompute(i,c);
        c->initialize();
        break;
    case computeBondsType:
        PatchMap::Object()->basePatchIDList(CkMyPe(),pids);
        c = new ComputeBonds(i,pids); // unknown delete
        map->registerCompute(i,c);
        c->initialize();
        break;
    case computeAnglesType:
        PatchMap::Object()->basePatchIDList(CkMyPe(),pids);
        c = new ComputeAngles(i,pids); // unknown delete
        map->registerCompute(i,c);
        c->initialize();
        break;
    case computeDihedralsType:
        PatchMap::Object()->basePatchIDList(CkMyPe(),pids);
        c = new ComputeDihedrals(i,pids); // unknown delete
        map->registerCompute(i,c);
        c->initialize();
        break;
    case computeImpropersType:
        PatchMap::Object()->basePatchIDList(CkMyPe(),pids);
        c = new ComputeImpropers(i,pids); // unknown delete
        map->registerCompute(i,c);
        c->initialize();
        break;
    case computeTholeType:
        PatchMap::Object()->basePatchIDList(CkMyPe(),pids);
        c = new ComputeThole(i,pids); // unknown delete
        map->registerCompute(i,c);
        c->initialize();
        break;
    case computeAnisoType:
        PatchMap::Object()->basePatchIDList(CkMyPe(),pids);
        c = new ComputeAniso(i,pids); // unknown delete
        map->registerCompute(i,c);
        c->initialize();
        break;
    case computeCrosstermsType:
        PatchMap::Object()->basePatchIDList(CkMyPe(),pids);
        c = new ComputeCrossterms(i,pids); // unknown delete
        map->registerCompute(i,c);
        c->initialize();
        break;
    case computeSelfExclsType:
        c = new ComputeSelfExcls(i,map->computeData[i].pids[0].pid);
        map->registerCompute(i,c);
        c->initialize();
        break;
    case computeSelfBondsType:
        c = new ComputeSelfBonds(i,map->computeData[i].pids[0].pid);
        map->registerCompute(i,c);
        c->initialize();
        break;
    case computeSelfAnglesType:
        c = new ComputeSelfAngles(i,map->computeData[i].pids[0].pid);
        map->registerCompute(i,c);
        c->initialize();
        break;
    case computeSelfDihedralsType:
        c = new ComputeSelfDihedrals(i,map->computeData[i].pids[0].pid);
        map->registerCompute(i,c);
        c->initialize();
        break;
    case computeSelfImpropersType:
        c = new ComputeSelfImpropers(i,map->computeData[i].pids[0].pid);
        map->registerCompute(i,c);
        c->initialize();
        break;
    case computeSelfTholeType:
        c = new ComputeSelfThole(i,map->computeData[i].pids[0].pid);
        map->registerCompute(i,c);
        c->initialize();
        break;
    case computeSelfAnisoType:
        c = new ComputeSelfAniso(i,map->computeData[i].pids[0].pid);
        map->registerCompute(i,c);
        c->initialize();
        break;
    case computeSelfCrosstermsType:
        c = new ComputeSelfCrossterms(i,map->computeData[i].pids[0].pid);
        map->registerCompute(i,c);
        c->initialize();
        break;
#ifdef DPMTA
    case computeDPMTAType:
        c = new ComputeDPMTA(i); // unknown delete
        map->registerCompute(i,c);
        c->initialize();
        break;
#endif
#ifdef DPME
    case computeDPMEType:
        c = computeDPMEObject = new ComputeDPME(i,this); // unknown delete
        map->registerCompute(i,c);
        c->initialize();
        break;
#endif
    case optPmeType:
        c = new OptPmeCompute(i); // unknown delete
        map->registerCompute(i,c);
        c->initialize();
        break;
    case computePmeType:
        c = new ComputePme(i,map->computeData[i].pids[0].pid); // unknown delete
        map->registerCompute(i,c);
        c->initialize();
        break;
    case computeEwaldType:
        c = computeEwaldObject = new ComputeEwald(i,this); // unknown delete
        map->registerCompute(i,c);
        c->initialize();
        break;
    case computeFullDirectType:
        c = new ComputeFullDirect(i); // unknown delete
        map->registerCompute(i,c);
        c->initialize();
        break;
    case computeGlobalType:
        c = computeGlobalObject = new ComputeGlobal(i,this); // unknown delete
        map->registerCompute(i,c);
        c->initialize();
        break;
    case computeStirType:
        c = new ComputeStir(i,map->computeData[i].pids[0].pid); // unknown delete
        map->registerCompute(i,c);
        c->initialize();
        break;
    case computeExtType:
        c = new ComputeExt(i); // unknown delete
        map->registerCompute(i,c);
        c->initialize();
        break;
    case computeGBISserType: //gbis serial
        c = new ComputeGBISser(i);
        map->registerCompute(i,c);
        c->initialize();
        break;
    case computeMsmSerialType: // MSM serial
        c = new ComputeMsmSerial(i);
        map->registerCompute(i,c);
        c->initialize();
        break;
#ifdef CHARM_HAS_MSA
    case computeMsmMsaType: // MSM parallel long-range part using MSA
        c = new ComputeMsmMsa(i);
        map->registerCompute(i,c);
        c->initialize();
        break;
#endif
    case computeMsmType: // MSM parallel
        c = new ComputeMsm(i);
        map->registerCompute(i,c);
        c->initialize();
        break;
    case computeEFieldType:
        c = new ComputeEField(i,map->computeData[i].pids[0].pid); // unknown delete
        map->registerCompute(i,c);
        c->initialize();
        break;
        /* BEGIN gf */
    case computeGridForceType:
        c = new ComputeGridForce(i,map->computeData[i].pids[0].pid);
        map->registerCompute(i,c);
        c->initialize();
        break;
        /* END gf */
    case computeSphericalBCType:
        c = new ComputeSphericalBC(i,map->computeData[i].pids[0].pid); // unknown delete
        map->registerCompute(i,c);
        c->initialize();
        break;
    case computeCylindricalBCType:
        c = new ComputeCylindricalBC(i,map->computeData[i].pids[0].pid); // unknown delete
        map->registerCompute(i,c);
        c->initialize();
        break;
    case computeTclBCType:
        c = new ComputeTclBC(i); // unknown delete
        map->registerCompute(i,c);
        c->initialize();
        break;
    case computeRestraintsType:
        c = new ComputeRestraints(i,map->computeData[i].pids[0].pid); // unknown delete
        map->registerCompute(i,c);
        c->initialize();
        break;
    case computeConsForceType:
        c = new ComputeConsForce(i,map->computeData[i].pids[0].pid);
        map->registerCompute(i,c);
        c->initialize();
        break;
    case computeConsTorqueType:
        c = new ComputeConsTorque(i,map->computeData[i].pids[0].pid);
        map->registerCompute(i,c);
        c->initialize();
        break;
    default:
        NAMD_bug("Unknown compute type in ComputeMgr::createCompute().");
        break;
    }
}

void registerUserEventsForAllComputeObjs()
{
#ifdef TRACE_COMPUTE_OBJECTS
    ComputeMap *map = ComputeMap::Object();
    PatchMap *pmap = PatchMap::Object();     
    char user_des[50];
    int p1, p2;
    int adim, bdim, cdim;
    int t1, t2;
    int x1, y1, z1, x2, y2, z2;
    int dx, dy, dz;
    for (int i=0; i<map->numComputes(); i++)
    {
        memset(user_des, 0, 50);
        switch ( map->type(i) )
        {
        case computeNonbondedSelfType:
            sprintf(user_des, "computeNonBondedSelfType_%d_pid_%d", i, map->pid(i,0));
            break;
        case computeLCPOType:
            sprintf(user_des, "computeLCPOType_%d_pid_%d", i, map->pid(i,0));
            break;
        case computeNonbondedPairType:
            adim = pmap->gridsize_a();
            bdim = pmap->gridsize_b();
            cdim = pmap->gridsize_c();
            p1 = map->pid(i, 0);
            t1 = map->trans(i, 0);
            x1 = pmap->index_a(p1) + adim * Lattice::offset_a(t1);
            y1 = pmap->index_b(p1) + bdim * Lattice::offset_b(t1);
            z1 = pmap->index_c(p1) + cdim * Lattice::offset_c(t1);
            p2 = map->pid(i, 1);
            t2 = map->trans(i, 1);
            x2 = pmap->index_a(p2) + adim * Lattice::offset_a(t2);
            y2 = pmap->index_b(p2) + bdim * Lattice::offset_b(t2);
            z2 = pmap->index_c(p2) + cdim * Lattice::offset_c(t2);
            dx = abs(x1-x2);
            dy = abs(y1-y2);
            dz = abs(z1-z2);
            sprintf(user_des, "computeNonBondedPairType_%d(%d,%d,%d)", i, dx,dy,dz);
            break;
        case computeExclsType:
            sprintf(user_des, "computeExclsType_%d", i);
            break;
        case computeBondsType:
            sprintf(user_des, "computeBondsType_%d", i);
            break;
        case computeAnglesType:
            sprintf(user_des, "computeAnglesType_%d", i);
            break;
        case computeDihedralsType:
            sprintf(user_des, "computeDihedralsType_%d", i);
            break;
        case computeImpropersType:
            sprintf(user_des, "computeImpropersType_%d", i);
            break;
        case computeTholeType:
            sprintf(user_des, "computeTholeType_%d", i);
            break;
        case computeAnisoType:
            sprintf(user_des, "computeAnisoType_%d", i);
            break;
        case computeCrosstermsType:
            sprintf(user_des, "computeCrosstermsType_%d", i);
            break;
        case computeSelfExclsType:
            sprintf(user_des, "computeSelfExclsType_%d", i);
            break;
        case computeSelfBondsType:
            sprintf(user_des, "computeSelfBondsType_%d", i);
            break;
        case computeSelfAnglesType:
            sprintf(user_des, "computeSelfAnglesType_%d", i);
            break;
        case computeSelfDihedralsType:
            sprintf(user_des, "computeSelfDihedralsType_%d", i);
            break;
        case computeSelfImpropersType:
            sprintf(user_des, "computeSelfImpropersType_%d", i);
            break;
        case computeSelfTholeType:
            sprintf(user_des, "computeSelfTholeType_%d", i);
            break;
        case computeSelfAnisoType:
            sprintf(user_des, "computeSelfAnisoType_%d", i);
            break;
        case computeSelfCrosstermsType:
            sprintf(user_des, "computeSelfCrosstermsType_%d", i);
            break;
#ifdef DPMTA
        case computeDPMTAType:
            sprintf(user_des, "computeDPMTAType_%d", i);
            break;
#endif
#ifdef DPME
        case computeDPMEType:
            sprintf(user_des, "computeDPMEType_%d", i);
            break;
#endif
        case computePmeType:
            sprintf(user_des, "computePMEType_%d", i);
            break;
        case computeEwaldType:
            sprintf(user_des, "computeEwaldType_%d", i);
            break;
        case computeFullDirectType:
            sprintf(user_des, "computeFullDirectType_%d", i);
            break;
        case computeGlobalType:
            sprintf(user_des, "computeGlobalType_%d", i);
            break;
        case computeStirType:
            sprintf(user_des, "computeStirType_%d", i);
            break;
        case computeExtType:
            sprintf(user_des, "computeExtType_%d", i);
            break;
        case computeEFieldType:
            sprintf(user_des, "computeEFieldType_%d", i);
            break;
            /* BEGIN gf */
        case computeGridForceType:
            sprintf(user_des, "computeGridForceType_%d", i);
            break;
            /* END gf */
        case computeSphericalBCType:
            sprintf(user_des, "computeSphericalBCType_%d", i);
            break;
        case computeCylindricalBCType:
            sprintf(user_des, "computeCylindricalBCType_%d", i);
            break;
        case computeTclBCType:
            sprintf(user_des, "computeTclBCType_%d", i);
            break;
        case computeRestraintsType:
            sprintf(user_des, "computeRestraintsType_%d", i);
            break;
        case computeConsForceType:
            sprintf(user_des, "computeConsForceType_%d", i);
            break;
        case computeConsTorqueType:
            sprintf(user_des, "computeConsTorqueType_%d", i);
            break;
        default:
            NAMD_bug("Unknown compute type in ComputeMgr::registerUserEventForAllComputeObjs().");
            break;
        }
	int user_des_len = strlen(user_des);
	char *user_des_cst = new char[user_des_len+1];
	memcpy(user_des_cst, user_des, user_des_len);
	user_des_cst[user_des_len] = 0;
	//Since the argument in traceRegisterUserEvent is supposed
	//to be a const string which will not be copied inside the
	//function when a new user event is created, user_des_cst 
	//has to be allocated in heap.
        int reEvenId = traceRegisterUserEvent(user_des_cst, TRACE_COMPOBJ_IDOFFSET+i);
	//printf("Register user event (%s) with id (%d)\n", user_des, reEvenId);
    }
#else
    return;
#endif
}

void
ComputeMgr::createComputes(ComputeMap *map)
{
    Node *node = Node::Object();
    SimParameters *simParams = node->simParameters;
    int myNode = node->myid();

    ComputeNonbondedUtil::select();

    if ( simParams->globalForcesOn && !myNode )
    {
        DebugM(4,"Mgr running on Node "<<CkMyPe()<<"\n");
        /* create a master server to allow multiple masters */
        masterServerObject = new GlobalMasterServer(this,
                PatchMap::Object()->numNodesWithPatches());

        /* create the individual global masters */
        // masterServerObject->addClient(new GlobalMasterTest());
        if (simParams->tclForcesOn)
            masterServerObject->addClient(new GlobalMasterTcl());
        if (simParams->IMDon && ! simParams->IMDignore)
            masterServerObject->addClient(new GlobalMasterIMD());

        if (simParams->SMDOn)
            masterServerObject->addClient(
                new GlobalMasterSMD(simParams->SMDk, simParams->SMDk2,
				    simParams->SMDVel,
                                    simParams->SMDDir, simParams->SMDOutputFreq,
                                    simParams->firstTimestep, simParams->SMDFile,
                                    node->molecule->numAtoms)
            );
            
        if (simParams->symmetryOn && 
          (simParams->firstTimestep < simParams->symmetryLastStep || 
          simParams->symmetryLastStep == -1))
            masterServerObject->addClient(new GlobalMasterSymmetry());    
        if (simParams->TMDOn)
            masterServerObject->addClient(new GlobalMasterTMD());
        if (simParams->miscForcesOn)
            masterServerObject->addClient(new GlobalMasterMisc());
        if ( simParams->freeEnergyOn )
            masterServerObject->addClient(new GlobalMasterFreeEnergy());
		if ( simParams->colvarsOn )
			masterServerObject->addClient(new GlobalMasterColvars());

    }

    if ( !myNode && simParams->IMDon && simParams->IMDignore ) {
      // GlobalMasterIMD constructor saves pointer to node->IMDOutput object
      new GlobalMasterIMD();
    }

#ifdef NAMD_CUDA
    bool deviceIsMine = ( cuda_device_pe() == CkMyPe() );
#endif

    for (int i=0; i < map->nComputes; i++)
    {
        if ( ! ( i % 100 ) )
        {
        }
#ifdef NAMD_CUDA
        switch ( map->type(i) )
        {
          case computeNonbondedSelfType:
          case computeNonbondedPairType:
            if ( ! deviceIsMine ) continue;
            if ( ! cuda_device_shared_with_pe(map->computeData[i].node) ) continue;
          break;
          case computeNonbondedCUDAType:
            if ( ! deviceIsMine ) continue;
          default:
            if ( map->computeData[i].node != myNode ) continue;
        }
#else
        if ( map->computeData[i].node != myNode ) continue;
#endif
        DebugM(1,"Compute " << i << '\n');
        DebugM(1,"  node = " << map->computeData[i].node << '\n');
        DebugM(1,"  type = " << map->computeData[i].type << '\n');
        DebugM(1,"  numPids = " << map->computeData[i].numPids << '\n');
        DebugM(1,"  numPidsAllocated = " << map->computeData[i].numPidsAllocated << '\n');
        for (int j=0; j < map->computeData[i].numPids; j++)
        {
            //      DebugM(1,"  pid " << map->computeData[i].pids[j] << '\n');
            if (!((j+1) % 6))
                DebugM(1,'\n');
        }
        DebugM(1,"\n---------------------------------------");
        DebugM(1,"---------------------------------------\n");

        createCompute(i, map);

    }

#ifdef NAMD_CUDA
    if ( computeNonbondedCUDAObject ) {
      computeNonbondedCUDAObject->assignPatches();
    }
#endif

}

#if 0
void ComputeMgr:: sendComputeGlobalConfig(ComputeGlobalConfigMsg *msg)
{
    (CProxy_ComputeMgr(CkpvAccess(BOCclass_group).computeMgr)).recvComputeGlobalConfig(msg);
}

void ComputeMgr:: recvComputeGlobalConfig(ComputeGlobalConfigMsg *msg)
{
    if ( computeGlobalObject )
    {
        computeGlobalObject->recvConfig(msg);
    }
    else if ( ! (PatchMap::Object())->numHomePatches() ) delete msg;
    else NAMD_die("ComputeMgr::computeGlobalObject is NULL!");
}
#endif

void ComputeMgr:: sendComputeGlobalData(ComputeGlobalDataMsg *msg)
{
    CProxy_ComputeMgr cm(CkpvAccess(BOCclass_group).computeMgr);
    cm[0].recvComputeGlobalData(msg);
}

void ComputeMgr:: recvComputeGlobalData(ComputeGlobalDataMsg *msg)
{
    if (masterServerObject)  // make sure it has been initialized
    {
        masterServerObject->recvData(msg);
    }
    else NAMD_die("ComputeMgr::masterServerObject is NULL!");
}

void ComputeMgr:: sendComputeGlobalResults(ComputeGlobalResultsMsg *msg)
{
    (CProxy_ComputeMgr(CkpvAccess(BOCclass_group).computeMgr)).recvComputeGlobalResults(msg);
}

void ComputeMgr:: recvComputeGlobalResults(ComputeGlobalResultsMsg *msg)
{
    if ( computeGlobalObject )
    {
        CmiEnableUrgentSend(1);
        computeGlobalObject->recvResults(msg);
        CmiEnableUrgentSend(0);
    }
    else if ( ! (PatchMap::Object())->numHomePatches() ) delete msg;
    else NAMD_die("ComputeMgr::computeGlobalObject is NULL!");
}

/*
 * Begin Ewald messages
 */
void ComputeMgr:: sendComputeEwaldData(ComputeEwaldMsg *msg)
{
    if (computeEwaldObject)
    {
        int node = computeEwaldObject->getMasterNode();
        CProxy_ComputeMgr cm(CkpvAccess(BOCclass_group).computeMgr);
        cm[node].recvComputeEwaldData(msg);
    }
    else if (!PatchMap::Object()->numHomePatches())
    {
        CkPrintf("skipping message on Pe(%d)\n", CkMyPe());
        delete msg;
    }
    else NAMD_die("ComputeMgr::computeEwaldObject is NULL!");
}

void ComputeMgr:: recvComputeEwaldData(ComputeEwaldMsg *msg)
{
    if (computeEwaldObject)
        computeEwaldObject->recvData(msg);
    else NAMD_die("ComputeMgr::computeEwaldObject in recvData is NULL!");
}

void ComputeMgr:: sendComputeEwaldResults(ComputeEwaldMsg *msg)
{
    (CProxy_ComputeMgr(CkpvAccess(BOCclass_group).computeMgr)).recvComputeEwaldResults(msg);
}

void ComputeMgr::recvComputeEwaldResults(ComputeEwaldMsg *msg)
{
    if (computeEwaldObject) {
        CmiEnableUrgentSend(1);
        computeEwaldObject->recvResults(msg);
        CmiEnableUrgentSend(0);
    }
    else if ( ! (PatchMap::Object())->numHomePatches() ) delete msg;
    else NAMD_die("ComputeMgr::computeEwaldObject in recvResults is NULL!");
}

void ComputeMgr:: sendComputeDPMEData(ComputeDPMEDataMsg *msg)
{
    if ( computeDPMEObject )
    {
#ifdef DPME
        int node = computeDPMEObject->getMasterNode();
        CProxy_ComputeMgr cm(CkpvAccess(BOCclass_group).computeMgr);
        cm.recvComputeDPMEData(msg,node);
#endif
    }
    else if ( ! (PatchMap::Object())->numHomePatches() ) delete msg;
    else NAMD_die("ComputeMgr::computeDPMEObject is NULL!");
}

void ComputeMgr:: recvComputeDPMEData(ComputeDPMEDataMsg *msg)
{
    if ( computeDPMEObject )
    {
#ifdef DPME
        computeDPMEObject->recvData(msg);
#endif
    }
    else if ( ! (PatchMap::Object())->numHomePatches() ) delete msg;
    else NAMD_die("ComputeMgr::computeDPMEObject is NULL!");
}

void ComputeMgr:: sendComputeDPMEResults(ComputeDPMEResultsMsg *msg, int node)
{
    CProxy_ComputeMgr cm(CkpvAccess(BOCclass_group).computeMgr);
    cm[node].recvComputeDPMEResults(msg);
}

void ComputeMgr:: recvComputeDPMEResults(ComputeDPMEResultsMsg *msg)
{
    if ( computeDPMEObject )
    {
#ifdef DPME
        computeDPMEObject->recvResults(msg);
#endif
    }
    else if ( ! (PatchMap::Object())->numHomePatches() ) delete msg;
    else NAMD_die("ComputeMgr::computeDPMEObject is NULL!");
}

void ComputeMgr::recvComputeConsForceMsg(ComputeConsForceMsg *msg)
{
    Molecule *m = Node::Object()->molecule;
    delete [] m->consForceIndexes;
    delete [] m->consForce;
    int n = msg->aid.size();
    if (n > 0)
    {
        m->consForceIndexes = new int32[m->numAtoms];
        m->consForce = new Vector[n];
        int i;
        for (i=0; i<m->numAtoms; i++) m->consForceIndexes[i] = -1;
        for (i=0; i<msg->aid.size(); i++)
        {
            m->consForceIndexes[msg->aid[i]] = i;
            m->consForce[i] = msg->f[i];
        }
    }
    else
    {
        m->consForceIndexes = NULL;
        m->consForce = NULL;
    }
    delete msg;
}

void ComputeMgr::sendYieldDevice(int pe) {
    CProxy_ComputeMgr cm(CkpvAccess(BOCclass_group).computeMgr);
    cm[pe].recvYieldDevice(CkMyPe());
}

void ComputeMgr::recvYieldDevice(int pe) {
#ifdef NAMD_CUDA
    computeNonbondedCUDAObject->recvYieldDevice(pe);
#endif
}

void ComputeMgr::sendBuildCudaForceTable() {
    CProxy_ComputeMgr cm(CkpvAccess(BOCclass_group).computeMgr);
    int pe = CkNodeFirst(CkMyNode());
    int end = pe + CkNodeSize(CkMyNode());
    for( ; pe != end; ++pe ) {
      cm[pe].recvBuildCudaForceTable();
    }
}

#ifdef NAMD_CUDA
  void build_cuda_force_table();
#endif

void ComputeMgr::recvBuildCudaForceTable() {
#ifdef NAMD_CUDA
    build_cuda_force_table();
#endif
}

class NonbondedCUDASlaveMsg : public CMessage_NonbondedCUDASlaveMsg {
public:
  int index;
  ComputeNonbondedCUDA *master;
};

void ComputeMgr::sendCreateNonbondedCUDASlave(int pe, int index) {
  NonbondedCUDASlaveMsg *msg = new NonbondedCUDASlaveMsg;
  msg->master = computeNonbondedCUDAObject;
  msg->index = index;
  CProxy_ComputeMgr cm(CkpvAccess(BOCclass_group).computeMgr);
  cm[pe].recvCreateNonbondedCUDASlave(msg);
}

void ComputeMgr::recvCreateNonbondedCUDASlave(NonbondedCUDASlaveMsg *msg) {
#ifdef NAMD_CUDA
  new ComputeNonbondedCUDA(msg->master->cid,this,msg->master,msg->index);
#endif
}

void ComputeMgr::sendNonbondedCUDASlaveReady(int pe, int np, int ac, int seq) {
  CProxy_ComputeMgr cm(CkpvAccess(BOCclass_group).computeMgr);
  cm[pe].recvNonbondedCUDASlaveReady(np,ac,seq);
}

void ComputeMgr::recvNonbondedCUDASlaveReady(int np, int ac, int seq) {
  for ( int i=0; i<np; ++i ) {
    computeNonbondedCUDAObject->patchReady(-1,ac,seq);
  }
}

void ComputeMgr::sendNonbondedCUDASlaveEnqueue(ComputeNonbondedCUDA *c, int pe, int seq, int prio, int ws) {
  if ( ws == 2 && c->localHostedPatches.size() == 0 ) return;
  LocalWorkMsg *msg = ( ws == 1 ? c->localWorkMsg : c->localWorkMsg2 );
  msg->compute = c;
  int type = c->type();
  int cid = c->cid;
  SET_PRIORITY(msg,seq,prio);
  CProxy_WorkDistrib wdProxy(CkpvAccess(BOCclass_group).workDistrib);
  wdProxy[pe].enqueueCUDA(msg);
}

#include "ComputeMgr.def.h"

