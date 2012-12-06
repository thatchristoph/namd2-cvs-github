/* 
 * This class implements a dummy load balancer. One application is to provide a distributed
 * load balancer by means of the hybrid approach, where the first level doesn't balance
 * the load at all.
 */ 

#ifndef _NAMDDUMMYLB_H_
#define _NAMDDUMMYLB_H_

#include <CentralLB.h>
#include "NamdDummyLB.decl.h"

#include "Node.h"
#include "PatchMap.h"
#include "SimParameters.h"
#include "RefineOnly.h"
#include "Alg7.h"
#include "AlgRecBisection.h"
#include "InfoStream.h"
#include "TorusLB.h"
#include "RefineTorusLB.h"

class NamdDummyLB : public CentralLB {

public:
  NamdDummyLB();
  NamdDummyLB(CkMigrateMessage *);
#if CHARM_VERSION > 60301
  void work(LDStats* stats);
#else
  void work(LDStats* stats, int n_pes);
#endif

private:
  CmiBool QueryBalanceNow(int step);
  CmiBool QueryDumpData();

};

void CreateNamdDummyLB();
NamdDummyLB *AllocateNamdDummyLB();

#endif /* _NAMDDUMMYLB_H_ */
