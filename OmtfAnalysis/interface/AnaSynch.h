#ifndef UserCode_L1RpcTriggerAnalysis_AnaSynch_H
#define UserCode_L1RpcTriggerAnalysis_AnaSynch_H

class TObjArray;
class TH1D;
class MuonObj;
class EventObj;

#include <vector>
#include <map>
#include <string>
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "CondFormats/RPCObjects/interface/LinkBoardElectronicIndex.h"
#include "DataFormats/RPCDigi/interface/RPCRawSynchro.h"
#include "UserCode/OmtfAnalysis/interface/RPCLinkSynchroStat.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"

namespace edm { class EventSetup; class Run; class ParameterSet; }


class AnaSynch {
public:
  AnaSynch(const edm::ParameterSet& cfg,  edm::ConsumesCollector cColl);
  void init(TObjArray& histos);
  void run( const EventObj* event, const MuonObj* muon, const  RPCRawSynchro::ProdItem & synchro);
  void beginRun(const edm::Run& ru, const edm::EventSetup& es);
  void endJob();
private:
  RPCLinkSynchroStat theSynchroStat;
};

#endif

