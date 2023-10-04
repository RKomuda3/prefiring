#ifndef SynchroSelectorMuon_H
#define SynchroSelectorMuon_H

#include "UserCode/OmtfAnalysis/interface/SynchroSelector.h"

class TObjArray;
class TH1F;
class TrajectoryStateOnSurface;

#include "FWCore/Framework/interface/ConsumesCollector.h"


class SynchroSelectorMuon : public SynchroSelector {
public:
   SynchroSelectorMuon(){}
   SynchroSelectorMuon(const edm::ParameterSet&, TObjArray&,  edm::ConsumesCollector cColl );
   SynchroSelectorMuon(const edm::ParameterSet&, edm::ConsumesCollector cColl );
   void initHistos(TObjArray& histos);
   virtual ~SynchroSelectorMuon(){}
   virtual bool takeIt(const RPCDetId & det, const edm::Event&ev, const edm::EventSetup& es);
   bool checkTraj( TrajectoryStateOnSurface & aTSOS, const RPCDetId & det, 
                   const edm::Event&ev, const edm::EventSetup& es);
private:
   TH1F * hDxy ;
};
#endif 
