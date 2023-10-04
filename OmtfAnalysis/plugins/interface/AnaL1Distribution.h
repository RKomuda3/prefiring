#ifndef UserCode_OmtfAnalysis_AnaL1Distribution_H
#define UserCode_OmtfAnalysis_AnaL1Distribution_H

class TObjArray;
class TH2D;
class MuonObj;
class MuonObjColl;
class L1ObjColl;
class L1Obj;
class EventObj;
class TrackObj;

#include <vector>
#include <map>
#include <string>
#include "UserCode/OmtfAnalysis/interface/Utilities.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

namespace edm { class ParameterSet;}

class AnaL1Distribution {
public:
  AnaL1Distribution(const edm::ParameterSet & cfg);
  void init(TObjArray& histos);
  void run( const EventObj* ev, const MuonObjColl * muons, const TrackObj *track, const L1ObjColl * l1Objs);
  bool debug;

private:
   edm::ParameterSet theCfg;
   std::map< std::string, TH2D* > theHistoMap;
};

#endif
