#include "UserCode/OmtfAnalysis/interface/SynchroCountsGrabber.h"

#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESTransientHandle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESWatcher.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "CondFormats/RPCObjects/interface/RPCEMap.h"
#include "CondFormats/DataRecord/interface/RPCEMapRcd.h"
#include "CondFormats/RPCObjects/interface/RPCReadOutMapping.h"
#include "CondFormats/RPCObjects/interface/LinkBoardSpec.h"
#include "DataFormats/MuonDetId/interface/RPCDetId.h"

#include "DataFormats/MuonReco/interface/Muon.h"

#include "Geometry/RPCGeometry/interface/RPCGeometry.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"

#include "DataFormats/RPCDigi/interface/RPCRawSynchro.h"
#include "TrackingTools/PatternTools/interface/Trajectory.h"
#include "DataFormats/RPCRecHit/interface/RPCRecHitCollection.h"


namespace {
  edm::EDGetTokenT<RPCRawSynchro::ProdItem> theSynchroCounts_Tag;
  edm::EDGetTokenT<std::vector<Trajectory> > theRefittedGlbMu_Tag;
  edm::EDGetTokenT<RPCRecHitCollection> theRPCRecHitCollection_Tag;
  edm::ESGetToken<RPCEMap, RPCEMapRcd> theReadoutMappingToken;
  edm::ESGetToken<RPCGeometry, MuonGeometryRecord> theRPCGeometryToken;
}

SynchroCountsGrabber::SynchroCountsGrabber(const edm::ParameterSet& cfg, edm::ConsumesCollector cColl)
 : trackAtSurface(cColl), theCabling(0), 
   theSelector(cfg.getParameter<edm::ParameterSet>("synchroSelector"),cColl), 
   deltaR_MuonToDetUnit_cutoff(cfg.getParameter<double>("deltaR_MuonToDetUnit_cutoff")), 
   checkInside(cfg.getParameter<bool>("checkInside")),
   theNoSynchroWarning(false)
{ 
  theSynchroCounts_Tag = cColl.consumes<RPCRawSynchro::ProdItem> (cfg.getParameter<edm::InputTag>("rawSynchroTag") );
//  theRefittedGlbMu_Tag = cColl.consumes<std::vector<Trajectory> >(edm::InputTag("refittedMuons","Refitted"));
  theRPCRecHitCollection_Tag = cColl.consumes<RPCRecHitCollection>(edm::InputTag("rpcRecHits"));
  theReadoutMappingToken = cColl.esConsumes<RPCEMap, RPCEMapRcd>(); 
  theRPCGeometryToken = cColl.esConsumes<RPCGeometry, MuonGeometryRecord>();
}

SynchroCountsGrabber::~SynchroCountsGrabber()
{
  if (theNoSynchroWarning) std::cout <<" **** SynchroCountsGrabber      **** WARNING - NoSynchroWarning was set!" << std::endl;
  delete theCabling;
}

RPCRawSynchro::ProdItem SynchroCountsGrabber::counts(const edm::Event &ev, const edm::EventSetup &es)
{
  RPCRawSynchro::ProdItem result;

  if (theMapWatcher.check(es)) {
    delete theCabling;
    edm::ESTransientHandle<RPCEMap> readoutMapping =  es.getTransientHandle(theReadoutMappingToken);
    theCabling = readoutMapping->convert();
  //LogTrace("") << "SynchroCountsGrabber - record has CHANGED!!, read map, VERSION: " << theCabling->version();
    std::cout    << "SynchroCountsGrabber - record has CHANGED!!, read map, VERSION: " << theCabling->version() << std::endl;
  }

  edm::Handle<RPCRawSynchro::ProdItem> synchroCounts;
  ev.getByLabel("muonRPCDigis",synchroCounts);
//  std::cout <<"synchroCounts SIZE: " << synchroCounts->size() << std::endl;
  
  if (!synchroCounts.isValid()) {
    theNoSynchroWarning = true;
    return result;
  }

  const RPCGeometry & rpcGeometry = es.getData(theRPCGeometryToken);

//  trackAtSurface.prepare(theMuon, ev,es);

  for(RPCRawSynchro::ProdItem::const_iterator it = synchroCounts->begin(); it != synchroCounts->end(); ++it) {
    const LinkBoardElectronicIndex & path = it->first;
    const  std::vector<FebConnectorSpec> & febs = theCabling->location(path)->febs();
    std::map<uint32_t,bool> dets;
    for (std::vector<FebConnectorSpec>::const_iterator iif = febs.begin(); iif != febs.end(); ++iif) dets[iif->rawId()] = true;
    for ( std::map<uint32_t,bool>::const_iterator im = dets.begin(); im != dets.end(); ++im) {
      RPCDetId rpcDet(im->first);
      const GeomDet *geomDet = rpcGeometry.idToDet(rpcDet);
      GlobalPoint detPosition = geomDet->position();
      if (deltaR(theMuon->eta(), theMuon->phi(), detPosition.eta(), detPosition.phi()) > deltaR_MuonToDetUnit_cutoff) continue;
      TrajectoryStateOnSurface stateAtDet = trackAtSurface.atDetFromClose(theMuon, rpcDet,detPosition,es);
      if (!stateAtDet.isValid()) continue;
      if (checkInside && !(geomDet->surface().bounds().inside(stateAtDet.localPosition()))) continue;
      if (!theSelector.checkTraj(stateAtDet, rpcDet, ev, es)) continue;
      result.push_back(*it);
    }
  }
//  std::cout <<"synchroCounts SIZE: " << synchroCounts->size() <<" Result size: "<< result.size() << std::endl;
  return result;
}

RPCRawSynchro::ProdItem SynchroCountsGrabber::counts(const edm::Event &ev, const edm::EventSetup &es, float eta, float phi)
{
  RPCRawSynchro::ProdItem result;

  if (theMapWatcher.check(es)) {
    delete theCabling;
//    edm::ESTransientHandle<RPCEMap> readoutMapping;
//    es.get<RPCEMapRcd>().get(readoutMapping);
    edm::ESTransientHandle<RPCEMap> readoutMapping =  es.getTransientHandle(theReadoutMappingToken);
    theCabling = readoutMapping->convert();
    //LogTrace("") << "SynchroCountsGrabber - record has CHANGED!!, read map, VERSION: " << theCabling->version();
    std::cout << "SynchroCountsGrabber - record has CHANGED!!, read map, VERSION: " << theCabling->version() << std::endl;
  }

  edm::Handle<RPCRawSynchro::ProdItem> synchroCounts;
  ev.getByLabel("muonRPCDigis",synchroCounts);
  std::cout <<"synchroCounts SIZE (all): " << synchroCounts->size() << std::endl;
  
  if (!synchroCounts.isValid()) {
    theNoSynchroWarning = true;
    return result;
  }

//  edm::ESHandle<RPCGeometry> rpcGeometry;
//  es.get<MuonGeometryRecord>().get(rpcGeometry);
  const RPCGeometry & rpcGeometry = es.getData(theRPCGeometryToken);


  for(RPCRawSynchro::ProdItem::const_iterator it = synchroCounts->begin(); it != synchroCounts->end(); ++it) {
    bool takeIt = false;
    const LinkBoardElectronicIndex & path = it->first;
    const  std::vector<FebConnectorSpec> & febs = theCabling->location(path)->febs();
    std::map<uint32_t,bool> dets;
    for (std::vector<FebConnectorSpec>::const_iterator iif = febs.begin(); iif != febs.end(); ++iif) dets[iif->rawId()] = true;
    for ( std::map<uint32_t,bool>::const_iterator im = dets.begin(); im != dets.end(); ++im) {
      RPCDetId rpcDet(im->first);
      const GeomDet *geomDet = rpcGeometry.idToDet(rpcDet);
      GlobalPoint detPosition = geomDet->position();
      if (deltaR(eta,phi, detPosition.eta(), detPosition.phi()) > deltaR_MuonToDetUnit_cutoff) continue;
      takeIt = true;
    }
    if (takeIt) result.push_back(*it);
  }
  std::cout <<"synchroCounts SIZE (sel): " << result.size() << std::endl;
  return result;
}

