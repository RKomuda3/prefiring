import FWCore.ParameterSet.Config as cms
import os
import sys
#import commands

process = cms.Process("Analysis")
process.source = cms.Source("EmptySource")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1))

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.Geometry.GeometryExtended2017Reco_cff')
#process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
#from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run3_data', '')
#process.GlobalTag.globaltag = '113X_dataRun3_Express_v2'


process.omtfAnalysis = cms.EDAnalyzer("OmtfTreeAnalysis",
  histoFileName = cms.string("omtfAnalysis.root"),
  treeFileNames = cms.vstring(
#        "omtfTree.root"
#  "../jobs/Run2023B_JetMET/crab_JetMET0_jobVer3/omtfTree.root",
#  "../jobs/Run2023B_JetMET/crab_JetMET1_jobVer3/omtfTree.root",
#  "../jobs/Run2023C_JetMET/crab_JetMET0_jobVer3/omtfTree.root",
#  "../jobs/Run2023C_JetMET/crab_JetMET1_jobVer3/omtfTree.root",
   "/afs/cern.ch/work/k/konec/CMSSW_13_1_0.ana/src/UserCode/OmtfAnalysis/jobs/Run2023B_Muon/crab_Muon0_jobVer3/omtfTree.root",
   "/afs/cern.ch/work/k/konec/CMSSW_13_1_0.ana/src/UserCode/OmtfAnalysis/jobs/Run2023B_Muon/crab_Muon1_jobVer3/omtfTree.root",
   "/afs/cern.ch/work/k/konec/CMSSW_13_1_0.ana/src/UserCode/OmtfAnalysis/jobs/Run2023C_Muon/crab_Muon0_jobVer3/omtfTree.root",
   "/afs/cern.ch/work/k/konec/CMSSW_13_1_0.ana/src/UserCode/OmtfAnalysis/jobs/Run2023C_Muon/crab_Muon1_jobVer3/omtfTree.root",
  ),
  filterByAnaEvent = cms.bool(True),
  anaEvent = cms.PSet(
   skipRuns = cms.vuint32(),
#   onlyRuns = cms.vuint32(299149),
  ),

#   anaSecMuSel = cms.PSet(
#     triggMuon = cms.PSet (
#       requireCharge      = cms.int32(1),
###     requireEtaSign     = cms.int32(1),
###      requireOutsideOMTF  = cms.bool(True),
#       requireTight        = cms.bool(True),
#       requireUnique       = cms.bool(True),
#       requireHLT          = cms.bool(True),
#       requireIsoForHLTIso = cms.bool(True),
#       minAcceptMuPtVal    = cms.double(24.),
#       minMatchStations    = cms.uint32(2),
#       maxMuEtaVal         = cms.double(2.4),
#       minAcceptL1PtVal    = cms.double(22.0),
#       maxL1DeltaR         = cms.double(0.1),
#       minAccepL1Q         = cms.int32(12), 
#     ),
#    probeMuon = cms.PSet (
#       requireCharge    = cms.int32(-1),
###     requireEtaSign   = cms.int32(-1),
###     requireInsideOMTF = cms.bool(True),
#       requireUnique     = cms.bool(True),
#       requireLoose      = cms.bool(True),
#       requireMedium     = cms.bool(True),
#       requireTight      = cms.bool(False),
#       requireTkIso      = cms.bool(False),
#       requirePfIso      = cms.bool(False),
#       minTrgMuDeltaR    = cms.double(0.6),
#       maxMuEtaVal       = cms.double(2.4),
#     ),
#   ), 

  filterByAnaMuonDistribution = cms.bool(False),
  anaMuonDistribution = cms.PSet (
    requireUnique  = cms.bool(True),
    requireGlobal  = cms.bool(True),
    requireInner   = cms.bool(True),
    requireOuter   = cms.bool(False),
    requireLoose   = cms.bool(True),
    requireMedium  = cms.bool(True),
    requireTight   = cms.bool(False),
    requireTkIso   = cms.bool(False),
    requirePFIso   = cms.bool(False),
    chi2Norm  = cms.double(2.),
    ptMin     = cms.double(3.0),
    etaMax    = cms.double(2.4),
    minNumberOfMatchedStations = cms.uint32(2),
    minNumberTkHits = cms.uint32(6),
    minNumberRpcHits = cms.uint32(0),
    minNumberDtCscHits = cms.uint32(0),
    minNumberRpcDtCscHits = cms.uint32(2),
  ),

  filterByAnaMenu = cms.bool(True),
  anaMenu = cms.PSet( #OR of conditions for L1 and separately for HLT
#   acceptL1_Names  = cms.vstring("any"),
#   acceptL1_Names  = cms.vstring("L1_IsolatedBunch","L1_FirstBunchInTrain","L1_FirstCollisionInOrbit","L1_FirstCollisionInTrain"),
    acceptL1_Names  = cms.vstring("L1_FirstBunchInTrain"),
    acceptHLT_Names = cms.vstring("any"),
#   acceptHLT_Names = cms.vstring("HLT_AK8PF","HLT_Jet","HLT_HT","HLT_PFJet","HLT_PFHT"),
  ),

   anaDataEmul =  cms.PSet( bxMin=cms.int32(-3), bxMax=cms.int32(4)),
   anaEff =  cms.PSet(),     
#
# NEXT: to do (?) propagete probe to AnaTime and check prefire for that muon only
# NOW: do not requre SecMuSel and filterBySecMu because triggering muon is independently found in AnaTime 
# FIXME: possible add another method for NEXT 
   anaTime = cms.PSet( 
     maxPtForDistributions = cms.double(21.99), minPtForDistributions = cms.double(0.), 
     requireOtherMuTrg = cms.bool(True), requireOtherMuTrgL1Pt = cms.double(22), 
     deltaRMatching = cms.double(0.3)
   ),
   anaSynch = cms.PSet(),
   anaDiMu = cms.PSet(),
   anaL1Distribution = cms.PSet(),
)

process.p = cms.Path(process.omtfAnalysis)
