import FWCore.ParameterSet.Config as cms
import commands
import os

from Configuration.StandardSequences.Eras import eras
process = cms.Process('OmtfTree',eras.Run2_2016)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

#
# For processing single files insert lines with 'file:/PATH/FILE.root'
# alernatively one can use 'root://xrootd.unl.edu//store/express/Run2015A.....root'
# or                       '/store/express/Run2015A/ExpressPhysics/FEVT/...root'
# (there is 255 file limit though). Can be empty for crab.
#
process.source = cms.Source("PoolSource", 
#fileNames = cms.untracked.vstring('file:/afs/cern.ch/work/k/konec/data/runs/run274002-Cosmic_12D87234-1823-E611-9A10-02163E014572.root'),
fileNames = cms.untracked.vstring( 
                                   'file:/afs/cern.ch/work/k/konec/data/runs/run274094_32779D03-C323-E611-BE75-02163E011821.root',
                                   '/store/express/Run2016B/ExpressPhysics/FEVT/Express-v2/000/274/094/00000/0070F70B-BD23-E611-A5C3-02163E011A28.root',
                                   '/store/express/Run2016B/ExpressPhysics/FEVT/Express-v2/000/274/094/00000/06E08AB1-C723-E611-AD43-02163E011C18.root',
                                   'file:/afs/cern.ch/work/k/konec/data/runs/run274157-384_12A67A1E-E324-E611-831F-02163E012658.root',
                                  ),
#skipEvents =  cms.untracked.uint32(220)
)

#
# import of standard configurations
#
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.Geometry.GeometryExtended2016Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.load('Configuration.StandardSequences.RawToDigi_Data_cff')
process.load('EventFilter.L1TRawToDigi.emtfStage2Digis_cfi')

#
# set proper GlobalTag
#
#process.GlobalTag.globaltag = 'FT_R_53_V10::All' #rereco 2012ABC
#process.GlobalTag.globaltag = 'FT_R_53_V18::All' #rereco 2012ABC
#process.GlobalTag.globaltag = 'FT_R_53_V21::All' #rereco 2012D
#process.GlobalTag.globaltag  = 'GR_E_V46::All' #rereco
#process.GlobalTag.globaltag  = 'GR_P_V50::All' #rereco
#process.GlobalTag.globaltag  = 'GR_P_V54::All' #rereco
#process.GlobalTag.globaltag  = 'GR_P_V56::All' #rereco
#process.GlobalTag.globaltag  = '74X_dataRun2_Prompt_v0'
#process.GlobalTag.globaltag  = '74X_dataRun2_Express_v0'
#process.GlobalTag.globaltag  = 'auto:run2_data'
#process.GlobalTag.globaltag  = '80X_dataRun2_Prompt_v8'
#from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_data', '')
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')



#
# message logger
#
process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.debugModules.append('rpcunpacker')
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1000)
process.MessageLogger.suppressWarning  = cms.untracked.vstring('Geometry', 'AfterSource','L1T','L1GlobalTriggerRawToDigi')
process.options = cms.untracked.PSet( wantSummary=cms.untracked.bool(False))

#
# fixes to R2D and emulation
#
#process.RawToDigi.remove(process.siPixelDigis)
#process.RawToDigi.remove(process.siStripDigis)
#process.RawToDigi.remove(process.ecalDigis)
#process.RawToDigi.remove(process.ecalPreshowerDigis)
#process.RawToDigi.remove(process.hcalDigis)
#process.RawToDigi.remove(process.castorDigis)
#process.RawToDigi.remove(process.tcdsDigis)
#process.RawToDigi.remove(process.caloStage1Digis)
#process.RawToDigi.remove(process.caloStage1FinalDigis)
#process.RawToDigi.remove(process.caloStage1LegacyFormatDigis)
#process.RawToDigi.remove(process.caloStage2Digis)
#process.RawToDigi.remove(process.gtStage2Digis)
#process.RawToDigi.remove(process.gtDigis)
#process.RawToDigi.remove(process.muonDTDigis)
#process.RawToDigi.remove(process.scalersRawToDigi)
#process.RawToDigi.remove(process.gctDigis)


#OMTF ESProducer. Fills CondFormats from XML files.
process.omtfParamsSource = cms.ESSource( "EmptyESSource",
     recordName = cms.string('L1TMuonOverlapParamsRcd'),
    iovIsRunNotTime = cms.bool(True),
    firstValid = cms.vuint32(1)
)

process.omtfParams = cms.ESProducer( "L1TMuonOverlapParamsESProducer",
     patternsXMLFiles = cms.VPSet( cms.PSet(patternsXMLFile = cms.FileInPath("L1Trigger/L1TMuon/data/omtf_config/Patterns_0x0003.xml")),),
     configXMLFile = cms.FileInPath("L1Trigger/L1TMuon/data/omtf_config/hwToLogicLayer_0x0003.xml"),
)

###OMTF emulator configuration
process.omtfEmulator = cms.EDProducer("L1TMuonOverlapTrackProducer",
#  srcDTPh =  cms.InputTag('simTwinMuxDigis'),
#  srcDTTh =  cms.InputTag('simTwinMuxDigis'),
  srcDTPh =  cms.InputTag('bmtfDigis'),
  srcDTTh =  cms.InputTag('bmtfDigis'),
#  srcDTPh =  cms.InputTag('bmtfDigis:PhiDigis'),
#  srcDTTh =  cms.InputTag('bmtfDigis:TheDigis'),
#  srcDTPh = cms.InputTag('simDtTriggerPrimitiveDigis'),
#  srcDTTh = cms.InputTag('simDtTriggerPrimitiveDigis'),
  srcCSC = cms.InputTag('csctfDigis'),
#  srcCSC = cms.InputTag('emtfStage2Digis'),
#  srcCSC = cms.InputTag('muonCSCDigis','MuonCSCCorrelatedLCTDigi'),
#  srcCSC = cms.InputTag('simCscTriggerPrimitiveDigis','MPCSORTED'),
  srcRPC = cms.InputTag('muonRPCDigis'),
#  srcRPC = cms.InputTag('simMuonRPCDigis'),
  dumpResultToXML = cms.bool(False),
  dumpDetailedResultToXML = cms.bool(False),
  XMLDumpFileName = cms.string("TestEvents.xml"),
  dumpGPToXML = cms.bool(True),
  readEventsFromXML = cms.bool(False),
  eventsXMLFiles = cms.vstring("TestEvents.xml"),
  dropRPCPrimitives = cms.bool(False),
  dropDTPrimitives = cms.bool(False),
  dropCSCPrimitives = cms.bool(False),
)

process.raw2digi_step = cms.Path(process.muonRPCDigis+process.csctfDigis+process.bmtfDigis+process.emtfStage2Digis+process.gmtStage2Digis)
#process.raw2digi_step = cms.Path(process.RawToDigi+process.emtfStage2Digis)
process.omtf_step = cms.Path(process.omtfEmulator)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.schedule = cms.Schedule(process.raw2digi_step, process.omtf_step, process.endjob_step)

#
# re-emulate L1T
#
#from L1Trigger.Configuration.customiseReEmul import L1TReEmulFromRAW
#process = L1TReEmulFromRAW(process)
#process.L1TReEmul.remove(process.simCaloStage2Digis)
#process.L1TReEmul.remove(process.simCaloStage2Layer1Digis)
#process.L1TReEmul.remove(process.simGmtCaloSumDigis)
#process.L1TReEmul.remove(process.simGmtStage2Digis)
#process.L1TReEmul.remove(process.simEcalTriggerPrimitiveDigis)

#process.simOmtfDigis.srcDTPh = cms.InputTag('BMTFStage2Digis')
#process.simOmtfDigis.srcDTTh = cms.InputTag('BMTFStage2Digis')
#process.simOmtfDigis.dumpResultToXML = cms.bool(True)

process.omtfTree = cms.EDAnalyzer("OmtfTreeMaker",
  histoFileName = cms.string("omtfHelper.root"),
  treeFileName = cms.string("omtfTree.root"),
#   omtfEmulSrc = cms.InputTag('simOmtfDigis','OMTF'),
  omtfEmulSrc = cms.InputTag('omtfEmulator','OMTF'),
  omtfDataSrc = cms.InputTag('gmtStage2Digis','OMTF'),
)

process.OmtfTree = cms.Path(process.omtfTree)
process.schedule.append(process.OmtfTree)

#print process.dumpPython();
print process.schedule



