#include "UserCode/OmtfAnalysis/interface/AnaTime.h"

#include "TObjArray.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TEfficiency.h"
#include "UserCode/OmtfDataFormats/interface/TrackObj.h"
#include "UserCode/OmtfDataFormats/interface/MuonObj.h"
#include "UserCode/OmtfDataFormats/interface/MuonObjColl.h"
#include "UserCode/OmtfDataFormats/interface/EventObj.h"
#include "UserCode/OmtfDataFormats/interface/L1ObjColl.h"
#include "UserCode/OmtfAnalysis/interface/Utilities.h"
#include "DataFormats/Math/interface/deltaR.h"


#include <ostream>
#include <iostream>
#include <cmath>
#include <bitset>


namespace { 
  TH1D *hTimeOmtf_A,  *hTimeBmtf_A,   *hTimeEmtf_A;  
  TH1D *hTimeOmtf_Q,  *hTimeBmtf_Q,   *hTimeEmtf_Q;  
  TH1D *hTimeOmtf_M,  *hTimeBmtf_M,   *hTimeEmtf_M;  
  TH1D *hTimeOmtf_QM, *hTimeBmtf_QM,  *hTimeEmtf_QM;  
  TH1D *hTimeOmtf_W,  *hTimeBmtf_W,  *hTimeEmtf_W;  
  TH1D *hTimeOmtf_QW, *hTimeBmtf_QW, *hTimeEmtf_QW;  
  TH1D *hTimeOmtf_emu_A, *hTimeOmtf_emu_Q, *hTimeOmtf_emu_M, *hTimeOmtf_emu_QM, *hTimeOmtf_emu_W, *hTimeOmtf_emu_QW;

  TH2D *hTimeBmtfOmtf, *hTimeOmtfEmtf, *hTimeOmtfOmtf_E;

  TH1D *hTimeDeltaR_Q, *hTimeDeltaR_QW, 
       *hTimeDeltaR_Q_B1, *hTimeDeltaR_Q_B2, *hTimeDeltaR_Q_B3, 
       *hTimeDeltaR_QW_B1, *hTimeDeltaR_QW_B2, *hTimeDeltaR_QW_B3;
  TH1D *hTimeLayers;

  TH2D *hTimeOmtfTrackDPhiT, *hTimeOmtfTrackDPhiM, *hTimeOmtfTrackDEtaT, *hTimeOmtfTrackDEtaM;
  TH1D *hTimeOmtfTrackDRM;

  TH2D *hTimeOmtfTrackBXT, *hTimeOmtfTrackBXM;
  TH1D *hTimeOmtfTrackBX0, *hTimeOmtfTrackBX1;
  TH2D *hTimeOmtfDrTrackMuon;

  TEfficiency *hTimeEffPt_BMTF, *hTimeEffPt_EMTF, *hTimeEffPt_OMTF, *hTimeEffPt_OMTF_emu;
  TEfficiency *hTimeEta_Pt0, *hTimeEta_Pt10, *hTimeEta_Pt22;


/////////////////////////////////////moja praca
  TEfficiency *hTimePrefireEta;
  TEfficiency *hTimePrefireEtaOMTF;
  TEfficiency *hTimePrefireEtaBMTF;
  TEfficiency *hTimePrefireEtaEMTF;
  TEfficiency *hTimePrefireEta1OMTF;
  TEfficiency *hTimePrefireEta1BMTF;
  TEfficiency *hTimePrefireEta1EMTF;
  TEfficiency *hTimePrefireEta1;
  TEfficiency *hTimePrefireEtahitpattern;
  TEfficiency *hTimeEffPt1;

  TH2D *hDeltaphieta010;
  TH2D *hDeltaphieta1022;
  TH2D *hDeltaphieta22inf;
  TH2D *hDeltaphieta010_025;
  TH2D *hDeltaphieta1022_025;
  TH2D *hDeltaphieta22inf_025;

  TEfficiency *hTimeEffPt1_BMTF, *hTimeEffPt1_EMTF, *hTimeEffPt1_OMTF, *hTimeEffPt1_OMTF_emu;

  TH2D *hQualitybx0;
  //TH2D *hPtbx;

  TH1D * hTimePt22;
  TH2D * hTimePtbinning;
  TH2D * hTimePtbinninghitpattern;

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

}

AnaTime::AnaTime(const edm::ParameterSet& cfg)
  : debug(false), theCfg (cfg)
{}

void AnaTime::init(TObjArray& histos)
{
  hTimeBmtf_A = new TH1D("hTimeBmtf_A","hTimeBmtf_A",7,-3.5,3.5); histos.Add(hTimeBmtf_A); 
  hTimeBmtf_Q = new TH1D("hTimeBmtf_Q","hTimeBmtf_Q",7,-3.5,3.5); histos.Add(hTimeBmtf_Q); 
  hTimeBmtf_M = new TH1D("hTimeBmtf_M","hTimeBmtf_M",7,-3.5,3.5); histos.Add(hTimeBmtf_M); 
  hTimeBmtf_QM = new TH1D("hTimeBmtf_QM","hTimeBmtf_QM",7,-3.5,3.5); histos.Add(hTimeBmtf_QM); 
  hTimeBmtf_W = new TH1D("hTimeBmtf_W","hTimeBmtf_W",7,-3.5,3.5); histos.Add(hTimeBmtf_W); 
  hTimeBmtf_QW = new TH1D("hTimeBmtf_QW","hTimeBmtf_QW",7,-3.5,3.5); histos.Add(hTimeBmtf_QW); 

  hTimeEmtf_A = new TH1D("hTimeEmtf_A","hTimeEmtf_A",7,-3.5,3.5); histos.Add(hTimeEmtf_A); 
  hTimeEmtf_Q = new TH1D("hTimeEmtf_Q","hTimeEmtf_Q",7,-3.5,3.5); histos.Add(hTimeEmtf_Q); 
  hTimeEmtf_M = new TH1D("hTimeEmtf_M","hTimeEmtf_M",7,-3.5,3.5); histos.Add(hTimeEmtf_M); 
  hTimeEmtf_QM = new TH1D("hTimeEmtf_QM","hTimeEmtf_QM",7,-3.5,3.5); histos.Add(hTimeEmtf_QM); 
  hTimeEmtf_W = new TH1D("hTimeEmtf_W","hTimeEmtf_W",7,-3.5,3.5); histos.Add(hTimeEmtf_W); 
  hTimeEmtf_QW = new TH1D("hTimeEmtf_QW","hTimeEmtf_QW",7,-3.5,3.5); histos.Add(hTimeEmtf_QW); 

  hTimeOmtf_A = new TH1D("hTimeOmtf_A","hTimeOmtf_A",7,-3.5,3.5); histos.Add(hTimeOmtf_A); 
  hTimeOmtf_Q = new TH1D("hTimeOmtf_Q","hTimeOmtf_Q",7,-3.5,3.5); histos.Add(hTimeOmtf_Q); 
  hTimeOmtf_M = new TH1D("hTimeOmtf_M","hTimeOmtf_M",7,-3.5,3.5); histos.Add(hTimeOmtf_M); 
  hTimeOmtf_QM = new TH1D("hTimeOmtf_QM","hTimeOmtf_QM",7,-3.5,3.5); histos.Add(hTimeOmtf_QM); 
  hTimeOmtf_W = new TH1D("hTimeOmtf_W","hTimeOmtf_W",7,-3.5,3.5); histos.Add(hTimeOmtf_W); 
  hTimeOmtf_QW = new TH1D("hTimeOmtf_QW","hTimeOmtf_QW",7,-3.5,3.5); histos.Add(hTimeOmtf_QW); 

  hTimeOmtf_emu_A = new TH1D("hTimeOmtf_emu_A","hTimeOmtf_emu_A",7,-3.5,3.5); histos.Add(hTimeOmtf_emu_A); 
  hTimeOmtf_emu_Q = new TH1D("hTimeOmtf_emu_Q","hTimeOmtf_emu_Q",7,-3.5,3.5); histos.Add(hTimeOmtf_emu_Q); 
  hTimeOmtf_emu_M = new TH1D("hTimeOmtf_emu_M","hTimeOmtf_emu_M",7,-3.5,3.5); histos.Add(hTimeOmtf_emu_M); 
  hTimeOmtf_emu_QM = new TH1D("hTimeOmtf_emu_QM","hTimeOmtf_emu_QM",7,-3.5,3.5); histos.Add(hTimeOmtf_emu_QM); 
  hTimeOmtf_emu_W = new TH1D("hTimeOmtf_emu_W","hTimeOmtf_emu_W",7,-3.5,3.5); histos.Add(hTimeOmtf_emu_W); 
  hTimeOmtf_emu_QW = new TH1D("hTimeOmtf_emu_QW","hTimeOmtf_emu_QW",7,-3.5,3.5); histos.Add(hTimeOmtf_emu_QW); 

  hTimeDeltaR_Q = new TH1D("hTimeDeltaR_Q","hTimeDeltaR_Q",   100,0.,4.);  histos.Add(hTimeDeltaR_Q);
  hTimeDeltaR_QW = new TH1D("hTimeDeltaR_QW","hTimeDeltaR_QW",100,0.,4.);  histos.Add(hTimeDeltaR_QW);
  hTimeDeltaR_Q_B3 = new TH1D("hTimeDeltaR_Q_B3","hTimeDeltaR_Q_B3",100,0.,4.);  histos.Add(hTimeDeltaR_Q_B3);
  hTimeDeltaR_Q_B2 = new TH1D("hTimeDeltaR_Q_B2","hTimeDeltaR_Q_B2",100,0.,4.);  histos.Add(hTimeDeltaR_Q_B2);
  hTimeDeltaR_Q_B1 = new TH1D("hTimeDeltaR_Q_B1","hTimeDeltaR_Q_B1",100,0.,4.);  histos.Add(hTimeDeltaR_Q_B1);
  hTimeDeltaR_QW_B3 = new TH1D("hTimeDeltaR_QW_B3","hTimeDeltaR_QW_B3",100,0.,4.);  histos.Add(hTimeDeltaR_QW_B3);
  hTimeDeltaR_QW_B2 = new TH1D("hTimeDeltaR_QW_B2","hTimeDeltaR_QW_B2",100,0.,4.);  histos.Add(hTimeDeltaR_QW_B2);
  hTimeDeltaR_QW_B1 = new TH1D("hTimeDeltaR_QW_B1","hTimeDeltaR_QW_B1",100,0.,4.);  histos.Add(hTimeDeltaR_QW_B1);



  unsigned int nOmtfLayers =  omtfUtilities::layerNames.size();
  hTimeLayers =  new TH1D("hTimeLayers","hTimeLayers", nOmtfLayers,-0.5,nOmtfLayers-0.5); histos.Add(hTimeLayers);
  for (unsigned int ibin=1; ibin<=nOmtfLayers; ibin++) hTimeLayers->GetXaxis()->SetBinLabel(ibin,omtfUtilities::layerNames[ibin-1].c_str());



  hTimeBmtfOmtf = new TH2D("hTimeBmtfOmtf","hTimeBmtfOmtf",5,-2.5,2.5, 5,-2.5,2.5); histos.Add(hTimeBmtfOmtf);
  hTimeOmtfEmtf = new TH2D("hTimeOmtfEmtf","hTimeOmtfEmtf",5,-2.5,2.5, 5,-2.5,2.5); histos.Add(hTimeOmtfEmtf);
  hTimeOmtfOmtf_E= new TH2D("hTimeOmtfOmtf_E","hTimeOmtfOmtf_E",5,-2.5,2.5, 5,-2.5,2.5); histos.Add(hTimeOmtfOmtf_E);

  hTimeOmtfTrackDPhiT = new TH2D("hTimeOmtfTrackDPhiT","hTimeOmtfTrackDPhiT",50,0.,25., 50, -1.,1.); histos.Add(hTimeOmtfTrackDPhiT);
  hTimeOmtfTrackDPhiM = new TH2D("hTimeOmtfTrackDPhiM","hTimeOmtfTrackDPhiM",50,0.,25., 50, -1.,1.); histos.Add(hTimeOmtfTrackDPhiM);
  hTimeOmtfTrackDEtaT = new TH2D("hTimeOmtfTrackDEtaT","hTimeOmtfTrackDEtaT",50,0.,25., 50, -1.,1.); histos.Add(hTimeOmtfTrackDEtaT);
  hTimeOmtfTrackDEtaM = new TH2D("hTimeOmtfTrackDEtaM","hTimeOmtfTrackDEtaM",50,0.,25., 50, -1.,1.); histos.Add(hTimeOmtfTrackDEtaM);
//  hTimeOmtfTrackDRM   = new TH2D("hTimeOmtfTrackDRM",  "hTimeOmtfTrackDRM",  50,0.,25., 50, -1.,1.); histos.Add(hTimeOmtfTrackDRM);
  hTimeOmtfTrackDRM   = new TH1D("hTimeOmtfTrackDRM",  "hTimeOmtfTrackDRM",  50, 0.,2.); histos.Add(hTimeOmtfTrackDRM);
  hTimeOmtfTrackBXT = new TH2D("hTimeOmtfTrackBXT", "hTimeOmtfTrackBXT", 50,0.,25.,4, -1.,3.); histos.Add(hTimeOmtfTrackBXT);
  hTimeOmtfTrackBXM = new TH2D("hTimeOmtfTrackBXM", "hTimeOmtfTrackBXM", 50,0.,25.,4, -1.,3.); histos.Add(hTimeOmtfTrackBXM);
  hTimeOmtfTrackBX0 = new TH1D("hTimeOmtfTrackBX0", "hTimeOmtfTrackBX0", 50,0.,25.); histos.Add(hTimeOmtfTrackBX0);
  hTimeOmtfTrackBX1 = new TH1D("hTimeOmtfTrackBX1", "hTimeOmtfTrackBX1", 50,0.,25.); histos.Add(hTimeOmtfTrackBX1);
  hTimeOmtfDrTrackMuon = new TH2D("hTimeOmtfDrTrackMuon","hTimeOmtfDrTrackMuon",50,0.,25.,50, -1.,1.); histos.Add(hTimeOmtfDrTrackMuon);

  double xbins[]={0.,4., 8., 12., 16., 22., 30.};
  hTimeEffPt_BMTF = new TEfficiency("hTimeEffPt_BMTF","BMTF: bx=-1/(bx=-1 or bx=0); L1 p_{T}; fraction",6,xbins); histos.Add(hTimeEffPt_BMTF);
  hTimeEffPt_EMTF = new TEfficiency("hTimeEffPt_EMTF","EMTF: bx=-1/(bx=-1 or bx=0); L1 p_{T}; fraction",6,xbins); histos.Add(hTimeEffPt_EMTF);
  hTimeEffPt_OMTF = new TEfficiency("hTimeEffPt_OMTF","OMTF: bx=-1/(bx=-1 or bx=0); L1 p_{T}; fraction",6,xbins); histos.Add(hTimeEffPt_OMTF);
  hTimeEffPt_OMTF_emu = new TEfficiency("hTimeEffPt_OMTF_emu","hTimeEffPt_OMTF_emu: BX=-1,-2/(BX==-1,-2 or BX==0); L1 p_{T}; fraction",6,xbins); histos.Add(hTimeEffPt_OMTF_emu);

  double etas[]={-2.4,-2.0,-1.6,-1.24,-0.83,-0.5,-0.15,0.15,0.5,0.83,1.24,1.6,2.0,2.4};
//  hTimeEta_Pt0  = new TEfficiency( "hTimeEta_Pt0", "hTimeEta_Pt0",24, -2.4, 2.4); histos.Add(hTimeEta_Pt0);
//  hTimeEta_Pt10 = new TEfficiency("hTimeEta_Pt10","hTimeEta_Pt10",24, -2.4, 2.4); histos.Add(hTimeEta_Pt10);
//  hTimeEta_Pt22 = new TEfficiency("hTimeEta_Pt22","hTimeEta_Pt22",24, -2.4, 2.4); histos.Add(hTimeEta_Pt22);
  hTimeEta_Pt0  = new TEfficiency( "hTimeEta_Pt0", "hTimeEta_Pt0",13,etas); histos.Add(hTimeEta_Pt0);
  hTimeEta_Pt10 = new TEfficiency("hTimeEta_Pt10","hTimeEta_Pt10",13,etas); histos.Add(hTimeEta_Pt10);
  hTimeEta_Pt22 = new TEfficiency("hTimeEta_Pt22","hTimeEta_Pt22",13,etas); histos.Add(hTimeEta_Pt22);

//////////////////////////////  /////////////////////////////////////// moja praca
  hTimePrefireEta = new TEfficiency("hTimePrefireEta","hTimePrefireEta; #eta ;p",50,-2.5,2.5); histos.Add(hTimePrefireEta);
  hTimePrefireEtahitpattern = new TEfficiency("hTimePrefireEtahitpattern","hTimePrefireEtahitpattern; #eta ;p",50,-2.5,2.5); histos.Add(hTimePrefireEtahitpattern);
  hTimePrefireEtaOMTF = new TEfficiency("hTimePrefireEtaOMTF","hTimePrefireEtaOMTF; #eta ;p",50,-2.5,2.5); histos.Add(hTimePrefireEtaOMTF);
  hTimePrefireEtaBMTF = new TEfficiency("hTimePrefireEtaBMTF","hTimePrefireEtaBMTF; #eta ;p",50,-2.5,2.5); histos.Add(hTimePrefireEtaBMTF);
  hTimePrefireEtaEMTF = new TEfficiency("hTimePrefireEtaEMTF","hTimePrefireEtaEMTF; #eta ;p",50,-2.5,2.5); histos.Add(hTimePrefireEtaEMTF);
  hTimePrefireEta1 = new TEfficiency("hTimePrefireEta1","hTimePrefireEta1; #eta ;p",50,-2.5,2.5); histos.Add(hTimePrefireEta1);
  hTimePrefireEta1OMTF = new TEfficiency("hTimePrefireEta1OMTF","hTimePrefireEta1OMTF; #eta ;p",50,-2.5,2.5); histos.Add(hTimePrefireEta1OMTF);
  hTimePrefireEta1BMTF = new TEfficiency("hTimePrefireEta1BMTF","hTimePrefireEta1BMTF; #eta ;p",50,-2.5,2.5); histos.Add(hTimePrefireEta1BMTF);
  hTimePrefireEta1EMTF = new TEfficiency("hTimePrefireEta1EMTF","hTimePrefireEta1EMTF; #eta ;p",50,-2.5,2.5); histos.Add(hTimePrefireEta1EMTF);
  hTimeEffPt1 = new TEfficiency("hTimeEffPt1 ","hTimeEffPt1 : ); L1 p_{T}; fraction",6,xbins); histos.Add(hTimeEffPt1 );


  hDeltaphieta010=new TH2D("hDeltaphieta010","hDeltaphieta010; #Delta #eta; #Delta #phi",10,0.,2.5, 10,0.,2*M_PI); histos.Add(hDeltaphieta010);
  hDeltaphieta1022=new TH2D("hDeltaphieta1022","hDeltaphieta1022; #Delta #eta; #Delta #phi",10,0.,2.5, 10,0.,2*M_PI); histos.Add(hDeltaphieta1022);
  hDeltaphieta22inf=new TH2D("hDeltaphieta22inf","hDeltaphieta22inf; #Delta #eta; #Delta #phi",10,0.,2.5, 10,0.,2*M_PI); histos.Add(hDeltaphieta22inf);
  hDeltaphieta010_025=new TH2D("hDeltaphieta010_025","hDeltaphieta010_025; #Delta #eta; #Delta #phi",10,0.,0.25, 10,0.,0.25); histos.Add(hDeltaphieta010_025);
  hDeltaphieta1022_025=new TH2D("hDeltaphieta1022_025","hDeltaphieta1022_025; #Delta #eta; #Delta #phi",10,0.,0.25, 10,0.,0.25); histos.Add(hDeltaphieta1022_025);
  hDeltaphieta22inf_025=new TH2D("hDeltaphieta22inf_025","hDeltaphieta22inf_025; #Delta #eta; #Delta #phi",10,0.,0.25, 10,0.,0.25); histos.Add(hDeltaphieta22inf_025);


  hQualitybx0= new TH2D("hQualitybx0","OMTF quality at bx0; Quality; p_{T}",13,0.,13., 3,0.,3.); histos.Add(hQualitybx0);
  hTimePt22= new TH1D("hTimePt22","OMTF bx0 and bx-1; bx; counts",2,0.,2.);  histos.Add(hTimePt22);
  hTimePtbinning= new TH2D("hTimePtbinning","OMTF bx0 and bx-1; bx; p_{T}",2,0.,2., 3,0.,3.); histos.Add(hTimePtbinning);
  hTimePtbinninghitpattern= new TH2D("hTimePtbinninghitpattern","OMTF bx0 and bx-1; bx; p_{T}",2,0.,2., 3,0.,3.); histos.Add(hTimePtbinninghitpattern);

  hTimeEffPt1_BMTF = new TEfficiency("hTimeEffPt1_BMTF","BMTF: bx=-1/(bx=-1 or bx=0); L1 p_{T}; fraction",6,xbins); histos.Add(hTimeEffPt1_BMTF);
  hTimeEffPt1_EMTF = new TEfficiency("hTimeEffPt1_EMTF","EMTF: bx=-1/(bx=-1 or bx=0); L1 p_{T}; fraction",6,xbins); histos.Add(hTimeEffPt1_EMTF);
  hTimeEffPt1_OMTF = new TEfficiency("hTimeEffPt1_OMTF","OMTF: bx=-1/(bx=-1 or bx=0); L1 p_{T}; fraction",6,xbins); histos.Add(hTimeEffPt1_OMTF);
  hTimeEffPt1_OMTF_emu = new TEfficiency("hTimeEffPt1_OMTF_emu","hTimeEffPt1_OMTF_emu: BX=-1,-2/(BX==-1,-2 or BX==0); L1 p_{T}; fraction",6,xbins); histos.Add(hTimeEffPt1_OMTF_emu);


  // float xbin[]={0,10,22,1000};
  //hPtbx = new TH2D("hPtbx","OMTF quality at bx0; Quality; p_{T}",xbin, 2,-1.,0.); histos.Add(hPtbx);
}

void AnaTime::run(const EventObj* ev, const MuonObjColl *muonColl, const TrackObj* track, const L1ObjColl * l1Objs)
{
  //std::cout<<__LINE__<<""<<ev->run<<" number: "<<ev->id<<std::endl;
  //if(ev->run!=366832 || ev->id!=218059960) return;
  //std::cout<<"*****************************"<<std::endl;
  std::vector<L1Obj::TYPE> mtfs= {L1Obj::BMTF, L1Obj::OMTF, L1Obj::EMTF, L1Obj::OMTF_emu};

  //
  // all triggers
  //
  const std::vector<L1Obj> & l1mtfs = *l1Objs;
  const std::vector<MuonObj> & muons = *muonColl;

  bool printdeb = false;
  double deltaRMatching = theCfg.getParameter<double>("deltaRMatching");

  //
  // find if muon has corresponding triggering L1 -> muonsExt
  //
  std::vector<std::pair<MuonObj,bool> > muonsExt;
  for (const MuonObj & muon : muons) {
    bool hasTriggeringL1 = false;
    if (!muon.isValid()) continue;
    for (const auto & l1mtf : l1mtfs) {
      if (!l1mtf.isValid()) continue;
      if (l1mtf.q < 12) continue;
      if (l1mtf.ptValue() < theCfg.getParameter<double>("requireOtherMuTrgL1Pt")) continue;
      if (l1mtf.bx != 0) continue;
      double deltaR = reco::deltaR( l1mtf.etaValue(), l1mtf.phiValue(), muon.l1Eta, muon.l1Phi);
      if (deltaR < deltaRMatching) hasTriggeringL1 = true;
    }
    muonsExt.push_back(std::make_pair(muon,hasTriggeringL1));
  }    


//
// for each L1 
//


  std::bitset<18>   checkDT(std::string("000000000000111111"));
  std::bitset<18> checkREST(std::string("111111111111000000"));
  L1Obj omtfBXm1;
  bool omtfBXm1b=false;
  for (const auto & l1mtf : l1mtfs) {
    bool matched  = false;
    bool matchedW = false;
    bool qualOK = (l1mtf.q >= 12);
    bool ptOK = (l1mtf.ptValue()< theCfg.getParameter<double>("maxPtForDistributions") && l1mtf.ptValue() >= theCfg.getParameter<double>("minPtForDistributions") );
    bool hasRequiredTrigger = theCfg.exists("requireOtherMuTrg") ? !theCfg.getParameter<bool>("requireOtherMuTrg") : true;


    for (const auto & muonExt : muonsExt) {
      const MuonObj & muon = muonExt.first;
      bool muon_associatedTriggeringL1 = muonExt.second; 
      if (!muon.isValid()) continue;

      double deltaR = reco::deltaR( l1mtf.etaValue(), l1mtf.phiValue(), muon.l1Eta, muon.l1Phi);
//    double deltaRW = reco::deltaR( l1mtf.etaValue(), l1mtf.phiValue(), -muon.l1Eta, muon.l1Phi+M_PI/2.);
      double deltaRW = reco::deltaR(-l1mtf.etaValue(), l1mtf.phiValue()+M_PI/2., muon.l1Eta, muon.l1Phi);
      if (qualOK && ptOK && l1mtf.type==L1Obj::OMTF) {
        hTimeDeltaR_Q->Fill(deltaR); 
        hTimeDeltaR_QW->Fill(deltaRW); 
        if (l1mtf.bx == -3) hTimeDeltaR_Q_B3->Fill(deltaR); 
        if (l1mtf.bx == -2) hTimeDeltaR_Q_B2->Fill(deltaR); 
        if (l1mtf.bx == -1) hTimeDeltaR_Q_B1->Fill(deltaR); 
        if (l1mtf.bx == -3) hTimeDeltaR_QW_B3->Fill(deltaRW); 
        if (l1mtf.bx == -2) hTimeDeltaR_QW_B2->Fill(deltaRW); 
        if (l1mtf.bx == -1) hTimeDeltaR_QW_B1->Fill(deltaRW); 
        if ( deltaR < 0.2 && l1mtf.bx < 0 && l1mtf.type==L1Obj::OMTF) {
//      if ( deltaR < 0.2 && l1mtf.bx == 0 && l1mtf.type==L1Obj::OMTF) {
          printdeb = false;
          std::bitset<18> hitLayers(l1mtf.hits);
//          for (unsigned int hitLayer=0; hitLayer<18;hitLayer++) if(hitLayers[hitLayer]) hTimeLayers->Fill(hitLayer);
          for(unsigned int hitLayer=0;hitLayer<18;hitLayer++){
            int check=1;
            if((check<<hitLayer)&l1mtf.hits) hTimeLayers->Fill(hitLayer);
          }
        }
      }
      if (deltaR  < deltaRMatching) matched=true; 
      if (deltaRW < deltaRMatching) matchedW=true;
      if (deltaR  > 2.*deltaRMatching && (muon.isMatchedHlt || muon.isMatchedIsoHlt) && muon_associatedTriggeringL1) hasRequiredTrigger = true; 
    }
    if (!hasRequiredTrigger) continue;

    TH1D *hA, *hQ, *hM, *hQM, *hW, *hQW; 
    hA=hQ=hM=hQM=hW=hQW=0; 
    TEfficiency * hE=0;
    TEfficiency * hE1=0;
    switch (l1mtf.type) {
        case (L1Obj::BMTF) : hA=hTimeBmtf_A; hQ=hTimeBmtf_Q; hM=hTimeBmtf_M;  hQM=hTimeBmtf_QM;  hW=hTimeBmtf_W;  hQW=hTimeBmtf_QW; hE=hTimeEffPt_BMTF; hE1=hTimeEffPt1_BMTF; break;
        case (L1Obj::EMTF) : hA=hTimeEmtf_A; hQ=hTimeEmtf_Q; hM=hTimeEmtf_M;  hQM=hTimeEmtf_QM;  hW=hTimeEmtf_W;  hQW=hTimeEmtf_QW; hE=hTimeEffPt_EMTF; hE1=hTimeEffPt1_EMTF; break;
        case (L1Obj::OMTF) : hA=hTimeOmtf_A; hQ=hTimeOmtf_Q; hM=hTimeOmtf_M;  hQM=hTimeOmtf_QM;  hW=hTimeOmtf_W;  hQW=hTimeOmtf_QW; hE=hTimeEffPt_OMTF; hE1=hTimeEffPt1_OMTF; break;
        case (L1Obj::OMTF_emu) : hA=hTimeOmtf_emu_A; hQ=hTimeOmtf_emu_Q; hM=hTimeOmtf_emu_M;  hQM=hTimeOmtf_emu_QM;  hW=hTimeOmtf_emu_W;  hQW=hTimeOmtf_emu_QW; hE=hTimeEffPt_OMTF_emu;hE1=hTimeEffPt1_OMTF_emu; break;
        default: ;
    }
    if (hA!=0 && ptOK) {
      hA->Fill(l1mtf.bx); 
      if (qualOK) hQ->Fill(l1mtf.bx); 
      if (matched) hM->Fill(l1mtf.bx);  
      if (qualOK && matched) hQM->Fill(l1mtf.bx);  
      if (matchedW) hW->Fill(l1mtf.bx);  
      if (qualOK && matchedW) hQW->Fill(l1mtf.bx);  
    }
//  bool pref = (l1mtf.bx == -1 || l1mtf.bx == -2);
    bool pref = (l1mtf.bx == -1);
    bool pref1 = (l1mtf.bx == -1);
    if (qualOK && matched && (pref || l1mtf.bx == 0)) {
      for(const auto & l1mtf2 : l1mtfs){
        if(l1mtf2.bx==0 && l1mtf2.type==L1Obj::OMTF && l1mtf.type==L1Obj::OMTF && l1mtf.position==l1mtf2.position && l1mtf.iProcessor==l1mtf2.iProcessor && (std::abs(reco::deltaPhi(l1mtf.phiValue(),l1mtf2.phiValue()))<0.09)){
            pref1=false;
         //   std::cout<<"CZEMU TU WCHODZISZ"<<std::endl;
        }
      }
      double ptValue = l1mtf.ptValue()> 25. ? 25. : l1mtf.ptValue();
      if (hE) { hE->Fill(pref, ptValue); }

      double ptValue1 = l1mtf.ptValue()> 25. ? 25. : l1mtf.ptValue();
      if (hE1) { hE1->Fill(pref1, ptValue1); }
//      if (pref && ptValue>=22 && l1mtf.type==L1Obj::OMTF) printdeb=true;
//    if (l1mtf.type==L1Obj::uGMT) {
      if (l1mtf.type==L1Obj::OMTF || l1mtf.type==L1Obj::BMTF || l1mtf.type==L1Obj::EMTF) {
        if (l1mtf.ptValue()<10) hTimeEta_Pt0->Fill(pref, l1mtf.etaValue());
        else if (l1mtf.ptValue()<22.) hTimeEta_Pt10->Fill(pref, l1mtf.etaValue());
        else hTimeEta_Pt22->Fill(pref, l1mtf.etaValue());
      }
    }

    bool checkprefire=false;
    if(qualOK && l1mtf.bx==-1){
      checkprefire=true;
    }
    if(l1mtf.ptValue()>10 && (l1mtf.type==L1Obj::OMTF || l1mtf.type==L1Obj::BMTF || l1mtf.type==L1Obj::EMTF)) {
      hTimePrefireEta->Fill(checkprefire,l1mtf.etaValue());
      if(l1mtf.type==L1Obj::OMTF) hTimePrefireEtaOMTF->Fill(checkprefire,l1mtf.etaValue());
      if(l1mtf.type==L1Obj::BMTF) hTimePrefireEtaBMTF->Fill(checkprefire,l1mtf.etaValue());
      if(l1mtf.type==L1Obj::EMTF) hTimePrefireEtaEMTF->Fill(checkprefire,l1mtf.etaValue());

    }

    bool checkprefirehitpattern=false;
    bool checkprefire1=false;
    std::bitset<18> l10bit(l1mtf.hits);
    std::bitset<18> l10checkDT=(checkDT & l10bit);
    std::bitset<18> l10checkREST=(checkREST & l10bit);
    if(qualOK && (l1mtf.type==L1Obj::OMTF || l1mtf.type==L1Obj::BMTF || l1mtf.type==L1Obj::EMTF)  && l1mtf.bx==-1){
      checkprefire1=true;
      checkprefirehitpattern=true;
    //  std::cout<<"CZEMU TU  nieeeeeee WCHODZISZ"<<std::endl;
    }
    for(const auto & l1mtf2 : l1mtfs){
      std::bitset<18> l11bit(l1mtf2.hits);
      std::bitset<18> l11checkDT=(checkDT & l11bit);
      std::bitset<18> l11checkREST=(checkREST & l11bit);
      if(l1mtf2.bx==0 && l1mtf2.type==L1Obj::OMTF && l1mtf.type==L1Obj::OMTF && l1mtf.position==l1mtf2.position && l1mtf.iProcessor==l1mtf2.iProcessor && (std::abs(reco::deltaPhi(l1mtf.phiValue(),l1mtf2.phiValue()))<0.09)){
          checkprefire1=false;

       //   std::cout<<"CZEMU TU WCHODZISZ"<<std::endl;
      }
      if(l1mtf2.bx==0 && l1mtf2.type==L1Obj::OMTF && l1mtf.type==L1Obj::OMTF && l1mtf.position==l1mtf2.position && l1mtf.iProcessor==l1mtf2.iProcessor){
        if(l10checkDT.count()>0 && l10checkREST.count()==0 && l11checkDT.count()>0 && l11checkREST.count()>0)checkprefirehitpattern=false;
      }
    }
    if(l1mtf.ptValue()>10 && (l1mtf.type==L1Obj::OMTF || l1mtf.type==L1Obj::BMTF || l1mtf.type==L1Obj::EMTF)){

      hTimePrefireEta1->Fill(checkprefire1,l1mtf.etaValue());
      if(l1mtf.type==L1Obj::OMTF) hTimePrefireEta1OMTF->Fill(checkprefire1,l1mtf.etaValue());
      if(l1mtf.type==L1Obj::BMTF) hTimePrefireEta1BMTF->Fill(checkprefire1,l1mtf.etaValue());
      if(l1mtf.type==L1Obj::EMTF) hTimePrefireEta1EMTF->Fill(checkprefire1,l1mtf.etaValue());
      hTimePrefireEtahitpattern->Fill(checkprefirehitpattern,l1mtf.etaValue());

    }




//    if(checkprefire1!=checkprefire3 && std::abs(l1mtf.etaValue())>1.0 && std::abs(l1mtf.etaValue())<1.1){
//
//              std::cout<<"???????????????????????????????????????????????????????????????"<<std::endl;
//              std::cout<<""<<ev->run<<" number: "<<ev->id<<std::endl;
//              std::cout<<"//////////////////////////////////////////////////////////////////////////"<<std::endl;
//              std::cout<<*l1Objs<<std::endl;
//              std::cout<<"HISTO"<<std::endl;
//              std::cout<<l1mtf<<std::endl;
//              std::cout<<"WARTOSC checkfire: "<<checkprefire1<<" Wartosc eta: "<<l1mtf.etaValue()<<std::endl;
//              std::cout<<"//////////////////////////////////////////////////////////////////////////"<<std::endl;
//              std::cout<<"CHECKPREFIRE:  "<< checkprefire<<std::endl;
//              std::cout<<"CHECKPREFIRE1:  "<< checkprefire1<<std::endl;
//              std::cout<<"CHECKPREFIRE3:  "<< checkprefire3<<std::endl;
//    }




//   int counter=0;
   for(const auto & l1mtf2 : l1mtfs){
     double deltaphi=0;
     double deltaeta=0;
 //      bool qualOK2 = (l1mtf2.q >= 12);
     if(l1mtf.q >= 12  && l1mtf.type==L1Obj::OMTF && l1mtf2.type==L1Obj::OMTF && l1mtf.bx==-1 && l1mtf2.bx==0 ){
       if(l1mtf.ptValue()<10){
         deltaphi=std::abs(reco::deltaPhi(l1mtf.phiValue(),l1mtf2.phiValue()));
         deltaeta=std::abs(l1mtf.etaValue()-l1mtf2.etaValue());
         hDeltaphieta010->Fill(deltaeta,deltaphi);
         if(l1mtf.position==l1mtf2.position && l1mtf.iProcessor==l1mtf2.iProcessor)hDeltaphieta010_025->Fill(deltaeta,deltaphi);
       }
       if(l1mtf.ptValue()>=10 && l1mtf.ptValue()<22 ){
         deltaphi=std::abs(reco::deltaPhi(l1mtf.phiValue(),l1mtf2.phiValue()));
         deltaeta=std::abs(l1mtf.etaValue()-l1mtf2.etaValue());
         hDeltaphieta1022->Fill(deltaeta,deltaphi);
         if(l1mtf.position==l1mtf2.position && l1mtf.iProcessor==l1mtf2.iProcessor) hDeltaphieta1022_025->Fill(deltaeta,deltaphi);
       }
       if(l1mtf.ptValue()>22){
         deltaphi=std::abs(reco::deltaPhi(l1mtf.phiValue(),l1mtf2.phiValue()));
         deltaeta=std::abs(l1mtf.etaValue()-l1mtf2.etaValue());
         hDeltaphieta22inf->Fill(deltaeta,deltaphi);
         if(l1mtf.position==l1mtf2.position && l1mtf.iProcessor==l1mtf2.iProcessor)hDeltaphieta22inf_025->Fill(deltaeta,deltaphi);
       }

     }
     if(l1mtf.q >= 12 && l1mtf.type==L1Obj::OMTF && l1mtf2.type==L1Obj::OMTF && l1mtf.bx==-1 && l1mtf2.bx==0 && l1mtf.position==l1mtf2.position && l1mtf.iProcessor==l1mtf2.iProcessor){
      if(l1mtf.ptValue()<10) hQualitybx0->Fill(l1mtf2.q,0.);
      if(l1mtf.ptValue()>=10 && l1mtf.ptValue()<22) hQualitybx0->Fill(l1mtf2.q,1.);
      if(l1mtf.ptValue()>=22) hQualitybx0->Fill(l1mtf2.q,2.);
//      if(l1mtf.ptValue()>22){
//        if(counter==1 || counter==2){
//        std::cout<<"//////////////////////////////////////////////////////////////////////////"<<std::endl;
//        std::cout<<*l1Objs<<std::endl;
//        std::cout<<"//////////////////////////////////////////////////////////////////////////"<<std::endl;
//        }
//        counter++;
//       }
   }
   }
   L1Obj tmp;
   if(l1mtf.q>=12 && l1mtf.bx==-1 && l1mtf.type==L1Obj::OMTF){
     tmp=l1mtf;
   }
   if(tmp.ptValue()>omtfBXm1.ptValue()){
     omtfBXm1=tmp;
     omtfBXm1b=true;
   }


//    if(qualOK && l1mtf.ptValue()> 22 && (l1mtf.bx==-1)){
//
//        std::cout<<"Tutaj jest event"<<std::endl;
//        std::cout << *ev << std::endl;
//        std::cout<<"//////////////////"<<std::endl;
//
//        std::cout<<"Tutaj jest muonColl"<<std::endl;
//        std::cout << *muonColl << std::endl;
//        std::cout<<"//////////////////"<<std::endl;
//
//        std::cout<<"Tutaj jest l1ObjColl"<<std::endl;
//        std::cout << *l1Objs << std::endl;
//        std::cout<<"//////////////////"<<std::endl;
//
//    }


  }


  L1Obj omtfBX0;
  bool omtfBX0b=false;
  for (const auto & l1mtf2 : l1mtfs){
        if(l1mtf2.bx==0 && l1mtf2.type==L1Obj::OMTF){
          omtfBX0=l1mtf2;
          omtfBX0b=true;
        }
  }

  if(omtfBXm1b==1 && omtfBX0b==0){
    if(omtfBXm1.ptValue()<10) hTimePtbinning->Fill(0.,0.);
    if(omtfBXm1.ptValue()>=10 && omtfBXm1.ptValue()<22) hTimePtbinning->Fill(0.,1.);
    if(omtfBXm1.ptValue()>=22) hTimePtbinning->Fill(0.,2.);

  }
  if(omtfBXm1b==1 && omtfBX0b==1 && omtfBXm1.position==omtfBX0.position && omtfBXm1.iProcessor==omtfBX0.iProcessor){
    if(omtfBXm1.ptValue()<10) hTimePtbinning->Fill(1.,0.);
    if(omtfBXm1.ptValue()>=10 && omtfBXm1.ptValue()<22) hTimePtbinning->Fill(1.,1.);
    if(omtfBXm1.ptValue()>=22) hTimePtbinning->Fill(1.,2.);
  }

  std::bitset<18> l10bit(omtfBXm1.hits);
  std::bitset<18> l10checkDT=(checkDT & l10bit);
  std::bitset<18> l10checkREST=(checkREST & l10bit);
  std::bitset<18> l11bit(omtfBX0.hits);
  std::bitset<18> l11checkDT=(checkDT & l11bit);
  std::bitset<18> l11checkREST=(checkREST & l11bit);
  if(omtfBXm1b==1 && omtfBX0b==0){
    if(omtfBXm1.ptValue()<10) hTimePtbinninghitpattern->Fill(0.,0.);
    if(omtfBXm1.ptValue()>=10 && omtfBXm1.ptValue()<22) hTimePtbinninghitpattern->Fill(0.,1.);
    if(omtfBXm1.ptValue()>=22) hTimePtbinninghitpattern->Fill(0.,2.);

  }
  if(omtfBXm1b==1 && omtfBX0b==1 && omtfBXm1.position==omtfBX0.position && omtfBXm1.iProcessor==omtfBX0.iProcessor){
    if(l10checkDT.count()>0 && l10checkREST.count()==0 && l11checkDT.count()>0 && l11checkREST.count()>0){
      if(omtfBXm1.ptValue()<10) hTimePtbinninghitpattern->Fill(1.,0.);
      if(omtfBXm1.ptValue()>=10 && omtfBXm1.ptValue()<22) hTimePtbinninghitpattern->Fill(1.,1.);
      if(omtfBXm1.ptValue()>=22) hTimePtbinninghitpattern->Fill(1.,2.);
    }
  }


  //
  // coincidence between triggers.
  //
  for (const auto & l1mtf_1 : l1mtfs) {
    for (const auto & l1mtf_2 : l1mtfs) {
      double deltaEta = l1mtf_1.etaValue()-l1mtf_2.etaValue();
      double deltaPhi = reco::deltaPhi( l1mtf_1.phiValue(),  l1mtf_2.phiValue());

      if ( (fabs(deltaEta) > 0.2) || (fabs(deltaPhi) >0.05) ) continue;
      if (l1mtf_1.type==L1Obj::BMTF && l1mtf_2.type ==  L1Obj::OMTF && l1mtf_2.q>=12) {
         hTimeBmtfOmtf->Fill(l1mtf_1.bx,l1mtf_2.bx);

//           std::cout <<" deta: "<<deltaEta<<" dphi: "<<deltaPhi<<std::endl;
//           std::cout << l1mtf_2<<std::endl;
//           std::cout << l1mtf_1<<std::endl;
//         }
      }
      if (l1mtf_1.type==L1Obj::OMTF && l1mtf_1.q>=12  && l1mtf_2.type ==  L1Obj::EMTF) {
         hTimeOmtfEmtf->Fill(l1mtf_1.bx,l1mtf_2.bx);
//         if (l1mtf_1.bx <0 && l1mtf_2.bx==0) { printdeb=true;
//           std::cout <<" deta: "<<deltaEta<<" dphi: "<<deltaPhi<<std::endl;
//           std::cout << l1mtf_1<<std::endl;
//           std::cout << l1mtf_2<<std::endl;
//         }
      }
      if (l1mtf_1.type==L1Obj::OMTF && l1mtf_2.type ==  L1Obj::OMTF_emu) {
         hTimeOmtfOmtf_E->Fill(l1mtf_1.bx,l1mtf_2.bx);
      }
    }
  }

  if (printdeb) {
    std::cout <<"-------- PREFIRE debug, event: "<<*ev<<std::endl; 
    if (muonColl) std::cout << *muonColl << std::endl;
    if (l1Objs)  std::cout << *l1Objs<< std::endl;
  } 


//  if (muonColl){
//    std::cout<<"/////////////////////"<<std::endl;
//    std::cout <<"Tu jest muonColl: "<< *muonColl << std::endl;
//    std::cout<<"/////////////////////"<<std::endl;
//  }
//  if (l1Objs)  {
//    std::cout <<"Tu jest l1ObjColl: "<< *l1Objs << std::endl;
//    std::cout<<"/////////////////////"<<std::endl;
//  }
//  if (ev) { std::cout <<"Tu jest EventObj: "<< *ev << std::endl;
//  std::cout<<"/////////////////////"<<std::endl;
//  }
}
