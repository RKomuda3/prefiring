#include "UserCode/OmtfAnalysis/interface/AnaEff.h"
#include "TProfile.h"
#include "TObjArray.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TTree.h"
#include "TFile.h"

#include "UserCode/OmtfDataFormats/interface/EventObj.h"
#include "UserCode/OmtfDataFormats/interface/MuonObj.h"
#include "UserCode/OmtfDataFormats/interface/L1Obj.h"
#include "UserCode/OmtfDataFormats/interface/L1ObjColl.h"
#include "UserCode/OmtfAnalysis/interface/Utilities.h"
#include "DataFormats/Math/interface/deltaR.h"

#include <cmath>
#include <vector>
#include <ostream>
#include <bitset>
#include <iostream>

namespace {

  const static unsigned int nPtCuts= 4;
  const double ptCuts[nPtCuts] ={ 0., 10., 16., 25. }; 

  std::string reg[3]={"_Bar","_Ove","_End"};

  std::map< std::string, TH1D* > hm;
  TH1D  *hEffEtaOMTFn, *hEffEtaOMTFn_D, *hEffEtaOMTFp, *hEffEtaOMTFp_D;
  TH1D  *hEffEtaAll, *hEffEtaAll_D;
  TH1D  *hEffDeltaR, *hEffDeltaEta, *hEffDeltaPhi;
  TH2D  *hEffEtaMuVsEtauGMT,*hEffEtaMuVsEtaOMTF; 
  TH1D  *hEffRunAver, *hEffRunAverWeighted;

  struct BestL1Obj : public L1Obj {
    BestL1Obj() : deltaR(9999.) {}
    BestL1Obj(const L1Obj & l1, const  MuonObj *muon) : L1Obj(l1) , deltaR(9999.) {
      if (l1.isValid() && muon) deltaR = reco::deltaR( l1.etaValue(),l1.phiValue(), muon->l1Eta, muon->l1Phi);
//    if (l1.isValid() && muon) deltaR = reco::deltaR( l1.etaValue(),l1.phiValue(), muon->eta(), muon->phi());
    }
    static constexpr double theDeltaRCut = 0.4;
    bool fired(double ptCut=0., int qMin=12, double matchingDeltaR=BestL1Obj::theDeltaRCut) const {
      double epsilon = 1.e-5;
      return (q>=qMin) && (std::fabs(ptValue())+epsilon >= ptCut) && (deltaR <= matchingDeltaR); 
    }
    bool isBetter(const BestL1Obj & best) const {
      if (type != best.type || q<12) return false; 
      //if (deltaR-best.deltaR < -1.e-6) return true;
      if (deltaR > theDeltaRCut) return false;
      if ( ptValue() > best.ptValue() ) return true;
      return false; 
    }
    double deltaR;
    friend std::ostream & operator<< (std::ostream &out, const BestL1Obj &o) {
      out << (L1Obj)o << " --> BestL1Obj : "<<o.deltaR<<" fired: "<<o.fired();
      return out;
    }
  };

}




void AnaEff::resume(TObjArray& histos)
{
  TGraphErrors * hGraphRun = new TGraphErrors();
  hGraphRun->SetName("hGraphEffRun");
  histos.Add(hGraphRun);
  std::vector<unsigned int> runs = theRunMap.runs();
  hGraphRun->Set(runs.size());
  double epsil = 1.e-6;
  for (unsigned int iPoint = 0; iPoint < runs.size(); iPoint++) {
    unsigned int run = runs[iPoint];
    RunEffMap::EffAndErr effAndErr = theRunMap.effAndErr(run);
    std::cout <<" RUN: "<<run <<" eff: "<< effAndErr 
              <<" stat: "<<theRunMap.stat(run).first<<"/"<<theRunMap.stat(run).second<<std::endl;
    hGraphRun->SetPoint(iPoint, run, effAndErr.eff());
    hGraphRun->SetPointError(iPoint, 0., effAndErr.effErr());
    hEffRunAver->Fill(effAndErr.eff()-epsil);
    hEffRunAverWeighted->Fill(effAndErr.eff()-epsil, theRunMap.stat(run).first);
  }
  std::cout << "Eff run summary: " << theRunMap.effAndErr() << std::endl;
   
}


void AnaEff::init(TObjArray& histos)
{
  int nBins = 60;
  double omin = 0.75;
  double omax = 1.35;
  hEffEtaOMTFn   = new TH1D("hEffEtaOMTFn",  "hEffEtaOMTFn",   nBins, -omax, -omin); histos.Add(hEffEtaOMTFn);
  hEffEtaOMTFn_D = new TH1D("hEffEtaOMTFn_D","hEffEtaOMTFn_D", nBins, -omax, -omin); histos.Add(hEffEtaOMTFn_D);
  hEffEtaOMTFp   = new TH1D("hEffEtaOMTFp",  "hEffEtaOMTFp",   nBins,  omin, omax); histos.Add(hEffEtaOMTFp);
  hEffEtaOMTFp_D = new TH1D("hEffEtaOMTFp_D","hEffEtaOMTFp_D", nBins,  omin, omax); histos.Add(hEffEtaOMTFp_D);

  hEffRunAver = new TH1D("hEffRunAver","hEffRunAver", 50, 0., 1.); histos.Add(hEffRunAver); hEffRunAver->Sumw2();
  hEffRunAverWeighted = new TH1D("hEffRunAverWeighted","hEffRunAverWeighted", 50, 0., 1.); histos.Add(hEffRunAverWeighted); hEffRunAverWeighted->Sumw2();

  hEffEtaAll     =  new TH1D("hEffEtaAll",    "hEffEtaAll",   96,   -2.4, 2.4); histos.Add(hEffEtaAll); 
  hEffEtaAll_D   =  new TH1D("hEffEtaAll_D",  "hEffEtaAll_D", 96,   -2.4, 2.4); histos.Add(hEffEtaAll_D); 

  hEffDeltaR     = new TH1D( "hEffDeltaR", "hEffDeltaR", 60, 0., 0.6);  histos.Add(hEffDeltaR);
  hEffDeltaEta     = new TH1D( "hEffDeltaEta", "hEffDeltaEta", 60, -0.3, 0.3);  histos.Add(hEffDeltaEta);
  hEffDeltaPhi     = new TH1D( "hEffDeltaPhi", "hEffDeltaPhi", 60, -0.3, 0.3);  histos.Add(hEffDeltaPhi);

  hEffEtaMuVsEtauGMT = new TH2D("hEffEtaMuVsEtauGMT","hEffEtaMuVsEtauGMT",   80, 0.3,1.8, 80, 0.3,1.8); histos.Add(hEffEtaMuVsEtauGMT);
  hEffEtaMuVsEtaOMTF = new TH2D("hEffEtaMuVsEtaOMTF","hEffEtaMuVsEtaOMTF",   80, 0.3,1.8, 80, 0.3,1.8); histos.Add(hEffEtaMuVsEtaOMTF);
/*
  hEfficMuPt_D = new TH1D("hEfficMuPt_D","hEfficMuPt_D", L1PtScale::nPtBins, L1PtScale::ptBins); histos.Add(hEfficMuPt_D);
  hEfficRpcNoCut_N = new TH1D("hEfficRpcNoCut_N","hEfficRpcNoCut_N", L1PtScale::nPtBins, L1PtScale::ptBins);  histos.Add(hEfficRpcNoCut_N);
  hEfficRpcPtCut_N = new TH1D("hEfficRpcPtCut_N","hEfficRpcPtCut_N", L1PtScale::nPtBins, L1PtScale::ptBins);  histos.Add(hEfficRpcPtCut_N);
*/

  std::string  base("hEff");
  const unsigned int nOpt = 5;
//  std::string opt[nOpt]={"_uGmtPtCut","_BmtfPtCut","_OmtfPtCut","_EmtfPtCut"};
  std::string opt[nOpt]={"_uGmt","_Omtf","_Bmtf","_Emtf","_OmtfQ4"};
  for (unsigned int ir=0; ir<3; ++ir) {
    std::string name=base+"_PtCutDenom"+reg[ir];
    TH1D *h= new TH1D(name.c_str(),name.c_str(), L1PtScale::nPtBins, L1PtScale::ptBins); h->Sumw2(); 
    histos.Add(h); theHistoMap[name]=h;
    //for (unsigned int iopt=0; iopt< nOpt; ++iopt) {
    for (unsigned int iopt=0; iopt< 2; ++iopt) {
    for (unsigned int icut=0; icut<nPtCuts; ++icut) {
      std::ostringstream str;
      str << base << opt[iopt]<<"PtCut"<< ptCuts[icut]<<reg[ir];
      TH1D *h= new TH1D(str.str().c_str(),str.str().c_str(), L1PtScale::nPtBins, L1PtScale::ptBins); h->Sumw2(); 
      h->SetXTitle("muon p_{T} [GeV/c]  "); h->SetYTitle("events");
      histos.Add(h); theHistoMap[str.str()]=h;
    }
    }
  }

//for (unsigned int icut=2; icut<nPtCuts; ++icut) {
  {
    unsigned int icut=2;
    std::ostringstream strEtaDenom; strEtaDenom << base << "_EtaDenom"<< ptCuts[icut];
//  TH1D *hD = new TH1D(strEtaDenom.str().c_str(),strEtaDenom.str().c_str(), L1RpcEtaScale::nEtaBins, L1RpcEtaScale::etaBins); hD->Sumw2();
    TH1D *hD = new TH1D(strEtaDenom.str().c_str(),strEtaDenom.str().c_str(), 96, -2.4, 2.4); hD->Sumw2();
    hD->SetXTitle("muon pseudorapidity"); hD->SetYTitle("events"); histos.Add(hD); theHistoMap[strEtaDenom.str()]=hD;
    for (unsigned int iopt=0; iopt <nOpt; ++iopt) {
      std::ostringstream strEtaCut  ; strEtaCut << base << "_EtaCut"<< ptCuts[icut]<<opt[iopt];
//    TH1D *hN = new TH1D(strEtaCut.str().c_str(),strEtaCut.str().c_str(), L1RpcEtaScale::nEtaBins, L1RpcEtaScale::etaBins); hN->Sumw2();
      TH1D *hN = new TH1D(strEtaCut.str().c_str(),strEtaCut.str().c_str(), 96, -2.4, 2.4); hN->Sumw2();
      hN->SetXTitle("muon pseudorapidity"); hN->SetYTitle("events"); histos.Add(hN); theHistoMap[strEtaCut.str()]=hN;
    }
  }

}

void AnaEff::run(  const EventObj* event, const MuonObj* muon, const L1ObjColl *l1Coll)
{
  if (!muon) return;
  double etaMu = muon->eta();
  double ptMu  = muon->pt();  
//std::cout <<" MUON: " << *muon << std::endl;
  if (ptMu < 0.1) return;
  if (!l1Coll || l1Coll->getL1Objs().size()==0) return;
//  if (muon && l1Coll) std::cout <<*event << std::endl << *muon<< std::endl<<*l1Coll<<std::endl<<std::endl;

  //
  // best (closest)  MTF L1Obj to muon
  //
  L1Obj::TYPE typeuGMT = L1Obj::uGMT;
  L1Obj::TYPE typeBMTF = L1Obj::BMTF;
  L1Obj::TYPE typeOMTF = L1Obj::OMTF;
  L1Obj::TYPE typeEMTF = L1Obj::EMTF;
  //L1Obj::TYPE typeOMTF = L1Obj::OMTF_emu;
  BestL1Obj bestuGMT; bestuGMT.type=typeuGMT;
  BestL1Obj bestBMTF; bestBMTF.type=typeBMTF;
  BestL1Obj bestOMTF; bestOMTF.type=typeOMTF;
  BestL1Obj bestEMTF; bestEMTF.type=typeEMTF;
//  std::vector<L1Obj> l1s = l1Coll->selectByBx(0,0);
  std::vector<L1Obj> l1s0 = l1Coll->selectByBx(-1,-1);
  std::vector<L1Obj> l1s1 = l1Coll->selectByBx(0,0);
  std::vector<L1Obj> l1s ;

//  std::cout<<"PRZERWAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"<<std::endl;
//  bool wasprefire=false;
//  std::cout<<"HITYYYYYYYYYYYYYYYYYYYYYYYYYYY"<<std::endl;

  std::bitset<18>   checkDT(std::string("000000000000111111"));
  std::bitset<18> checkREST(std::string("111111111111000000"));
  for (const auto & l10 : l1s0) {
    bool checkprefire=false;
    std::bitset<18> l10bit(l10.hits);
    std::bitset<18> l10checkDT=(checkDT & l10bit);
    std::bitset<18> l10checkREST=(checkREST & l10bit);

    for (const auto & l11 : l1s1) {
       std::bitset<18> l11bit(l11.hits);
       std::bitset<18> l11checkDT=(checkDT & l11bit);
       std::bitset<18> l11checkREST=(checkREST & l11bit);
       if(l10.q >= 12 && l10.type==L1Obj::OMTF && l11.type==L1Obj::OMTF && l10.position==l11.position && l10.iProcessor==l11.iProcessor  && (std::abs(reco::deltaPhi(l10.phiValue(),l11.phiValue()))<0.09)){
//         std::cout<<"//////////////////////////////////////////////////////////////////////////"<<std::endl;
//         std::cout<<l10<<std::endl;
//         std::cout<<"??????????????????????????????????????????????????????????????????????????"<<std::endl;
//         std::cout<<l11<<std::endl;
//         std::cout<<"??????????????????????????????????????????????????????????????????????????"<<std::endl;

         if(l10checkDT.count()>0 && l10checkREST.count()==0){
             if(l11checkDT.count()>0 && l11checkREST.count()>0) std::cout<<"";
         }
           checkprefire=true;
//         wasprefire=true;
       }
    }
//    std::cout<<"H/////////////////////////////"<<std::endl;
//    std::cout<<l10.hits<<std::endl;
//    std::cout<<"K/////////////////////////////"<<std::endl;
    if(checkprefire==true) l1s.push_back(l10);

  }
//  std::cout<<"KONIECHITOOOOOOOOOOOOOOOOOOOOW"<<std::endl;
//    if(wasprefire){
//      std::cout<<"POCZAAAAAAAATEEEEEEKKKKKK"<<std::endl;
//      if (muon && l1Coll ) std::cout <<*event << std::endl << *muon<< std::endl<<*l1Coll<<std::endl<<std::endl;
//      std::cout<<"PRZERWAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"<<std::endl;
//      for(const auto &ls :l1ss){
//        std::cout<<ls<<std::endl;
//
//      }
//      std::cout<<"Koooooooooooooooooooooooooooooooooooooooooooooooooooniec"<<std::endl;
//    }



  for (const auto & l1 : l1s) {

    BestL1Obj cand(l1,muon);

    if (cand.isBetter(bestuGMT)) bestuGMT=cand;
    if (cand.isBetter(bestBMTF)) bestBMTF=cand;
    if (cand.isBetter(bestOMTF)) bestOMTF=cand;
    if (cand.isBetter(bestEMTF)) bestEMTF=cand;
  }

  //
  // controls for bestOMTF
  //
  if (bestOMTF.isValid() && fabs(etaMu) > 0.9 && fabs(etaMu) < 1.2) {
    hEffDeltaR->Fill(bestOMTF.deltaR); 
    hEffDeltaPhi->Fill(reco::deltaPhi(bestOMTF.phiValue(),muon->l1Phi)); 
    hEffDeltaEta->Fill(bestOMTF.etaValue()-muon->l1Eta); 
  } 

  //
  // OMTF efficiency history 
  //
  if ( ptMu > 7. && fabs(etaMu) < 1.24 && fabs(etaMu) > 0.83 ) theRunMap.addEvent(event->run, bestOMTF.fired()); 

  //
  // efficiency with ptCut for each region
  //
  unsigned int iregion;
  if ( fabs(etaMu) < 0.83) iregion = 0;
  else if ( fabs(etaMu) < 1.24) iregion = 1;
  else iregion = 2;

  theHistoMap["hEff_PtCutDenom"+reg[iregion]]->Fill(ptMu);
  for (unsigned int icut=0; icut < nPtCuts; icut++) { 
    double threshold = ptCuts[icut];
    if ( bestOMTF.fired(threshold) ) {
       std::ostringstream strPt;  strPt  << "hEff_OmtfPtCut"<<  ptCuts[icut]<<reg[iregion];
       theHistoMap[strPt.str()]->Fill(ptMu);
    } 
    if ( bestuGMT.fired(threshold) ) {
       std::ostringstream strPt;  strPt  << "hEff_uGmtPtCut"<<  ptCuts[icut]<<reg[iregion];
       theHistoMap[strPt.str()]->Fill(ptMu);
    } 
  }

  //
  // efficiency vs Eta
  //
  if (ptMu > 7.) {
    TH1D* h_N = (etaMu > 0 ) ? hEffEtaOMTFp : hEffEtaOMTFn; 
    TH1D* h_D = (etaMu > 0 ) ? hEffEtaOMTFp_D : hEffEtaOMTFn_D; 
    h_D->Fill(etaMu); 
    if (bestOMTF.fired()) h_N->Fill(etaMu);
    hEffEtaAll_D->Fill(etaMu);
    if (bestuGMT.fired()) hEffEtaAll->Fill(etaMu);
  }

  //
  // other histos
  //
  if ( fabs(etaMu) > 0.3 && fabs(etaMu) < 1.74 ) {
    if (bestuGMT.fired()) hEffEtaMuVsEtauGMT->Fill(fabs(muon->l1Eta), fabs(bestuGMT.etaValue()));
    if (bestOMTF.pt>0)    hEffEtaMuVsEtaOMTF->Fill(fabs(muon->l1Eta), fabs(bestOMTF.etaValue()));
  }

  double threshold = ptCuts[2];
  double ptMuPlat = 1.5*threshold; 
  if (ptMu >= ptMuPlat) {
    std::ostringstream strEtaDenom;  strEtaDenom  << "hEff_EtaDenom"<<  threshold;
    theHistoMap[strEtaDenom.str()]->Fill(muon->eta()); 

    if ( bestuGMT.fired(threshold) ) {
       std::ostringstream strEtaCut;  strEtaCut  << "hEff_EtaCut"<<  threshold <<"_uGmt";
       theHistoMap[strEtaCut.str()]->Fill(muon->eta()); 
    }
    if ( bestBMTF.fired(threshold) ) {
       std::ostringstream strEtaCut;  strEtaCut  << "hEff_EtaCut"<<  threshold <<"_Bmtf";
       theHistoMap[strEtaCut.str()]->Fill(muon->eta()); 
    }
    if ( bestOMTF.fired(threshold) ) {
       std::ostringstream strEtaCut;  strEtaCut  << "hEff_EtaCut"<<  threshold <<"_Omtf";
       theHistoMap[strEtaCut.str()]->Fill(muon->eta()); 
    }
    if ( bestEMTF.fired(threshold) ) {
       std::ostringstream strEtaCut;  strEtaCut  << "hEff_EtaCut"<<  threshold <<"_Emtf";
       theHistoMap[strEtaCut.str()]->Fill(muon->eta()); 
    }
//    if ( fabs(muon->eta()) < 0.6 && !bestuGMT.fired(threshold) && bestBMTF.fired(threshold) ) {
//    std::cout <<" Muon: " << *muon << std::endl;
//    std::cout <<*l1Coll << std::endl;
//    std::cout <<" BestuGMT: " << bestuGMT << std::endl;
//    std::cout <<" BestBMTF: " << bestBMTF<<" fired: "<< bestBMTF.fired(threshold)<< std::endl;
//    }

  }
}	
