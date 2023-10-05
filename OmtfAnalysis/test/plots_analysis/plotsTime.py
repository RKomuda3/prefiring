#!/usr/bin/python

import sys
import math
from ROOT import *


def plot1(hIn):
#  h.Print('all')
  h = hIn.Clone()
  integral = max(1.,h.Integral())
  h.Scale(1./integral)
  h.SetMaximum(1.2)
#  h.SetMinimum(1./integral/1.2)
  h.SetMinimum(5.e-4)
  h.SetStats(0);
  t=TLatex(); t.SetNDC(0); t.SetTextSize(0.04); t.SetTextColor(2)
  entries="#Ev: {0:1.3E}".format(integral)
  h.DrawCopy('hist e')
  t.DrawLatex(0.78, 0.8, entries)
  return

def cTimeLayers(canvas):
  c = TCanvas("cTimeLayers","cTimeLayers",1000,600)
  canvas.Add(c)

  h=gROOT.FindObject("hTimeLayers")
  h.SetStats(0)
  h.SetMinimum(0.)
  h.DrawCopy()
  h.SetFillColor(2)
  h.SetFillStyle(3345)
  h.GetXaxis().SetRange(1,6)
  h.DrawCopy('same')

  h.SetFillColor(4)
  h.SetFillStyle(3344)
  h.GetXaxis().SetRange(7,10)
  h.DrawCopy('same')

  h.SetFillColor(8)
  h.SetFillStyle(3354)
  h.GetXaxis().SetRange(11,18)
  h.DrawCopy('same')

#  line = TLine()
#  line.SetLineColor(2)
#  line.DrawLine(5.5,0.5,5.5,vMax)
  c.Update()
  return


def cTimeMtfs(canvas, what1, what2, what3):
  c = TCanvas("cTimeMtfs" + what1+what2+what3, "cTimeMTFs" + what1+what2+what3, 1200, 400)
  canvas.Add(c)
  c.Divide(3)

  pad1 = c.cd(1)
  pad1.SetLogy()
  hTime1 = gROOT.FindObject("hTime" + what1).Clone("hTime" + what1+"Copy")
  plot1(hTime1)

  pad2 = c.cd(2)
  pad2.SetLogy()
  hTime2 = gROOT.FindObject("hTime"+what2)
  plot1(hTime2)
 
  pad3 = c.cd(3)
  pad3.SetLogy()
  hTime3 = gROOT.FindObject("hTime"+what3)
  plot1(hTime3)
 
  c.Update()
  return

def cTimeEffPt(canvas):
  c=TCanvas("cTimeEffPt","cTimeEffPt",1200,400)
  canvas.Add(c)
  c.Divide(3)
  hNull=TH1D("hNull","hNull; p_{T} [GeV];fraction",6,0.,29.9);
  hNull.SetMaximum(1.5)
  hNull.SetMinimum(1.e-5)
  hNull.SetStats(0)

  pad1 = c.cd(1)
  pad1.SetLogy()
  pad1.SetTicky()   
  pad1.SetLeftMargin(0.13)
  hEffPt1 = gROOT.FindObject("hTimeEffPt_BMTF").Clone("hTimeEffPt_BMTF_Copy")
  hNull.SetTitle(hEffPt1.GetTitle())
  hNull.DrawCopy()
  hEffPt1.Draw('same')

  pad2 = c.cd(2)
  pad2.SetLogy()
  pad2.SetTicky()
  pad2.SetLeftMargin(0.13)
  hEffPt2 = gROOT.FindObject("hTimeEffPt_OMTF").Clone("hTimeEffPt_OMTF_Copy")
  hNull.SetTitle(hEffPt2.GetTitle())
  hNull.DrawCopy()
  hEffPt2.Draw('same')

  pad3 = c.cd(3)
  pad3.SetLogy()
  pad3.SetTicky()
  pad3.SetLeftMargin(0.13)
  hEffPt3 = gROOT.FindObject("hTimeEffPt_EMTF").Clone("hTimeEffPt_EMTF_Copy")
  hNull.SetTitle(hEffPt3.GetTitle())
  hNull.DrawCopy()
  hEffPt3.Draw('same')
  c.Update()
  return


def cTimeMtf(canvas, what):
  c = TCanvas("cTimeMtf" + what, "cTimeMTF" + what, 1100, 500)
  canvas.Add(c)
  c.Divide(2)

  pad1 = c.cd(1)
  pad1.SetLogy()
  hTime_A = gROOT.FindObject("hTime" + what+"_A").Clone("hTime" + what+"_A")
  integral = max(1.,hTime_A.Integral())
  hTime_A.Scale(1./integral)
  hTime_A.SetMaximum(1.2)
#  hTime_A.SetMinimum(1./integral/1.2)
  hTime_A.SetMinimum(5.e-5)
  hTime_A.SetStats(0)
  hTime_A.SetTitle(what)
  hTime_A.SetLineColor(8)
  t=TLatex(); t.SetNDC(1); t.SetTextSize(0.04); t.SetTextColor(8)
  entries="#Ev: {0:1.3E}".format(integral)
  hTime_A.DrawCopy('hist ')
  t.DrawLatex(0.68, 0.91, entries)


  hTime_M = gROOT.FindObject("hTime" + what+"_M").Clone("hTime" + what+"_M")
  hTime_M.Scale(1./integral)
  hTime_M.SetLineColor(9)
  hTime_M.SetFillColor(9)
  hTime_M.SetFillStyle(3545)
  hTime_M.DrawCopy('hist same')


  hTime_W = gROOT.FindObject("hTime" + what+"_W").Clone("hTime" + what+"_W")
  hTime_W.Scale(1./integral)
  hTime_W.SetLineColor(46)
  hTime_W.SetFillColor(46)
  hTime_W.SetFillStyle(3554)
  hTime_W.DrawCopy('hist same')

  legend1 = TLegend(0.9, max(5.e-2,1./integral), 2.2, 0.9,"","") #,what+","+entries,"")
  legend1.AddEntry(hTime_A, 'All Q')
  legend1.AddEntry(hTime_M, '#mu, All Q')
  legend1.AddEntry(hTime_W,  'wrong #mu, All Q')
  legend1.Draw()
  canvas.Add(legend1)

 
  pad2 = c.cd(2)
  pad2.SetLogy()

  hTime_Q = gROOT.FindObject("hTime" + what+"_Q").Clone("hTime" + what+"_Q")
  integral = max(1.,hTime_Q.Integral())
  hTime_Q.SetStats(0)
  hTime_Q.Scale(1./integral)
  hTime_Q.SetTitle(what+"_Q")
  hTime_Q.SetLineColor(8)
  hTime_Q.SetMaximum(1.2)
  hTime_Q.SetMinimum(5.e-5)
  hTime_Q.DrawCopy('hist')
  t=TLatex(); t.SetNDC(1); t.SetTextSize(0.04); t.SetTextColor(8)
  entries="#Ev: {0:1.3E}".format(integral)
  t.DrawLatex(0.68, 0.91, entries)

  hTime_QM = gROOT.FindObject("hTime" + what+"_QM").Clone("hTime" + what+"_QM")
  hTime_QM.Scale(1./integral)
  hTime_QM.SetLineColor(9)
  hTime_QM.SetFillColor(9)
  hTime_QM.SetFillStyle(3545)
  hTime_QM.Print("all")
  hTime_QM.DrawCopy('hist same')

  hTime_QW = gROOT.FindObject("hTime" + what+"_QW").Clone("hTime" + what+"_QW")
  hTime_QW.Scale(1./integral)
  hTime_QW.SetLineColor(46)
  hTime_QW.SetFillColor(46)
  hTime_QW.SetFillStyle(3554)
  hTime_QW.DrawCopy('hist e same')

  legend = TLegend(0.9, max(5.e-2,1./integral), 2.2, 0.9,"","") #,what+","+entries,"")
  legend.AddEntry(hTime_Q, 'Q #geq 12')
  legend.AddEntry(hTime_QM, '#mu, Q #geq 12')
  legend.AddEntry(hTime_QW, 'wrong #mu, Q #geq 12')
  legend.Draw()
  canvas.Add(legend)

  c.Update()
  return


def cTimeMtfsCorr(canvas,what=''):
  c = TCanvas("cTimeMtfsCorr"+what,"cTimeMTFsCorr"+what,1200,400)
  canvas.Add(c)
  c.Divide(3)

  pad1 = c.cd(1)
  hTimeBmtfOmtf = gROOT.FindObject("hTimeBmtfOmtf"+what)
  hTimeBmtfOmtf.SetMinimum(0.5)
  hTimeBmtfOmtf.DrawCopy('col text')

  pad2 = c.cd(2)
  hTimeOmtfEmtf= gROOT.FindObject("hTimeOmtfEmtf"+what)
  hTimeOmtfEmtf.SetMinimum(0.5)
  hTimeOmtfEmtf.DrawCopy('col text')

  pad3 = c.cd(3)
  hTimeOmtfOmtf_E= gROOT.FindObject("hTimeOmtfOmtf_E"+what)
  hTimeOmtfOmtf_E.SetMinimum(0.5)
  hTimeOmtfOmtf_E.DrawCopy('col text')
  c.Update()
  return

def cTimeTrackPt(canvas) :
  c =  TCanvas("cTimeTrackPt","cTimeTrackPt",800,800) 
  canvas.Add(c)
  c.Divide(2,2)
#  c.cd(1); hTimeOmtfTrackDPhiT.DrawCopy('box')
  c.cd(1);  hTimeOmtfTrackDRM.DrawCopy()
  c.cd(2); hTimeOmtfTrackDPhiM.DrawCopy('box')
#  c.cd(3); hTimeOmtfTrackDEtaT.DrawCopy('box')
  c.cd(4); hTimeOmtfTrackDEtaM.DrawCopy('box')
  return

def cTimeTrackBX(canvas) :
  c =  TCanvas("cTimeTrackBX","cTimeTrackBX",800,800) 
  canvas.Add(c)
  c.Divide(2,2)

  c.cd(1); hTimeOmtfTrackBXT.DrawCopy('box')
  c.cd(2); hTimeOmtfTrackBXM.DrawCopy('box')
  c.cd(3); hTimeOmtfDrTrackMuon.DrawCopy('box')
  c.cd(4)
  hTimeOmtfTrackBX1.Divide(hTimeOmtfTrackBX1,hTimeOmtfTrackBX0);
  hTimeOmtfTrackBX1.DrawCopy()
  return

def cTimeDeltaR(canvas) :
  c = TCanvas("cTimeDeltR","cTimeDeltR",800,800)
  canvas.Add(c)
  hTimeDeltaR_Q=gROOT.FindObject("hTimeDeltaR_Q")
  c.SetLogy()
  hTimeDeltaR_Q.SetStats(0)
  hTimeDeltaR_Q.SetLineColor(4)
  hTimeDeltaR_Q.SetMinimum(0.9)
  hTimeDeltaR_Q.DrawCopy()
  hTimeDeltaR_QW=gROOT.FindObject("hTimeDeltaR_QW")
  hTimeDeltaR_QW.SetLineColor(2)
  hTimeDeltaR_QW.DrawCopy('same')
  for i in range(3) :
    hTimeDeltaR_Q_B=gROOT.FindObject("hTimeDeltaR_Q_B"+str(i+1))
    hTimeDeltaR_Q_B.SetLineColor(7)
    hTimeDeltaR_Q_B.DrawCopy('same')
    hTimeDeltaR_QW_B=gROOT.FindObject("hTimeDeltaR_QW_B"+str(i+1))
    hTimeDeltaR_QW_B.SetLineColor(6)
    hTimeDeltaR_QW_B.DrawCopy('same')
  legend = TLegend(2.0, 0.3*hTimeDeltaR_Q.GetMaximum(), 3.8, hTimeDeltaR_Q.GetMaximum(),"","")
  legend.AddEntry(hTimeDeltaR_Q, '#mu Q #geq 12')
  legend.AddEntry(hTimeDeltaR_QW, 'wrong #mu, Q #geq 12')
  legend.AddEntry(hTimeDeltaR_Q_B, '#mu, Q #geq 12, wrong BX')
  legend.AddEntry(hTimeDeltaR_QW_B, 'wrong #mu, Q #geq 12, wrong BX')
  legend.Draw()
  canvas.Add(legend)

  c.Update()
  return

def cTimeEta(canvas) :
  c=TCanvas("cTimeEta","cTimeEta",600,400)
  canvas.Add(c)
  c.SetLogy()
  c.SetTicky()
  frame=c.DrawFrame(-2.4, 1.e-6, 2.4, 1.02);
  frame.SetYTitle("fraction")
  frame.SetXTitle("\eta")
  h0 = gROOT.FindObject("hTimeEta_Pt0").Clone("hTimeEta_Pt0")
  h0.SetLineColor(2)
  h0.SetMarkerStyle(24)
  h0.SetMarkerColor(2)
  h10 = gROOT.FindObject("hTimeEta_Pt10").Clone("hTimeEta_Pt10")
  h10.SetLineColor(3)
  h10.SetMarkerStyle(25)
  h10.SetMarkerColor(3)
  h22 = gROOT.FindObject("hTimeEta_Pt22").Clone("hTimeEta_Pt22")
  h22.SetLineColor(4)
  h22.SetMarkerStyle(26)
  h22.SetMarkerColor(4)
  h0.Draw('same')
  h10.Draw('same')
  h22.Draw('same')
  c.Update()
  return 
  
  
####################################moja praca  
  
def chTimePrefireEta(canvas):
    c=TCanvas("chTimePrefireEta","chTimePrefireEta",700,400)
    canvas.Add(c)
    c.SetGridy(100)
    c.SetGridx(100)
    hhTimePrefireEta = gROOT.FindObject("hTimePrefireEta")
    hhTimePrefireEta.SetLineColor(1)
    hhTimePrefireEta.SetTitle("Prefiring \, for \, bx_{-1} \, p_{T} \geq 10 \, GeV; \eta [-]; bx-1/(bx-0 or bx-1)")
    hhTimePrefireEta.Draw()
    c.Update()
    graph = hhTimePrefireEta.GetPaintedGraph()
    graph.GetYaxis().SetNdivisions(808)
    graph.GetXaxis().SetNdivisions(808)
    graph.SetStats(1)
    c.Update()
    return

def chTimePrefireEta1(canvas):
    c=TCanvas("chTimePrefireEta1","chTimePrefireEta1",700,400)
    canvas.Add(c)
    c.SetGridy(100)
    c.SetGridx(100)
    hhTimePrefireEta1 = gROOT.FindObject("hTimePrefireEta1")
    hhTimePrefireEta1.SetLineColor(1)
    hhTimePrefireEta1.SetTitle("Prefiring \, with \, OMTF \, veto \, for  \, bx_{-1} \, p_{T} \geq \, 10 GeV; \eta [-]; bx-1/(bx-0 or bx-1)")
    hhTimePrefireEta1.Draw()
    c.Update()
    graph = hhTimePrefireEta1.GetPaintedGraph()
    graph.GetYaxis().SetNdivisions(808)
    graph.GetXaxis().SetNdivisions(808)
    graph.SetStats(1)
    c.Update()
    return


def cDeltaphieta(canvas):
  c = TCanvas("cDeltaphieta","cDeltaphieta",1300,400)
  canvas.Add(c)
  c.Divide(3)
  pad1 = c.cd(1)
  hDeltaphieta010 = gROOT.FindObject("hDeltaphieta010")
  hDeltaphieta010.SetTitle("OMTF \ between \, bx_{0} \, and \, bx_{-1} \ p_T<10 GeV; \Delta \eta; \Delta \phi ")
  hDeltaphieta010.SetStats(0)
  hDeltaphieta010.Draw('BOX,TEXT')
  
  
  pad2 = c.cd(2)
  hDeltaphieta1022 = gROOT.FindObject("hDeltaphieta1022")
  hDeltaphieta1022.SetStats(0)
  hDeltaphieta1022.SetTitle("OMTF \ between \, bx_{0} \, and \, bx_{-1} \ 10 GeV \leq p_T<22 GeV; \Delta \eta; \Delta \phi ")
  hDeltaphieta1022.Draw('BOX,TEXT')
  
  pad3 = c.cd(3)
  hDeltaphieta22inf = gROOT.FindObject("hDeltaphieta22inf")
  hDeltaphieta22inf.SetStats(0)
  hDeltaphieta22inf.SetTitle("OMTF \ between \, bx_{0} \, and \, bx_{-1} \ p_T \geq 22 GeV ; \Delta \eta; \Delta \phi ")
  hDeltaphieta22inf.Draw('BOX,TEXT')
  
  c.Update()
  return 
  
def cDeltaphieta025(canvas):
  c = TCanvas("cDeltaphieta025","cDeltaphieta025",1300,400)
  canvas.Add(c)
  c.Divide(3)
  pad1 = c.cd(1)
  hDeltaphieta010_025 = gROOT.FindObject("hDeltaphieta010_025")
  hDeltaphieta010_025.SetTitle("OMTF \ between \, bx_{0} \, and \, bx_{-1} \ p_T<10 GeV; \Delta \eta; \Delta \phi ")
  hDeltaphieta010_025.SetStats(0)

  hDeltaphieta010_025.DrawCopy('BOX,TEXT')
  
  pad2 = c.cd(2)
  hDeltaphieta1022_025 = gROOT.FindObject("hDeltaphieta1022_025")
  hDeltaphieta1022_025.SetTitle("OMTF \ between \, bx_{0} \, and \, bx_{-1} \ 10 GeV \leq p_T<22 GeV; \Delta \eta; \Delta \phi ")
  hDeltaphieta1022_025.SetStats(0)
  hDeltaphieta1022_025.Draw('BOX,TEXT')

  
  pad3 = c.cd(3)
  hDeltaphieta22inf_025 = gROOT.FindObject("hDeltaphieta22inf_025")
  hDeltaphieta22inf_025.SetTitle("OMTF \ between \, bx_{0} \, and \, bx_{-1} \ p_T \geq 22 GeV; \Delta \eta; \Delta \phi ")
  hDeltaphieta22inf_025.SetStats(0)
  hDeltaphieta22inf_025.DrawCopy('BOX,TEXT')

  
  c.Update()
  return  
  
  
def  chQualitybx0(canvas): 
    
    pad = TPad("pad", "Pad", 0, 0, 1, 1)
    pad.SetLeftMargin(0.2)     
    c=TCanvas("chQualitybx0","chQualitybx0",700,400)
    pad.Draw()
    pad.cd()

    canvas.Add(c)
    hQualitybx0=gROOT.FindObject("hQualitybx0")
    hQualitybx0.GetYaxis().SetBinLabel(1,"p_{t}<10 GeV")
    hQualitybx0.GetYaxis().SetBinLabel(2,"10 GeV $\\leq$ p_{t} $<$22 GeV")
    hQualitybx0.GetYaxis().SetBinLabel(3,"p_{t}$\\geq$22 GeV")
    hQualitybx0.Draw('BOX,TEXT')
    hQualitybx0.SetStats(0)
    
    pad.Update()

    c.Update()
    return

def  chTimePt22(canvas):
    c=TCanvas("chTimePt22","chTimePt22",700,400)
    c.SetLogy()
    canvas.Add(c)
    hTimePt22=gROOT.FindObject("hTimePt22")
    hTimePt22.Draw()
    hTimePt22.Scale(1/hTimePt22.Integral(),'WIDTH')
    
    c.Update()
    return
def  chTimePtbinning(canvas):
    
    pad = TPad("pad", "Pad", 0, 0, 1, 1)
    pad.SetLeftMargin(0.2) 
    c=TCanvas("chTimePtbinning","chTimePtbinning",800,500)
    
    pad.Draw()
    pad.cd()
    
    canvas.Add(c)
    hTimePtbinning=gROOT.FindObject("hTimePtbinning")
    hTimePtbinning.GetYaxis().SetBinLabel(1,"p_{t}<10 GeV")
    hTimePtbinning.GetYaxis().SetBinLabel(2,"10 GeV $\\leq$ p_{t} $<$22 GeV")
    hTimePtbinning.GetYaxis().SetBinLabel(3,"p_{t}$\\geq$22 GeV")
    hTimePtbinning.GetXaxis().SetBinLabel(2,"bx_{-1} and bx_{0}")
    hTimePtbinning.GetXaxis().SetBinLabel(1,"bx_{-1}")
    hTimePtbinning.SetTitle("OMTF")
    hTimePtbinning.GetXaxis().SetTitle("")
    hTimePtbinning.SetStats(0)
    hTimePtbinning.Draw('BOX,TEXT')
    
    c.Update()
    c.SetLeftMargin(0.8)
    pad.Update()
    c.Update()
    return

def cTimeEffPt1(canvas):
  c=TCanvas("cTimeEffPt1","cTimeEffPt1",1200,400)
  canvas.Add(c)
  c.Divide(3)
  hNull=TH1D("hNull","hNull; p_{T} [GeV];fraction",6,0.,29.9);
  hNull.SetMaximum(1.5)
  hNull.SetMinimum(1.e-5)
  hNull.SetStats(0)

  pad1 = c.cd(1)
  pad1.SetLogy()
  pad1.SetTicky()
  pad1.SetLeftMargin(0.13)
  hEffPt1 = gROOT.FindObject("hTimeEffPt1_BMTF").Clone("hTimeEffPt1_BMTF_Copy")
  hNull.SetTitle(hEffPt1.GetTitle())
  hNull.DrawCopy()
  hEffPt1.Draw('same')

  pad2 = c.cd(2)
  pad2.SetLogy()
  pad2.SetTicky()
  pad2.SetLeftMargin(0.13)
  hEffPt2 = gROOT.FindObject("hTimeEffPt1_OMTF").Clone("hTimeEffPt1_OMTF_Copy")
  hNull.SetTitle(hEffPt2.GetTitle())
  hNull.DrawCopy()
  hEffPt2.Draw('same')

  pad3 = c.cd(3)
  pad3.SetLogy()
  pad3.SetTicky()
  pad3.SetLeftMargin(0.13)
  hEffPt3 = gROOT.FindObject("hTimeEffPt1_EMTF").Clone("hTimeEffPt1_EMTF_Copy")
  hNull.SetTitle(hEffPt3.GetTitle())
  hNull.DrawCopy()
  hEffPt3.Draw('same')
  c.Update()
  return

def plotAll(canvas) :
#  cTimeMtfsCorr(canvas)
#  cTimeMtfs(canvas,'Bmtf_A','Omtf_A','Emtf_A')
#  cTimeMtfs(canvas,'Bmtf_Q','Omtf_Q','Emtf_Q')
#  cTimeMtfs(canvas,'Bmtf_M','Omtf_M','Emtf_M')
#  cTimeMtfs(canvas,'Bmtf_QM','Omtf_QM','Emtf_QM')
#  cTimeMtfs(canvas,'Omtf_Q','Omtf_QM','Omtf_emu_QM')
#  cTimeMtfs(canvas,'Bmtf_W','Omtf_W','Emtf_W')
#  cTimeMtfs(canvas,'Bmtf_QW','Omtf_QW','Emtf_QW')
#  cTimeMtf(canvas,'Bmtf')
#  cTimeMtf(canvas,'Omtf')
#  cTimeMtf(canvas,'Omtf_emu')
#  cTimeMtf(canvas,'Emtf')
#  cTimeTrackPt(canvas)
#  cTimeTrackBX(canvas)
#  cTimeDeltaR(canvas)
#  cTimeLayers(canvas)
   cTimeEffPt(canvas)
#  cTimeEta(canvas)
   chTimePrefireEta(canvas)
   chTimePrefireEta1(canvas)
   cDeltaphieta(canvas)
   cDeltaphieta025(canvas)
   chQualitybx0(canvas)
#  chTimePt22(canvas)
   chTimePtbinning(canvas)
   cTimeEffPt1(canvas)
   return


