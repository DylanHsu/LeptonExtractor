#include "TROOT.h"
#include "TPaveText.h"
#include "TStyle.h"
#include "TH1D.h"
#include "TF1.h"
#include "TTree.h"
#include "TMath.h"
#include "TFile.h"
#include "TLine.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TIterator.h"
#include "TComplex.h"
#include "TVirtualFFT.h"

#include <vector>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>

#ifndef __CINT__
#include "RooGlobalFunc.h"
#include "RooMsgService.h"
#endif
#include "RooAbsPdf.h"
#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooAddPdf.h"
#include "RooArgList.h"
#include "RooArgSet.h"
#include "RooFFTConvPdf.h"
#include "RooFitResult.h"
#include "RooPlot.h"
#include "fitterShape.h"

using namespace RooFit;
const int fitMassLo=80;
const int fitMassHi=100;
int range=fitMassHi-fitMassLo;
void findTransferFunction(
  string histName="pass_ptBin0_etaBin0",
  string genFileName="/home/dhsu/CMSSW_8_0_26_patch1/src/LeptonExtractor/transFuncStudy/DYJetsToLL_BaselineToMedium_GenMass_PtBinnedPlusInclusiveNLO_electronTnP_binnedHistos.root",
  string recoFileName="/home/dhsu/CMSSW_8_0_26_patch1/src/LeptonExtractor/transFuncStudy/DYJetsToLL_BaselineToMedium_ScaleSmearCorrections_PtBinnedPlusInclusiveNLO_electronTnP_binnedHistos.root"
) {
  TFile *genFile=TFile::Open(genFileName.c_str(),"READ"); assert(genFile && genFile->IsOpen());
  TH1D *genHist=(TH1D*)genFile->Get(histName.c_str())->Clone("genHist"); assert(genHist); genHist->SetDirectory(0); genFile->Close();
  TFile *recoFile=TFile::Open(recoFileName.c_str(),"READ"); assert(recoFile && recoFile->IsOpen());
  TH1D *recoHist=(TH1D*)recoFile->Get(histName.c_str())->Clone("recoHist"); assert(recoHist); recoHist->SetDirectory(0); recoFile->Close();
  
  assert(recoHist->GetBinWidth(1)==genHist->GetBinWidth(1));
  assert(recoHist->GetNbinsX()==genHist->GetNbinsX());
  int nbins = (fitMassHi-fitMassLo)/genHist->GetBinWidth(1);
  // gen,reco distributions in the time domain (mass)
  TH1D *hTimeGen  = new TH1D("hTimeGen" ,"hTimeGen" ,nbins,fitMassLo,fitMassHi); hTimeGen->Sumw2();
  TH1D *hTimeReco = new TH1D("hTimeReco","hTimeReco",nbins,fitMassLo,fitMassHi); hTimeReco->Sumw2();
  for(int i=1; i<=nbins; i++) {
    int j=genHist->FindBin(hTimeGen->GetBinCenter(i));
    hTimeGen->SetBinContent(i, genHist->GetBinContent(j)); hTimeGen->SetBinError(i, genHist->GetBinError(j));
    hTimeReco->SetBinContent(i, recoHist->GetBinContent(j)); hTimeReco->SetBinError(i, recoHist->GetBinError(j));
  }
  hTimeGen->Scale(1./hTimeGen->Integral()); hTimeReco->Scale(1./hTimeReco->Integral());
  
  // real and imaginary components
  vector<double> vFreqGenRe, vFreqGenIm, vFreqRecoRe, vFreqRecoIm;
  vFreqGenRe.reserve(nbins); vFreqGenIm.reserve(nbins);
  vFreqRecoRe.reserve(nbins); vFreqRecoIm.reserve(nbins);
  vector<TComplex> vFreqResZ; vFreqResZ.reserve(nbins);

  // gen,reco distributions in the frequency domain
  int r2rType=4;
  TH1 *hFreqGen=0, *hFreqReco=0; TVirtualFFT *fft;
  // r2r implementation
  //TVirtualFFT *fftFwdGen = TVirtualFFT::SineCosine(1,&nbins,&r2rType,"ex k"), // exhaustive, keep
  //            *fftFwdReco = TVirtualFFT::SineCosine(1,&nbins,&r2rType,"ex k"); 
  //for(int i=0; i<nbins; i++) {
  //  fftFwdGen->SetPoint(i, hTimeGen->GetBinContent(i+1));
  //  fftFwdReco->SetPoint(i, hTimeReco->GetBinContent(i+1));
  //}
  //fftFwdGen->Transform();
  //fftFwdReco->Transform();
  //hFreqGen =TH1::TransformHisto(fftFwdGen,hFreqGen, "MA"); hFreqGen->SetName("hFreqGen");
  //hFreqReco=TH1::TransformHisto(fftFwdReco,hFreqReco, "MA"); hFreqReco->SetName("hFreqReco");

  // r2c implementation

  hFreqGen=(TH1D*)hTimeGen->FFT(0, "Mag R2C EX"); hFreqGen->SetName("hFreqGen");
  fft=TVirtualFFT::GetCurrentTransform();
  for(int i=0; i<nbins; i++) fft->GetPointComplex(i, vFreqGenRe[i], vFreqGenIm[i]);
  hFreqReco=(TH1D*)hTimeReco->FFT(0, "Mag R2C EX"); hFreqReco->SetName("hFreqReco");
  fft=TVirtualFFT::GetCurrentTransform();
  for(int i=0; i<nbins; i++) fft->GetPointComplex(i, vFreqRecoRe[i], vFreqRecoIm[i]);
  TVirtualFFT *fftBack=TVirtualFFT::FFT(1,&nbins,"C2R M K");
  //TVirtualFFT *fftBack=TVirtualFFT::SineCosine(1,&nbins,&r2rType,"m k");
  
  TF1 *func=new TF1("func","[0]*exp(-0.5*((x-[1])/[2])**2)",0,nbins/2.);
  func->FixParameter(1,0.0);
  func->SetParameter(2,nbins/6.);
  hFreqReco->Fit(func,"MN0","",0,nbins/2.);
  double sigmaFreqReco=func->GetParameter(2);
  int fiveSigmaCutoff=(int)round(sigmaFreqReco*5.);
  delete func;

  for(int i=0; i<nbins; i++) {
    //fftBack->SetPoint(i, hFreqReco->GetBinContent(i+1)/hFreqGen->GetBinContent(i+1));
    // r2c implementation
    TComplex zReco(vFreqRecoRe[i],vFreqRecoIm[i]),
             zGen(vFreqGenRe[i], vFreqGenIm[i]);
    TComplex z(0,0);//TComplex z = zReco/zGen;
    if(i<fiveSigmaCutoff || i>=nbins-fiveSigmaCutoff) {
      z=zReco/zGen;
      fftBack->SetPointComplex(i, z); 
    }
    vFreqResZ[i]=z;
    printf("transfer z(%.2f,%.2f) = reco z(%.2f,%.2f) / gen z(%.2f,%.2f)\n", z.Re(), z.Im(), zReco.Re(), zReco.Im(), zGen.Re(), zGen.Im());
  } 
  fftBack->Transform();
  TH1 *resolutionFunction=0;
  resolutionFunction=TH1::TransformHisto(fftBack,resolutionFunction,"Re");
  //for(int i=1;i<=nbins; i++) resolutionFunction->SetBinContent(i,TMath::Abs(resolutionFunction->GetBinContent(i))); 
  
  // shift for plotting
  double xmin=-0.5/recoHist->GetBinWidth(1)*(1.-2./nbins);
  double xmax=0.5/recoHist->GetBinWidth(1);
  TH1D *hFreqGenScaled = new TH1D("hFreqGenScaled","hFreqGenScaled",nbins-1,xmin,xmax);
  TH1D *hFreqRecoScaled = new TH1D("hFreqRecoScaled","hFreqRecoScaled",nbins-1,xmin,xmax);
  TH1D *hFreqResScaled = new TH1D("hFreqResScaled","hFreqResScaled",nbins-1,xmin,xmax);
  TH1D *hResFuncShifted = new TH1D("hResFuncShifted","hResFuncShifted",nbins,-range/2.,range/2.);
  
  for(int i=2; i<=nbins; i++) { 
    int j;
    if(i<=nbins/2+1) j=i+nbins/2-2;
    else             j=i-nbins/2-1;
    hResFuncShifted->SetBinContent(j, resolutionFunction->GetBinContent(i));
    hFreqGenScaled->SetBinContent(j, hFreqGen->GetBinContent(i));
    hFreqRecoScaled->SetBinContent(j, hFreqReco->GetBinContent(i));
    hFreqResScaled->SetBinContent(j,vFreqResZ[i-1].Rho());
  }
  hFreqGenScaled->Scale(1./nbins);
  hFreqRecoScaled->Scale(1./nbins);
  hFreqResScaled->Scale(hFreqGenScaled->Integral()/hFreqResScaled->Integral());
  gStyle->SetOptStat(0); 
  TCanvas *cTimeDomain=new TCanvas("cTimeDomain","Time domain");
  hTimeGen->SetLineColor(kOrange-5); hTimeGen->SetFillColor(kOrange-5); hTimeGen->SetMarkerColor(kOrange-5);
  hTimeGen->SetFillStyle(3254);
  hTimeGen->SetLineWidth(2);
  hTimeGen->Draw();
  hTimeReco->SetLineColor(kViolet-1); hTimeReco->SetFillColor(kViolet-1); hTimeReco->SetMarkerColor(kViolet-1);
  hTimeReco->SetFillStyle(3254);
  hTimeReco->SetLineWidth(2);
  hTimeReco->Draw("SAME");
  TLegend *legendTime = new TLegend(.6,.6,.8,.8); 
  legendTime->AddEntry(hTimeGen , Form("Gen. mass (#mu=%.2f, #sigma=%.2f)",hTimeGen->GetMean(),hTimeGen->GetStdDev()), "lf");
  legendTime->AddEntry(hTimeReco, Form("Reco. mass (#mu=%.2f, #sigma=%.2f)",hTimeReco->GetMean(),hTimeReco->GetStdDev()), "lf");
  legendTime->SetFillColor(0); legendTime->Draw("SAME");
  TCanvas *cFreqDomain = new TCanvas("cFreqDomain","Frequency domain");
  hFreqGenScaled->SetLineColor(kOrange-5);
  hFreqGenScaled->SetLineWidth(2);
  hFreqGenScaled->Draw();
  hFreqRecoScaled->SetLineColor(kViolet-1);
  hFreqRecoScaled->SetLineWidth(2);
  hFreqRecoScaled->Draw("SAME");
  hFreqResScaled->SetLineColor(kBlue-6);
  hFreqResScaled->SetLineWidth(2);
  hFreqResScaled->Draw("SAME");
  hResFuncShifted->Scale(1./hResFuncShifted->Integral());
  TCanvas *cRes = new TCanvas("cRes","Resolution function");
  hResFuncShifted->Draw("SAME");

}

