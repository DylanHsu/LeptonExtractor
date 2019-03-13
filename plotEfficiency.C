#include "TROOT.h"
#include "TPaveText.h"
#include "TStyle.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TTree.h"
#include "TMath.h"
#include "TFile.h"
#include "TLine.h"
#include "TCanvas.h"
#include "TLegend.h"
#include <TRandom3.h>
#include <TEfficiency.h>

#include <vector>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
void plotEfficiency(std::string outputDir, std::string binFile) {
  system(Form("rm %s/*.junk 2>/dev/null",outputDir.c_str()));
  system(Form("for i in `ls %s/plots/fitres_pass* | sort -V`; do grep summary_ $i | tr '\n' ' ' | sed 's/summary_\\w*//g' | sed 's/[+-]//g' >> %s/passingEffs.junk; echo " " >> %s/passingEffs.junk; done",outputDir.c_str(),outputDir.c_str(),outputDir.c_str()));
  system(Form("for i in `ls %s/plots/fitres_fail* | sort -V`; do grep summary_ $i | tr '\n' ' ' | sed 's/summary_\\w*//g' | sed 's/[+-]//g' >> %s/failingEffs.junk; echo " " >> %s/failingEffs.junk; done",outputDir.c_str(),outputDir.c_str(),outputDir.c_str()));
  
  // parse bin file
  std::vector<double> fPtBinEdgesv, fEtaBinEdgesv, fPhiBinEdgesv, fNPVBinEdgesv, fJetsBinEdgesv, fMETBinEdgesv;
  // flags for |eta| and |phi| binning
  bool fDoAbsEta, fDoAbsPhi;
  // flags for binnings to compute efficiencies for
  bool fDoPt, fDoEta, fDoPhi, fDoEtaPt, fDoEtaPhi, fDoNPV, fDoJets, fDoMET;
  {
    std::ifstream ifs; ifs.open(binFile.c_str());  assert(ifs.is_open());
    std::string line;
    int state=0;
    int opts[8];
    while(getline(ifs,line)) {
      if(line[0]=='#') continue;
      if(line[0]=='%') { 
        state++; 
        continue; 
      }
      
      double edge;
      std::stringstream ss(line);
      if(state==0) {
        ss >> opts[0] >> opts[1] >> opts[2] >> opts[3] >> opts[4] >> opts[5] >> opts[6] >> opts[7];
        fDoPt     = (opts[0]==1);
        fDoEta    = (opts[1]==1);
        fDoPhi    = (opts[2]==1);
        fDoNPV    = (opts[3]==1);
        fDoJets   = (opts[4]==1);
        fDoMET    = (opts[5]==1);
        fDoEtaPt  = (opts[6]==1);
        fDoEtaPhi = (opts[7]==1);
      
      } else {
        ss >> edge;
        if     (state==1) { fPtBinEdgesv.push_back(edge);  }
        else if(state==2) { fDoAbsEta = (int(edge)==1); state++; }
        else if(state==3) { fEtaBinEdgesv.push_back(edge); }
        else if(state==4) { fDoAbsPhi = (int(edge)==1); state++; }
        else if(state==5) { fPhiBinEdgesv.push_back(edge); }
        else if(state==6) { fNPVBinEdgesv.push_back(edge); }
        else if(state==7) { fJetsBinEdgesv.push_back(edge); }
        else if(state==8) { fMETBinEdgesv.push_back(edge); }
      }
    }
    ifs.close();
  }
  //double ptBinEdges[fPtBinEdgesv.size()];   for(unsigned int i=0; i<fPtBinEdgesv.size();  i++) { ptBinEdges[i]  = fPtBinEdgesv[i];  }
  //double etaBinEdges[fEtaBinEdgesv.size()]; for(unsigned int i=0; i<fEtaBinEdgesv.size(); i++) { etaBinEdges[i] = fEtaBinEdgesv[i]; }
  
  TH2D *h_NsigPass = new TH2D("h_NsigPass","h_NsigPass",fEtaBinEdgesv.size()-1,fEtaBinEdgesv.data(),fPtBinEdgesv.size()-1,fPtBinEdgesv.data());
  TH2D *h_NsigFail = new TH2D("h_NsigFail","h_NsigFail",fEtaBinEdgesv.size()-1,fEtaBinEdgesv.data(),fPtBinEdgesv.size()-1,fPtBinEdgesv.data());
  TH2D *hEffEtaPt = new TH2D("hEffEtaPt","hEffEtaPt",fEtaBinEdgesv.size()-1,fEtaBinEdgesv.data(),fPtBinEdgesv.size()-1,fPtBinEdgesv.data());
  TH2D *hErrhEtaPt = new TH2D("hErrhEtaPt","hErrhEtaPt",fEtaBinEdgesv.size()-1,fEtaBinEdgesv.data(),fPtBinEdgesv.size()-1,fPtBinEdgesv.data());
  TH2D *hErrlEtaPt = new TH2D("hErrlEtaPt","hErrlEtaPt",fEtaBinEdgesv.size()-1,fEtaBinEdgesv.data(),fPtBinEdgesv.size()-1,fPtBinEdgesv.data());
  TH2D *hEdmPass = new TH2D("hEdmPass","hEdmPass",fEtaBinEdgesv.size()-1,fEtaBinEdgesv.data(),fPtBinEdgesv.size()-1,fPtBinEdgesv.data());
  TH2D *hEdmFail = new TH2D("hEdmFail","hEdmFail",fEtaBinEdgesv.size()-1,fEtaBinEdgesv.data(),fPtBinEdgesv.size()-1,fPtBinEdgesv.data());
  TH2D *hFitProbPass = new TH2D("hFitProbPass","hFitProbPass",fEtaBinEdgesv.size()-1,fEtaBinEdgesv.data(),fPtBinEdgesv.size()-1,fPtBinEdgesv.data());
  TH2D *hFitProbFail = new TH2D("hFitProbFail","hFitProbFail",fEtaBinEdgesv.size()-1,fEtaBinEdgesv.data(),fPtBinEdgesv.size()-1,fPtBinEdgesv.data());
  std::ifstream pass_ifs; pass_ifs.open(Form("%s/passingEffs.junk",outputDir.c_str())); assert(pass_ifs.is_open());
  std::ifstream fail_ifs; fail_ifs.open(Form("%s/failingEffs.junk",outputDir.c_str())); assert(fail_ifs.is_open());
  std::string inputLine;
  
  unsigned long int time_now = static_cast<unsigned long int>(time(NULL));
  unsigned int randomToySeed=(time_now-731178000); // random seed based on Dylan's age in seconds
  TRandom3 toymaker(randomToySeed);
  
  {// fill the histogram from the text file
   double NsigPass, NsigErrhPass, NsigErrlPass, edmPass, fitProbPass;
   double NsigFail, NsigErrhFail, NsigErrlFail, edmFail, fitProbFail;
   double NsigTotal;
   bool gotLinePass, gotLineFail;
   for(unsigned int iPt=1;iPt<fPtBinEdgesv.size(); iPt++) { for(unsigned int iEta=1; iEta<fEtaBinEdgesv.size(); iEta++) { 
    if(getline(pass_ifs,inputLine)) gotLinePass=true;
    std::stringstream ss_pass(inputLine);
    if(getline(fail_ifs,inputLine)) gotLineFail=true;
    std::stringstream ss_fail(inputLine);
    ss_pass >> NsigPass >> NsigErrhPass >> NsigErrlPass >> edmPass >> fitProbPass;
    ss_fail >> NsigFail >> NsigErrhFail >> NsigErrlFail >> edmFail >> fitProbFail;
    int mcBin=h_NsigPass->GetBin(iEta,iPt);
    double NsigTotal = NsigPass+NsigFail;
    double eff = NsigPass/NsigTotal;
    double effLo = (NsigPass - NsigErrlPass) / (NsigPass + NsigFail - NsigErrlPass + NsigErrhFail);
    double effHi = (NsigPass + NsigErrhPass) / (NsigPass + NsigFail + NsigErrlPass - NsigErrlFail);
    // conservatively estimate the variance of efficiency
    double effEstError = TMath::Max( fabs(effLo-eff), fabs(effHi-eff));
    TH1F *toyEffs = new TH1F("toyEffs","toyEffs", 100, TMath::Max(0.,eff-5*effEstError), TMath::Min(1.,eff+5*effEstError));
    //printf("toyEffs %.3f to %.3f\n", toyEffs->GetBinLowEdge(1), toyEffs->GetBinLowEdge(101));
    for(int iToy=0; iToy<1000; iToy++) {
      // Toys are Gaussian about 0 with sigma 1
      double toyPass=toymaker.Gaus(0,1);
      double toyFail=toymaker.Gaus(0,1);
      // Use the upper or lower error bar on the signal size based on the sign of the toys
      double toyNPass = NsigPass + toyPass * (toyPass>0? NsigErrhPass : NsigErrlPass);
      double toyNFail = NsigFail + toyFail * (toyFail>0? NsigErrhFail : NsigErrlFail);
      double toyEff = toyNPass/(toyNPass+toyNFail);
      //printf("toyEff #%d: %.3f\n", iToy, toyEff);
      toyEffs->Fill(toyEff);
    }
    double quantileProbs[3]={0.159,0.5,0.841};
    double theQuantiles[3];
    toyEffs->GetQuantiles(3, theQuantiles, quantileProbs);
    delete toyEffs;
    double errLo = TMath::Max(0.,eff - TMath::Max(0.,theQuantiles[0]));
    double errHi = TMath::Max(0.,TMath::Min(1.,theQuantiles[2]) - eff);
    h_NsigPass->SetBinContent(mcBin, NsigPass); h_NsigPass->SetBinError(mcBin, TMath::Max(NsigErrlPass,NsigErrhPass));
    h_NsigFail->SetBinContent(mcBin, NsigFail); h_NsigFail->SetBinError(mcBin, TMath::Max(NsigErrhFail,NsigErrlFail));
    hEffEtaPt->SetBinContent(mcBin, eff); hEffEtaPt->SetBinError(mcBin, TMath::Max(errLo,errHi));
    hErrhEtaPt->SetBinContent(mcBin, errHi); hErrlEtaPt->SetBinContent(mcBin, errLo);
    hEdmPass->SetBinContent(mcBin, TMath::Min(1.,TMath::Max(1e-6, edmPass))); hEdmFail->SetBinContent(mcBin, TMath::Min(1.,TMath::Max(1e-6, edmFail)));
    hFitProbPass->SetBinContent(mcBin, TMath::Max(1e-4,fitProbPass));    hFitProbFail->SetBinContent(mcBin, TMath::Max(1e-4,fitProbFail));
   }}
  }
  pass_ifs.close();
  fail_ifs.close();
  gStyle->SetPaintTextFormat("4.2f"); gStyle->SetOptStat(0); gStyle->SetPalette(kBlackBody);
  TH2D *theHist; TCanvas *theCanvas;
  theCanvas=new TCanvas("cEffEtaPt","cEffEtaPt"); theHist=hEffEtaPt;
  theCanvas->SetLogy();
  theHist->SetTitle("");
  theHist->SetMinimum(0);
  theHist->SetMaximum(1.);
  theHist->GetXaxis()->SetTitle("probe #eta");
  theHist->GetYaxis()->SetTitle("probe p_{T}");
  theHist->GetYaxis()->SetMoreLogLabels();
  theHist->GetYaxis()->SetRangeUser(10,100);
  //theHist->Draw("COLZ TEXT45");
  theHist->Draw("COLZ");
  theCanvas->Print(Form("%s/cEffEtaPt.png",outputDir.c_str()));
  theCanvas=new TCanvas("cFitProbPass","cFitProbPass"); theHist=hFitProbPass;
  theCanvas->SetLogy(); theCanvas->SetLogz();
  theHist->SetTitle("");
  theHist->GetXaxis()->SetTitle("probe #eta");
  theHist->GetYaxis()->SetTitle("probe p_{T}");
  theHist->GetYaxis()->SetMoreLogLabels();
  theHist->GetYaxis()->SetRangeUser(10,100);
  theHist->SetMinimum(1e-4);
  theHist->SetMaximum(1.);
  theHist->Draw("COLZ");
  theCanvas->Print(Form("%s/cFitProbPass.png",outputDir.c_str()));
  theCanvas=new TCanvas("cFitProbFail","cFitProbFail"); theHist=hFitProbFail;
  theCanvas->SetLogy(); theCanvas->SetLogz();
  theHist->SetTitle("");
  theHist->GetXaxis()->SetTitle("probe #eta");
  theHist->GetYaxis()->SetTitle("probe p_{T}");
  theHist->GetYaxis()->SetMoreLogLabels();
  theHist->GetYaxis()->SetRangeUser(10,100);
  theHist->SetMinimum(1e-4);
  theHist->SetMaximum(1.);
  theHist->Draw("COLZ");
  theCanvas->Print(Form("%s/cFitProbFail.png",outputDir.c_str()));
  theCanvas=new TCanvas("cEdmPass","cEdmPass"); theHist=hEdmPass;
  theCanvas->SetLogy(); theCanvas->SetLogz();
  theHist->SetTitle("");
  theHist->GetXaxis()->SetTitle("probe #eta");
  theHist->GetYaxis()->SetTitle("probe p_{T}");
  theHist->GetYaxis()->SetMoreLogLabels();
  theHist->GetYaxis()->SetRangeUser(10,100);
  theHist->SetMinimum(1e-6);
  theHist->SetMaximum(1.);
  theHist->Draw("COLZ");
  theCanvas->Print(Form("%s/cEdmPass.png",outputDir.c_str()));
  theCanvas=new TCanvas("cEdmFail","cEdmFail"); theHist=hEdmFail;
  theCanvas->SetLogy(); theCanvas->SetLogz();
  theHist->SetTitle("");
  theHist->GetXaxis()->SetTitle("probe #eta");
  theHist->GetYaxis()->SetTitle("probe p_{T}");
  theHist->GetYaxis()->SetMoreLogLabels();
  theHist->GetYaxis()->SetRangeUser(10,100);
  theHist->SetMinimum(1e-6);
  theHist->SetMaximum(1.);
  theHist->Draw("COLZ");
  theCanvas->Print(Form("%s/cEdmFail.png",outputDir.c_str()));

  TFile *rootfile=TFile::Open(Form("%s/eff.root",outputDir.c_str()),"recreate"); assert(rootfile);
  hEffEtaPt    ->Write();
  hErrhEtaPt   ->Write();
  hErrlEtaPt   ->Write();
  hEdmPass     ->Write();
  hEdmFail     ->Write();
  hFitProbPass ->Write();
  hFitProbFail ->Write();
  rootfile->Close();
  //system(Form("rm %s/*.junk 2>/dev/null",outputDir.c_str()));
}
