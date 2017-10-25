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
#include <TEfficiency.h>

#include <vector>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
const int fitMassLo=60;
const int fitMassHi=120;
void cutAndCount(std::string outputDir, std::string binFile, std::string templateFileName) {
  system(Form("mkdir -p %s",outputDir.c_str()));
  TFile *templateFile=TFile::Open(templateFileName.c_str(),"READ"); assert(templateFile);
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
  
  TH2D *h_NsigPass = new TH2D("h_NsigPass","h_NsigPass",fEtaBinEdgesv.size()-1,fEtaBinEdgesv.data(),fPtBinEdgesv.size()-1,fPtBinEdgesv.data());
  TH2D *h_NsigFail = new TH2D("h_NsigFail","h_NsigFail",fEtaBinEdgesv.size()-1,fEtaBinEdgesv.data(),fPtBinEdgesv.size()-1,fPtBinEdgesv.data());
  TH2D *hEffEtaPt = new TH2D("hEffEtaPt","hEffEtaPt",fEtaBinEdgesv.size()-1,fEtaBinEdgesv.data(),fPtBinEdgesv.size()-1,fPtBinEdgesv.data());
  TH2D *hErrhEtaPt = new TH2D("hErrhEtaPt","hErrhEtaPt",fEtaBinEdgesv.size()-1,fEtaBinEdgesv.data(),fPtBinEdgesv.size()-1,fPtBinEdgesv.data());
  TH2D *hErrlEtaPt = new TH2D("hErrlEtaPt","hErrlEtaPt",fEtaBinEdgesv.size()-1,fEtaBinEdgesv.data(),fPtBinEdgesv.size()-1,fPtBinEdgesv.data());
  
  {// get the histogram from the root file
   double NsigPass, NsigErrPass;
   double NsigFail, NsigErrFail;
   double NsigTotal, eff, effLo, effHi, errLo, errHi;
   bool gotLinePass, gotLineFail;
   for(unsigned int iPt=1;iPt<fPtBinEdgesv.size(); iPt++) { for(unsigned int iEta=1; iEta<fEtaBinEdgesv.size(); iEta++) { 
    NsigPass=0;
    NsigFail=0;
    NsigErrPass=0;
    NsigErrFail=0;
    TH1D *histPass=(TH1D*)templateFile->Get(Form("pass_ptBin%d_etaBin%d",iPt-1,iEta-1));
    TH1D *histFail=(TH1D*)templateFile->Get(Form("fail_ptBin%d_etaBin%d",iPt-1,iEta-1));
    assert(histPass); assert(histFail);
    unsigned massBin1=histPass->FindBin(fitMassLo+0.01);
    unsigned massBin2=histPass->FindBin(fitMassHi-0.01);
    for(unsigned iM=massBin1; iM<=massBin2; iM++) {
     NsigPass+=histPass->GetBinContent(iM); NsigFail+=histFail->GetBinContent(iM);
     NsigErrPass+=pow(histPass->GetBinError(iM),2); NsigErrFail+=pow(histFail->GetBinError(iM),2);
    }
    NsigErrPass=sqrt(NsigErrPass); NsigErrFail=sqrt(NsigErrFail);
    int nBin=h_NsigPass->GetBin(iEta,iPt);
    NsigTotal = NsigPass+NsigFail; eff = NsigPass/NsigTotal;
    effLo = (NsigPass - NsigErrPass) / (NsigPass + NsigFail - NsigErrPass + NsigErrFail);
    effHi = (NsigPass + NsigErrPass) / (NsigPass + NsigFail + NsigErrPass - NsigErrFail);
    errLo = eff - TMath::Max(0.,effLo);
    errHi = TMath::Min(1.,effHi) - eff;
    h_NsigPass->SetBinContent(nBin, NsigPass); h_NsigPass->SetBinError(nBin, NsigErrPass);
    h_NsigFail->SetBinContent(nBin, NsigFail); h_NsigFail->SetBinError(nBin, NsigErrFail);
    hEffEtaPt->SetBinContent(nBin, eff); hEffEtaPt->SetBinError(nBin, TMath::Max(errLo,errHi));
    hErrhEtaPt->SetBinContent(nBin, errHi); hErrlEtaPt->SetBinContent(nBin, errLo);
   }}
  }
  gStyle->SetPaintTextFormat("4.2f"); gStyle->SetOptStat(0);
  gStyle->SetPalette(kBlackBody);
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

  TFile *rootfile=TFile::Open(Form("%s/eff.root",outputDir.c_str()),"recreate"); assert(rootfile);
  hEffEtaPt    ->Write();
  hErrhEtaPt   ->Write();
  hErrlEtaPt   ->Write();
  rootfile->Close();
}
