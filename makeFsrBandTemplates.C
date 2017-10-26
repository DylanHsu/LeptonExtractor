#include "TROOT.h"
#include "TPaveText.h"
#include "TStyle.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TTree.h"
#include "TMath.h"
#include "TFile.h"
#include "TCanvas.h"
#include <TEfficiency.h>

#include <vector>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <boost/algorithm/string.hpp> // include Boost, a C++ library

const int fitMassLo=60;
const int fitMassHi=120;
void makeFsrBandTemplates(
  std::string templateFileName, 
  std::string genZmmFileName, 
  std::string genZeeFileName, 
  std::string binFile
) {
  TH1::AddDirectory(kFALSE);
  TFile *templateFile=TFile::Open(templateFileName.c_str(),"READ"); assert(templateFile);
  TFile *genZmmFile=TFile::Open(genZmmFileName.c_str(),"READ"); assert(genZmmFile);
  TFile *genZeeFile=TFile::Open(genZeeFileName.c_str(),"READ"); assert(genZeeFile);
  std::string fsrUpFileName=templateFileName; boost::replace_all(fsrUpFileName,".root","_fsrUp.root");
  std::string fsrDownFileName=templateFileName; boost::replace_all(fsrDownFileName,".root","_fsrDown.root");
  if(templateFileName==fsrUpFileName || templateFileName==fsrDownFileName) { assert(0); return;}
  TFile *fsrUpFile=TFile::Open(fsrUpFileName.c_str(),"recreate"); assert(fsrUpFile);
  TFile *fsrDownFile=TFile::Open(fsrDownFileName.c_str(),"recreate"); assert(fsrDownFile);
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
  TH2D *h_fsrChangePass = new TH2D("h_fsrChangePass","h_fsrChangePass",fEtaBinEdgesv.size()-1,fEtaBinEdgesv.data(),fPtBinEdgesv.size()-1,fPtBinEdgesv.data());
  TH2D *h_fsrChangeFail = new TH2D("h_fsrChangeFail","h_fsrChangeFail",fEtaBinEdgesv.size()-1,fEtaBinEdgesv.data(),fPtBinEdgesv.size()-1,fPtBinEdgesv.data());
  for(unsigned int iPt=1;iPt<fPtBinEdgesv.size(); iPt++) { 
    TH1D *genZmm=0; TH1D *genZee=0;
    //unsigned m1=10,m2=40;
    /*for(unsigned int iEta=1; iEta<fEtaBinEdgesv.size(); iEta++) { 
      if(!genZmm) { genZmm=(TH1D*)genZmmFile->Get(Form("pass_ptBin%d_etaBin%d",iPt-1,iEta-1)); assert(genZmm); }
      else genZmm->Add((TH1D*)genZmmFile->Get(Form("pass_ptBin%d_etaBin%d",iPt-1,iEta-1)));
      genZmm->Add((TH1D*)genZmmFile->Get(Form("fail_ptBin%d_etaBin%d",iPt-1,iEta-1)));
      if(!genZee) { genZee=(TH1D*)genZeeFile->Get(Form("pass_ptBin%d_etaBin%d",iPt-1,iEta-1)); assert(genZee); }
      else genZee->Add((TH1D*)genZeeFile->Get(Form("pass_ptBin%d_etaBin%d",iPt-1,iEta-1)));
      genZee->Add((TH1D*)genZeeFile->Get(Form("fail_ptBin%d_etaBin%d",iPt-1,iEta-1)));
    }
    //genZmm->Scale(1./genZmm->Integral(m1,m2));
    //genZee->Scale(1./genZee->Integral(m1,m2));
    */
    for(unsigned int iEta=1; iEta<fEtaBinEdgesv.size(); iEta++) { 
      // get passing and failing reconstructed templates
      TH1D *histPass=(TH1D*)templateFile->Get(Form("pass_ptBin%d_etaBin%d",iPt-1,iEta-1)); assert(histPass);
      TH1D *histFail=(TH1D*)templateFile->Get(Form("fail_ptBin%d_etaBin%d",iPt-1,iEta-1)); assert(histFail);
      unsigned m1=histPass->FindBin(fitMassLo+0.01);
      unsigned m2=histPass->FindBin(fitMassHi-0.01);
      unsigned mLowMassCutoff=histPass->FindBin(79.999);
      genZmm=(TH1D*)genZmmFile->Get(Form("pass_ptBin%d_etaBin%d",iPt-1,iEta-1)); assert(genZmm);
      genZmm->Add((TH1D*)genZmmFile->Get(Form("fail_ptBin%d_etaBin%d",iPt-1,iEta-1)));
      genZee=(TH1D*)genZeeFile->Get(Form("pass_ptBin%d_etaBin%d",iPt-1,iEta-1)); assert(genZee);
      genZee->Add((TH1D*)genZeeFile->Get(Form("fail_ptBin%d_etaBin%d",iPt-1,iEta-1)));
      genZmm->Scale(1./genZmm->Integral(m1,m2));
      genZee->Scale(1./genZee->Integral(m1,m2));
      // get the generator shapes. add passing and failing together
      //TH1D *genZmmPass=(TH1D*)genZmmFile->Get(Form("pass_ptBin%d_etaBin%d",iPt-1,iEta-1)); assert(genZmmPass);
      //TH1D *genZmm    =(TH1D*)genZmmFile->Get(Form("fail_ptBin%d_etaBin%d",iPt-1,iEta-1)); assert(genZmm);
      //genZmm->Add(genZmmPass); genZmm->Scale(1./genZmm->Integral(m1,m2));
      //TH1D *genZeePass=(TH1D*)genZeeFile->Get(Form("pass_ptBin%d_etaBin%d",iPt-1,iEta-1)); assert(genZeePass);
      //TH1D *genZee    =(TH1D*)genZeeFile->Get(Form("fail_ptBin%d_etaBin%d",iPt-1,iEta-1)); assert(genZee);
      //genZee->Add(genZeePass); genZee->Scale(1./genZee->Integral(m1,m2));
      // Take the uncertainty band as (0.9 genZmm + 0.1 genZee) / (1.0 genZmm);
      TH1D *fsrUncBand = (TH1D*)genZmm->Clone("fsrUncBand");
      //fsrUncBand->Scale(0.9); fsrUncBand->Add(genZee,.1); fsrUncBand->Divide(genZmm);
      fsrUncBand->Scale(0.0); fsrUncBand->Add(genZee,1.0); fsrUncBand->Divide(genZmm);
      // Do a cut off so there are no crazy values
      for(unsigned m=m1;m<=m2;m++) {
        if     (fsrUncBand->GetBinContent(m)>5  ) fsrUncBand->SetBinContent(m,5);
        else if(fsrUncBand->GetBinContent(m)<0.2) fsrUncBand->SetBinContent(m,0.2);
      }
      // Wag the tail
      TH1D *histPassFsrUp=(TH1D*)histPass->Clone("histPassFsrUp");
      TH1D *histPassFsrDown=(TH1D*)histPass->Clone("histPassFsrDown");
      TH1D *histFailFsrUp=(TH1D*)histFail->Clone("histFailFsrUp");
      TH1D *histFailFsrDown=(TH1D*)histFail->Clone("histFailFsrDown");
      histPassFsrUp->Multiply(fsrUncBand);
      histPassFsrDown->Divide(fsrUncBand);
      histFailFsrUp->Multiply(fsrUncBand);
      histFailFsrDown->Divide(fsrUncBand);
      //histPassFsrUp->Scale(histPass->Integral(m1,m2) / histPassFsrUp->Integral(m1,m2));
      //histPassFsrDown->Scale(histPass->Integral(m1,m2) / histPassFsrDown->Integral(m1,m2));
      //histFailFsrUp->Scale(histFail->Integral(m1,m2) / histFailFsrUp->Integral(m1,m2));
      //histFailFsrDown->Scale(histFail->Integral(m1,m2) / histFailFsrDown->Integral(m1,m2));
      int nb2d=h_fsrChangePass->GetBin(iEta,iPt);
      double fsrChangePass=(histPassFsrUp->Integral(m1,mLowMassCutoff) - histPass->Integral(m1,mLowMassCutoff))/histPass->Integral(m1,mLowMassCutoff);
      double fsrChangeFail=(histFailFsrUp->Integral(m1,mLowMassCutoff) - histFail->Integral(m1,mLowMassCutoff))/histFail->Integral(m1,mLowMassCutoff);
      //double fsrChangePass=fsrUncBand->Integral(m1,mLowMassCutoff) / (fsrUncBand->GetBinLowEdge(mLowMassCutoff+1)-fsrUncBand->GetBinLowEdge(m1));
      //double fsrChangeFail=fsrUncBand->Integral(m1,mLowMassCutoff) / (fsrUncBand->GetBinLowEdge(mLowMassCutoff+1)-fsrUncBand->GetBinLowEdge(m1));
      h_fsrChangePass->SetBinContent(nb2d, fsrChangePass);
      h_fsrChangeFail->SetBinContent(nb2d, fsrChangeFail);
      printf("iEta %d, iPt %d, fsrChangePass %f, fsrChangeFail %f\n",iEta-1,iPt-1,fsrChangePass,fsrChangeFail);
      fsrUpFile->cd();
      histPassFsrUp->Write(histPass->GetName()); histFailFsrUp->Write(histFail->GetName());
      fsrDownFile->cd();
      histPassFsrDown->Write(histPass->GetName()); histFailFsrDown->Write(histFail->GetName());
    }
  }
  gStyle->SetPalette(kCool);gStyle->SetOptStat(0);
  TCanvas *c_fsrChangePass=new TCanvas("c_fsrChangePass","c_fsrChangePass");
  h_fsrChangePass->SetMinimum(0);
  h_fsrChangePass->SetMaximum(.2);
  h_fsrChangePass->GetYaxis()->SetRangeUser(10,45);
  h_fsrChangePass->Draw("colz");
  c_fsrChangePass->Print("fsrChangePass.pdf");
  TCanvas *c_fsrChangeFail=new TCanvas("c_fsrChangeFail","c_fsrChangeFail");
  h_fsrChangeFail->SetMinimum(0);
  h_fsrChangeFail->SetMaximum(.2);
  h_fsrChangeFail->GetYaxis()->SetRangeUser(10,45);
  h_fsrChangeFail->Draw("colz");
  c_fsrChangeFail->Print("fsrChangeFail.pdf");
  fsrUpFile->cd();
  h_fsrChangePass->Write(); h_fsrChangeFail->Write();
  fsrDownFile->cd();
  h_fsrChangePass->Write(); h_fsrChangeFail->Write();
  templateFile->Close();
  genZmmFile->Close();
  genZeeFile->Close();
  fsrUpFile->Close();
  fsrDownFile->Close();
}
