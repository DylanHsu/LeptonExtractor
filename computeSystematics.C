#include <TROOT.h>
#include <TMath.h>
#include <TChain.h>
#include <TFile.h>
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include <vector>
#include <sstream>
#include <cassert>
#include <iomanip>
#include <fstream>
#include <TH1.h>
#include <TH1F.h>
#include <TH2.h>
#include <THStack.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TGraphAsymmErrors.h>
#include <TLegend.h>
#include <TFile.h>
#include <TCut.h>
#include <TColor.h>
#include <TPaletteAxis.h>
#include <iostream>

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void computeSystematics(
  TString nominalFileName,     // root file with nominal scale factors
  TString altFileName,         // root file with the alternative scale factors
  TString sfName,              // name of the scale factors histo in the root file
  TString systName,            // name of the systematic source
  bool replaceExisting=true,   // whether to replace the existing calculation (if not, take the max)
  TString plotTitle="",
  int width=2400,
  int height=1200,
  double xmin=0,
  double xmax=0,
  double ymin=0,
  double ymax=0,
  bool plotOnly=true
) {
  system("mkdir -p plots");
  TFile *altFile     = TFile::Open(altFileName, "read"); assert(altFile);
  TFile *nominalFile = TFile::Open(nominalFileName, "update"); assert(nominalFile);
  TH2D *nominalSF = (TH2D*)nominalFile->Get(sfName)->Clone("nominal sf"); assert(nominalSF); nominalSF->SetDirectory(0);
  TH2D *altSF     = (TH2D*)altFile    ->Get(sfName)->Clone("alt sf"    ); assert(altSF    ); altSF    ->SetDirectory(0);
  TString systHistoName = Form("%s_%s", sfName.Data(), systName.Data());
  TH2D *h_syst=0; TH1D *h_pulls=0;
  if(!replaceExisting) {
    h_pulls=(TH1D*)nominalFile->Get(Form("pulls_%s",systName.Data())); 
    h_syst=(TH2D*)nominalFile->Get(systHistoName);
    if(!h_syst || !h_pulls) { printf("Warning: existing systematics was not in file to take the maximum with, we have to make it from scratch\n"); replaceExisting=true; }
    else { h_pulls->SetDirectory(0); h_syst->SetDirectory(0); }
  }
  if(replaceExisting) { 
    h_syst=(TH2D*)nominalSF->Clone(systHistoName); h_syst->SetDirectory(0); h_syst->Reset(); h_syst->Clear(); h_syst->SetTitle(plotTitle); h_syst->SetName(systHistoName);
    h_pulls=new TH1D(Form("pulls_%s",systName.Data()),Form("Pulls"),20,-5,5);
  }
  assert(h_syst); assert(h_pulls);
  nominalSF->GetXaxis()->UnZoom(); nominalSF->GetYaxis()->UnZoom();
  if(!plotOnly || replaceExisting) { for(unsigned iEta=1; iEta<=(unsigned)nominalSF->GetNbinsX(); iEta++) { for(unsigned iPt=1; iPt<=(unsigned)nominalSF->GetNbinsY(); iPt++) {
    unsigned nb2d = h_syst->GetBin(iEta,iPt);
    double syst = TMath::Abs(nominalSF->GetBinContent(nb2d) - altSF->GetBinContent(nb2d));
    double pull = (altSF->GetBinContent(nb2d) - nominalSF->GetBinContent(nb2d)) / nominalSF->GetBinError(nb2d);
    //printf("iEta %d iPt %d nominal %f alt %f syst %f\n",iEta,iPt,nominalSF->GetBinContent(nb2d),altSF->GetBinContent(nb2d),syst);
    if(!replaceExisting) syst=TMath::Max( h_syst->GetBinContent(nb2d), syst );
    h_syst->SetBinContent(nb2d, syst);
    h_pulls->Fill(pull);
  }}}
  h_syst->GetXaxis()->SetTitle("#eta");
  h_syst->GetXaxis()->SetTitleOffset(1.1);
  h_syst->GetXaxis()->SetTitleSize(0.04);
  h_syst->GetXaxis()->SetLabelSize(0.04);
  h_syst->GetYaxis()->SetTitle("p_{T} [GeV]");
  h_syst->GetYaxis()->SetTitleOffset(1.1);
  h_syst->GetYaxis()->SetTitleSize(0.04);
  h_syst->GetYaxis()->SetLabelSize(0.04);
  h_syst->SetMinimum(-0.0001); h_syst->SetMaximum(0.05);
  gStyle->SetOptStat(0); gStyle->SetPalette(kCool); gStyle->SetPaintTextFormat("4.3f");
  TCanvas *c_syst=new TCanvas("c_syst","c_syst",width,height); h_syst->Draw("colz text45"); 
  if(!plotOnly) { h_syst->Write(systHistoName, TObject::kOverwrite); h_pulls->Write(h_pulls->GetName(), TObject::kOverwrite); }

  TLine *lines[1000]; int nLines=0;
  bool xRange=(xmin!=xmax); bool yRange=(ymin!=ymax);
  if(xRange) h_syst->GetXaxis()->SetRangeUser(xmin,xmax); else h_syst->GetXaxis()->SetRangeUser(h_syst->GetXaxis()->GetBinLowEdge(1), h_syst->GetXaxis()->GetBinLowEdge(h_syst->GetXaxis()->GetNbins()+1));
  if(yRange) h_syst->GetYaxis()->SetRangeUser(ymin,ymax); else h_syst->GetYaxis()->SetRangeUser(h_syst->GetYaxis()->GetBinLowEdge(1), h_syst->GetYaxis()->GetBinLowEdge(h_syst->GetYaxis()->GetNbins()+1));
  for(unsigned i=(xRange ? h_syst->GetXaxis()->FindBin(xmin+.001):1); i<(unsigned)(xRange ? h_syst->GetXaxis()->FindBin(xmax-.001) : h_syst->GetXaxis()->GetNbins()); i++) {
    lines[nLines]=new TLine(
      h_syst->GetXaxis()->GetBinUpEdge(i),   
      yRange? h_syst->GetYaxis()->GetBinUpEdge(h_syst->GetYaxis()->FindBin(ymin+.001)-1) : h_syst->GetYaxis()->GetBinUpEdge(0),
      h_syst->GetXaxis()->GetBinUpEdge(i),
      yRange? h_syst->GetYaxis()->GetBinUpEdge(h_syst->GetYaxis()->FindBin(ymax-.001)) : h_syst->GetYaxis()->GetBinUpEdge(h_syst->GetYaxis()->GetNbins())
    ); lines[nLines]->SetLineStyle(kDashed); nLines++;
  }
  for(unsigned j=(yRange ? h_syst->GetYaxis()->FindBin(ymin+.001):1); j<(unsigned)(yRange ? h_syst->GetYaxis()->FindBin(ymax-.001) : h_syst->GetYaxis()->GetNbins()); j++) {
    lines[nLines]=new TLine(
      xRange? h_syst->GetXaxis()->GetBinUpEdge(h_syst->GetXaxis()->FindBin(xmin+.001)-1) : h_syst->GetXaxis()->GetBinUpEdge(0),
      h_syst->GetYaxis()->GetBinUpEdge(j),   
      xRange? h_syst->GetXaxis()->GetBinUpEdge(h_syst->GetXaxis()->FindBin(xmax-.001)) : h_syst->GetXaxis()->GetBinUpEdge(h_syst->GetXaxis()->GetNbins()),
      h_syst->GetYaxis()->GetBinUpEdge(j)
    ); lines[nLines]->SetLineStyle(kDashed); nLines++;
  }
   
  if(xRange) h_syst->GetXaxis()->SetRangeUser(xmin,xmax); else h_syst->GetXaxis()->SetRangeUser(h_syst->GetXaxis()->GetBinLowEdge(1), h_syst->GetXaxis()->GetBinLowEdge(h_syst->GetXaxis()->GetNbins()+1));
  if(yRange) h_syst->GetYaxis()->SetRangeUser(ymin,ymax); else h_syst->GetYaxis()->SetRangeUser(h_syst->GetYaxis()->GetBinLowEdge(1), h_syst->GetYaxis()->GetBinLowEdge(h_syst->GetYaxis()->GetNbins()+1));
  h_syst->SetMarkerSize(1.1);
  h_syst->Draw("TEXT45 colz");
  c_syst->Update();
  TPaletteAxis *palette_axis = (TPaletteAxis*) h_syst->GetListOfFunctions()->FindObject("palette"); 
  palette_axis->SetLabelSize(0.02);
  c_syst->Update();
  for(int nl=0;nl<nLines;nl++) lines[nl]->Draw("SAME");
  c_syst->Print(Form("plots/%s_%.1f-%.1f_%.1f-%.1f.png",h_syst->GetName(),xmin,xmax,ymin,ymax));
  c_syst->Print(Form("plots/%s_%.1f-%.1f_%.1f-%.1f.pdf",h_syst->GetName(),xmin,xmax,ymin,ymax));
  gStyle->SetOptStat(110001111);
  TCanvas *c_pulls=new TCanvas("c_pulls","c_pulls",height,height);
  h_pulls->GetXaxis()->SetTitle("Difference / Stat. Unc.");
  h_pulls->SetLineColor(kViolet-4); h_pulls->SetFillColor(kViolet-4); h_pulls->SetFillStyle(1001); h_pulls->Draw("HIST");
  c_pulls->Print(Form("plots/pulls_%s.png",h_syst->GetName()));
  c_pulls->Print(Form("plots/pulls_%s.pdf",h_syst->GetName()));
  nominalFile->Close(); altFile->Close();
}
void computeAll() {
  computeSystematics("2018-02-09/trackToMediumMuonTightIso/sfMediumMuonBtoFNominal.root", "2018-02-09/trackToMediumMuonTightIso/sfMediumMuonBtoFAltBkg.root" , "scalefactors_Medium2016MuonBCDEF", "AltBkg"     , true , "Systematic from background modeling", 1200,500, 0,0, 10,120, false);  
  computeSystematics("2018-02-09/trackToMediumMuonTightIso/sfMediumMuonBtoFNominal.root", "2018-02-09/trackToMediumMuonTightIso/sfMediumMuonBtoFAltGen.root" , "scalefactors_Medium2016MuonBCDEF", "AltGen"     , true , "Systematic from generator choice"   , 1200,500, 0,0, 10,120, false);  
  computeSystematics("2018-02-09/trackToMediumMuonTightIso/sfMediumMuonBtoFNominal.root", "2018-02-09/trackToMediumMuonTightIso/sfMediumMuonBtoFAltTag.root" , "scalefactors_Medium2016MuonBCDEF", "AltTag"     , true , "Systematic from tag selection bias" , 1200,500, 0,0, 10,120, false);  
  computeSystematics("2018-02-09/trackToMediumMuonTightIso/sfMediumMuonBtoFNominal.root", "2018-02-09/trackToMediumMuonTightIso/sfMediumMuonBtoFfsrDown.root", "scalefactors_Medium2016MuonBCDEF", "FsrModeling", true , "Systematic from FSR modeling"       , 1200,500, 0,0, 10,120, false);  
  computeSystematics("2018-02-09/trackToMediumMuonTightIso/sfMediumMuonBtoFNominal.root", "2018-02-09/trackToMediumMuonTightIso/sfMediumMuonBtoFfsrUp.root"  , "scalefactors_Medium2016MuonBCDEF", "FsrModeling", false, "Systematic from FSR modeling"       , 1200,500, 0,0, 10,120, false);  
  computeSystematics("2018-02-09/trackToMediumMuonTightIso/sfMediumMuonBtoFNominal.root", "2018-02-09/trackToMediumMuonTightIso/sfMediumMuonBtoFResFunc.root", "scalefactors_Medium2016MuonBCDEF", "ResFunc"    , true , "Systematic from resolution modeling", 1200,500, 0,0, 10,120, false);  
  computeSystematics("2018-02-09/trackToMediumMuonTightIso/sfMediumMuonGtoHNominal.root", "2018-02-09/trackToMediumMuonTightIso/sfMediumMuonGtoHAltBkg.root" , "scalefactors_Medium2016MuonGH", "AltBkg"     , true , "Systematic from background modeling", 1200,500, 0,0, 10,120, false);  
  computeSystematics("2018-02-09/trackToMediumMuonTightIso/sfMediumMuonGtoHNominal.root", "2018-02-09/trackToMediumMuonTightIso/sfMediumMuonGtoHAltGen.root" , "scalefactors_Medium2016MuonGH", "AltGen"     , true , "Systematic from generator choice"   , 1200,500, 0,0, 10,120, false);  
  computeSystematics("2018-02-09/trackToMediumMuonTightIso/sfMediumMuonGtoHNominal.root", "2018-02-09/trackToMediumMuonTightIso/sfMediumMuonGtoHAltTag.root" , "scalefactors_Medium2016MuonGH", "AltTag"     , true , "Systematic from tag selection bias" , 1200,500, 0,0, 10,120, false);  
  computeSystematics("2018-02-09/trackToMediumMuonTightIso/sfMediumMuonGtoHNominal.root", "2018-02-09/trackToMediumMuonTightIso/sfMediumMuonGtoHfsrDown.root", "scalefactors_Medium2016MuonGH", "FsrModeling", true , "Systematic from FSR modeling"       , 1200,500, 0,0, 10,120, false);  
  computeSystematics("2018-02-09/trackToMediumMuonTightIso/sfMediumMuonGtoHNominal.root", "2018-02-09/trackToMediumMuonTightIso/sfMediumMuonGtoHfsrUp.root"  , "scalefactors_Medium2016MuonGH", "FsrModeling", false, "Systematic from FSR modeling"       , 1200,500, 0,0, 10,120, false);  
  computeSystematics("2018-02-09/trackToMediumMuonTightIso/sfMediumMuonGtoHNominal.root", "2018-02-09/trackToMediumMuonTightIso/sfMediumMuonGtoHResFunc.root", "scalefactors_Medium2016MuonGH", "ResFunc"    , true , "Systematic from resolution modeling", 1200,500, 0,0, 10,120, false);  
}
void computeAllEle() {
  computeSystematics("2018-02-27/gsfElectronToMedium/sfMediumElectronBtoFNominal.root", "2018-02-27/gsfElectronToMedium/sfMediumElectronBtoFAltBkg.root" , "scalefactors_MediumElectronBCDEF", "AltBkg"     , true , "Systematic from background modeling", 1200,500, 0,0, 10,100, false);  
  computeSystematics("2018-02-27/gsfElectronToMedium/sfMediumElectronBtoFNominal.root", "2018-02-27/gsfElectronToMedium/sfMediumElectronBtoFAltGen.root" , "scalefactors_MediumElectronBCDEF", "AltGen"     , true , "Systematic from generator choice"   , 1200,500, 0,0, 10,100, false);  
  computeSystematics("2018-02-27/gsfElectronToMedium/sfMediumElectronBtoFNominal.root", "2018-02-27/gsfElectronToMedium/sfMediumElectronBtoFAltTag.root" , "scalefactors_MediumElectronBCDEF", "AltTag"     , true , "Systematic from tag selection bias" , 1200,500, 0,0, 10,100, false);  
  computeSystematics("2018-02-27/gsfElectronToMedium/sfMediumElectronBtoFNominal.root", "2018-02-27/gsfElectronToMedium/sfMediumElectronBtoFfsrDown.root", "scalefactors_MediumElectronBCDEF", "FsrModeling", true , "Systematic from FSR modeling"       , 1200,500, 0,0, 10,100, false);  
  computeSystematics("2018-02-27/gsfElectronToMedium/sfMediumElectronBtoFNominal.root", "2018-02-27/gsfElectronToMedium/sfMediumElectronBtoFfsrUp.root"  , "scalefactors_MediumElectronBCDEF", "FsrModeling", false, "Systematic from FSR modeling"       , 1200,500, 0,0, 10,100, false);  
  computeSystematics("2018-02-27/gsfElectronToMedium/sfMediumElectronBtoFNominal.root", "2018-02-27/gsfElectronToMedium/sfMediumElectronBtoFResFunc.root", "scalefactors_MediumElectronBCDEF", "ResFunc"    , true , "Systematic from resolution modeling", 1200,500, 0,0, 10,100, false);  
  computeSystematics("2018-02-27/gsfElectronToMedium/sfMediumElectronGtoHNominal.root", "2018-02-27/gsfElectronToMedium/sfMediumElectronGtoHAltBkg.root" , "scalefactors_MediumElectronGH", "AltBkg"     , true , "Systematic from background modeling", 1200,500, 0,0, 10,100, false);  
  computeSystematics("2018-02-27/gsfElectronToMedium/sfMediumElectronGtoHNominal.root", "2018-02-27/gsfElectronToMedium/sfMediumElectronGtoHAltGen.root" , "scalefactors_MediumElectronGH", "AltGen"     , true , "Systematic from generator choice"   , 1200,500, 0,0, 10,100, false);  
  computeSystematics("2018-02-27/gsfElectronToMedium/sfMediumElectronGtoHNominal.root", "2018-02-27/gsfElectronToMedium/sfMediumElectronGtoHAltTag.root" , "scalefactors_MediumElectronGH", "AltTag"     , true , "Systematic from tag selection bias" , 1200,500, 0,0, 10,100, false);  
  computeSystematics("2018-02-27/gsfElectronToMedium/sfMediumElectronGtoHNominal.root", "2018-02-27/gsfElectronToMedium/sfMediumElectronGtoHfsrDown.root", "scalefactors_MediumElectronGH", "FsrModeling", true , "Systematic from FSR modeling"       , 1200,500, 0,0, 10,100, false);  
  computeSystematics("2018-02-27/gsfElectronToMedium/sfMediumElectronGtoHNominal.root", "2018-02-27/gsfElectronToMedium/sfMediumElectronGtoHfsrUp.root"  , "scalefactors_MediumElectronGH", "FsrModeling", false, "Systematic from FSR modeling"       , 1200,500, 0,0, 10,100, false);  
  computeSystematics("2018-02-27/gsfElectronToMedium/sfMediumElectronGtoHNominal.root", "2018-02-27/gsfElectronToMedium/sfMediumElectronGtoHResFunc.root", "scalefactors_MediumElectronGH", "ResFunc"    , true , "Systematic from resolution modeling", 1200,500, 0,0, 10,100, false);  
}
/*  TString nominalFileName,     // root file with nominal scale factors
  TString altFileName,         // root file with the alternative scale factors
  TString sfName,              // name of the scale factors histo in the root file
  TString systName,            // name of the systematic source
  bool replaceExisting=true,   // whether to replace the existing calculation (if not, take the max)
  TString plotTitle="",
  int width=2400,
  int height=1200,
  double xmin=0,
  double xmax=0,
  double ymin=0,
  double ymax=0
*/
