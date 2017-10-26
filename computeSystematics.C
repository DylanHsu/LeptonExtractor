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
  computeSystematics("2017-10-19/electrons/bkg_lpi_emu/scalefactors_ele_medium_2016.root", "2017-10-19/electrons/sig_fsrUp/scalefactors_ele_medium_2016.root", "scalefactors_Medium_Electron", "signalFsrTNP", true, "Absolute SF uncertainty from FSR modeling (Medium Electron selection)", 2400, 1200, 0, 0, 10, 40, false);
  computeSystematics("2017-10-19/electrons/bkg_lpi_emu/scalefactors_ele_medium_2016.root", "2017-10-19/electrons/sig_fsrDown/scalefactors_ele_medium_2016.root", "scalefactors_Medium_Electron", "signalFsrTNP", false, "Absolute SF uncertainty from FSR modeling (Medium Electron selection)", 2400, 1200, 0, 0, 10, 40, false);
  computeSystematics("2017-10-19/electrons/bkg_lpi_emu/scalefactors_ele_medium_2016.root", "2017-10-19/electrons/sig_fsrUp/scalefactors_ele_medium_2016.root", "scalefactors_Medium_Electron", "signalFsrTNP", false, "Absolute SF uncertainty from FSR modeling (Medium Electron selection)", 2400, 1200, 0, 0, 40, 100, true);
  //computeSystematics("2017-10-19/electrons/bkg_lpi_emu/scalefactors_ele_medium_2016.root", "2017-10-19/electrons/sig_fsrDown/scalefactors_ele_medium_2016.root", "scalefactors_Medium_Electron", "signalFsrTNP", false, "Absolute SF uncertainty from FSR modeling (Medium Electron selection)", 2400, 1200, 0, 0, 40, 100, true);
  computeSystematics("2017-10-19/electrons/bkg_lpi_emu/scalefactors_ele_medium_2016.root", "2017-10-19/electrons/sig_genmc/scalefactors_ele_medium_2016.root", "scalefactors_Medium_Electron", "signalResTNP", true, "Absolute SF uncertainty from resolution modeling (Medium Electron selection)", 2400, 1200, 0, 0, 10, 40, false);
  computeSystematics("2017-10-19/electrons/bkg_lpi_emu/scalefactors_ele_medium_2016.root", "2017-10-19/electrons/sig_genmc/scalefactors_ele_medium_2016.root", "scalefactors_Medium_Electron", "signalResTNP", true, "Absolute SF uncertainty from resolution modeling (Medium Electron selection)", 2400, 1200, 0, 0, 40, 100, true);
  computeSystematics("2017-10-19/electrons/bkg_lpi_emu/scalefactors_ele_medium_2016.root", "2017-10-19/electrons/nominal/scalefactors_ele_medium_2016.root", "scalefactors_Medium_Electron", "bkgModelTNP", true, "Absolute SF uncertainty from background modeling (Medium Electron selection)", 2400, 1200, 0, 0, 10, 40, false);
  computeSystematics("2017-10-19/electrons/bkg_lpi_emu/scalefactors_ele_medium_2016.root", "2017-10-19/electrons/nominal/scalefactors_ele_medium_2016.root", "scalefactors_Medium_Electron", "bkgModelTNP", true, "Absolute SF uncertainty from background modeling (Medium Electron selection)", 2400, 1200, 0, 0, 40, 100, true);
  computeSystematics("2017-10-19/electrons/bkg_lpi_emu/scalefactors_ele_medium_2016.root", "2017-10-19/electrons/alt_tag/scalefactors_ele_medium_2016.root", "scalefactors_Medium_Electron", "tagBiasTNP", true, "Absolute SF uncertainty from tag selection (Medium Electron selection)", 2400, 1200, 0, 0, 10, 40, false);
  computeSystematics("2017-10-19/electrons/bkg_lpi_emu/scalefactors_ele_medium_2016.root", "2017-10-19/electrons/alt_tag/scalefactors_ele_medium_2016.root", "scalefactors_Medium_Electron", "tagBiasTNP", true, "Absolute SF uncertainty from tag selection (Medium Electron selection)", 2400, 1200, 0, 0, 40, 100, true);
  computeSystematics("2017-10-19/electrons/bkg_lpi_emu/scalefactors_ele_medium_2016.root", "2017-10-19/electrons/mc_lo/scalefactors_ele_medium_2016_lo.root", "scalefactors_Medium_Electron", "generatorChoiceTNP", true, "Absolute SF unc. from MC generator (Medium Electron selection)", 2400, 1200, 0, 0, 10, 40, false);
  computeSystematics("2017-10-19/electrons/bkg_lpi_emu/scalefactors_ele_medium_2016.root", "2017-10-19/electrons/mc_lo/scalefactors_ele_medium_2016_lo.root", "scalefactors_Medium_Electron", "generatorChoiceTNP", true, "Absolute SF unc. from MC generator (Medium Electron selection)", 2400, 1200, 0, 0, 40, 100, true);

  computeSystematics("2017-10-19/muons/bkg_lpi/scalefactors_mu_medium_2016.root", "2017-10-19/muons/sig_fsrUp/scalefactors_mu_medium_2016.root", "scalefactors_Medium_Muon", "signalFsrTNP", true, "Absolute SF unc. from FSR modeling (Medium Muon selection)", 2400, 1200, 0, 0, 10, 40, false);
  computeSystematics("2017-10-19/muons/bkg_lpi/scalefactors_mu_medium_2016.root", "2017-10-19/muons/sig_fsrDown/scalefactors_mu_medium_2016.root", "scalefactors_Medium_Muon", "signalFsrTNP", false, "Absolute SF unc. from FSR modeling (Medium Muon selection)", 2400, 1200, 0, 0, 10, 40, false);
  //computeSystematics("2017-10-19/muons/bkg_lpi/scalefactors_mu_medium_2016.root", "2017-10-19/muons/sig_fsrUp/scalefactors_mu_medium_2016.root", "scalefactors_Medium_Muon", "signalFsrTNP", true, "Absolute SF unc. from FSR modeling (Medium Muon selection)", 2400, 1200, 0, 0, 40, 50, true);
  //computeSystematics("2017-10-19/muons/bkg_lpi/scalefactors_mu_medium_2016.root", "2017-10-19/muons/sig_fsrUp/scalefactors_mu_medium_2016.root", "scalefactors_Medium_Muon", "signalFsrTNP", true, "Absolute SF unc. from FSR modeling (Medium Muon selection)", 2400, 1200, 0, 0, 50, 100, true);
  computeSystematics("2017-10-19/muons/bkg_lpi/scalefactors_mu_medium_2016.root", "2017-10-19/muons/sig_fsrDown/scalefactors_mu_medium_2016.root", "scalefactors_Medium_Muon", "signalFsrTNP", false, "Absolute SF unc. from FSR modeling (Medium Muon selection)", 2400, 1200, 0, 0, 40, 50, true);
  computeSystematics("2017-10-19/muons/bkg_lpi/scalefactors_mu_medium_2016.root", "2017-10-19/muons/sig_fsrDown/scalefactors_mu_medium_2016.root", "scalefactors_Medium_Muon", "signalFsrTNP", false, "Absolute SF unc. from FSR modeling (Medium Muon selection)", 2400, 1200, 0, 0, 50, 100, true);
  computeSystematics("2017-10-19/muons/bkg_lpi/scalefactors_mu_medium_2016.root", "2017-10-19/muons/sig_genmc/scalefactors_mu_medium_2016.root", "scalefactors_Medium_Muon", "signalResTNP", true, "Absolute SF unc. from resolution modeling (Medium Muon selection)", 2400, 1200, 0, 0, 10, 40, false);
  computeSystematics("2017-10-19/muons/bkg_lpi/scalefactors_mu_medium_2016.root", "2017-10-19/muons/sig_genmc/scalefactors_mu_medium_2016.root", "scalefactors_Medium_Muon", "signalResTNP", true, "Absolute SF unc. from resolution modeling (Medium Muon selection)", 2400, 1200, 0, 0, 40, 50, false);
  computeSystematics("2017-10-19/muons/bkg_lpi/scalefactors_mu_medium_2016.root", "2017-10-19/muons/sig_genmc/scalefactors_mu_medium_2016.root", "scalefactors_Medium_Muon", "signalResTNP", true, "Absolute SF unc. from resolution modeling (Medium Muon selection)", 2400, 1200, 0, 0, 50, 100, true);
  computeSystematics("2017-10-19/muons/bkg_lpi/scalefactors_mu_medium_2016.root", "2017-10-19/muons/nominal/scalefactors_mu_medium_2016.root", "scalefactors_Medium_Muon", "bkgModelTNP", true, "Absolute SF unc. from background modeling (Medium Muon selection)", 2400, 1200, 0, 0, 10, 40, false);
  computeSystematics("2017-10-19/muons/bkg_lpi/scalefactors_mu_medium_2016.root", "2017-10-19/muons/nominal/scalefactors_mu_medium_2016.root", "scalefactors_Medium_Muon", "bkgModelTNP", true, "Absolute SF unc. from background modeling (Medium Muon selection)", 2400, 1200, 0, 0, 40, 50, false);
  computeSystematics("2017-10-19/muons/bkg_lpi/scalefactors_mu_medium_2016.root", "2017-10-19/muons/nominal/scalefactors_mu_medium_2016.root", "scalefactors_Medium_Muon", "bkgModelTNP", true, "Absolute SF unc. from background modeling (Medium Muon selection)", 2400, 1200, 0, 0, 50, 100, true);
  computeSystematics("2017-10-19/muons/bkg_lpi/scalefactors_mu_medium_2016.root", "2017-10-19/muons/alt_tag/scalefactors_mu_medium_2016.root", "scalefactors_Medium_Muon", "tagBiasTNP", true, "Absolute SF unc. from tag selection (Medium Muon selection)", 2400, 1200, 0, 0, 10, 40, false);
  computeSystematics("2017-10-19/muons/bkg_lpi/scalefactors_mu_medium_2016.root", "2017-10-19/muons/alt_tag/scalefactors_mu_medium_2016.root", "scalefactors_Medium_Muon", "tagBiasTNP", true, "Absolute SF unc. from tag selection (Medium Muon selection)", 2400, 1200, 0, 0, 40, 50, false);
  computeSystematics("2017-10-19/muons/bkg_lpi/scalefactors_mu_medium_2016.root", "2017-10-19/muons/alt_tag/scalefactors_mu_medium_2016.root", "scalefactors_Medium_Muon", "tagBiasTNP", true, "Absolute SF unc. from tag selection (Medium Muon selection)", 2400, 1200, 0, 0, 50, 100, true);
  computeSystematics("2017-10-19/muons/bkg_lpi/scalefactors_mu_medium_2016.root", "2017-10-19/muons/mc_lo/scalefactors_mu_medium_2016_lo.root", "scalefactors_Medium_Muon", "generatorChoiceTNP", true, "Absolute SF unc. from MC generator (Medium Muon selection)", 2400, 1200, 0, 0, 10, 40, false);
  computeSystematics("2017-10-19/muons/bkg_lpi/scalefactors_mu_medium_2016.root", "2017-10-19/muons/mc_lo/scalefactors_mu_medium_2016_lo.root", "scalefactors_Medium_Muon", "generatorChoiceTNP", true, "Absolute SF unc. from MC generator (Medium Muon selection)", 2400, 1200, 0, 0, 40, 50, false);
  computeSystematics("2017-10-19/muons/bkg_lpi/scalefactors_mu_medium_2016.root", "2017-10-19/muons/mc_lo/scalefactors_mu_medium_2016_lo.root", "scalefactors_Medium_Muon", "generatorChoiceTNP", true, "Absolute SF unc. from MC generator (Medium Muon selection)", 2400, 1200, 0, 0, 50, 100, true);

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
