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
#include <TRandom3.h>
#include <iostream>
Int_t mit_red  = 1861; 
Int_t mit_gray = 1862; 
//int marker_colors[] = {kBlack, kRed, kOrange+10, kOrange, kSpring+9, kGreen+3, kCyan-2, kBlue, kViolet+8, kMagenta+3};
int marker_colors[] = {kBlack, mit_red, mit_gray, kBlue, kMagenta+3};
int marker_styles[] = {20, 21, 22, 23, 20, 21, 22, 23, 20, 21};
Float_t ele_pt_bins[] = {10,20,40,80,200};
Float_t ele_eta_bins[] = {-2.5,-2.0,-1.566,-1.4442,-0.8,0,0.8,1.4442,1.566,2.0,2.5};
Int_t n_ele_pt_bins=4;
Int_t n_ele_eta_bins=10;
Float_t mu_pt_bins[] = {10,15,20,25,30,40,50,100,200,1000};
Float_t mu_eta_bins[] = {-2.4,-2.1,-1.2,-0.9,0,0.9,1.2,2.1,2.4};
Int_t n_mu_pt_bins=9;
Int_t n_mu_eta_bins=8;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void scale_factors(TString dataEffFile, TString mcEffFile, TString outputFileName, TString outputBasename, TString selectionName, TString outputDir, int width=2400, int height=1200, double xmin=0, double xmax=0, double ymin=0, double ymax=0){ 
  TFile *outputFile = new TFile(outputDir+"/"+outputFileName, "RECREATE");
  {
    TFile *f_data = TFile::Open(dataEffFile.Data(), "READ");
    TH2D *hEffData     = (TH2D*) f_data->Get("hEffEtaPt");  hEffData     ->SetDirectory(0); 
    TH2D *hErrorLoData = (TH2D*) f_data->Get("hErrlEtaPt"); hErrorLoData->SetDirectory(0); 
    TH2D *hErrorHiData = (TH2D*) f_data->Get("hErrhEtaPt"); hErrorHiData->SetDirectory(0); 
    hEffData->SetName(TString("eff_data_"+selectionName).Data());
    hEffData->SetTitle(TString("Efficiency for "+selectionName+" selection (Data)").Data());
    
    TFile *f_mc   = TFile::Open(mcEffFile.Data(), "READ");
    TH2D *hEffMC      = (TH2D*) f_mc->Get("hEffEtaPt");  hEffMC     ->SetDirectory(0); 
    TH2D *hErrorLoMC = (TH2D*) f_mc->Get("hErrlEtaPt"); hErrorLoMC->SetDirectory(0); 
    TH2D *hErrorHiMC = (TH2D*) f_mc->Get("hErrhEtaPt"); hErrorHiMC->SetDirectory(0); 
    hEffMC->SetName(TString("eff_mc_"+selectionName).Data());
    hEffMC->SetTitle(TString("Efficiency for "+selectionName+" selection (MC)").Data());
    
    // Divide Data/MC to get the scale factors
    TH2D *hSF = (TH2D*) hEffData->Clone(); hSF->SetDirectory(0); hSF->Clear(); hSF->Reset(); 
    hSF->SetName(TString("scalefactors_"+selectionName).Data()); 
    hSF->SetTitle(TString("Scale factors for "+selectionName+" selection").Data());

    // Propagate asymmetric statistical errors
    TH2D *hSFErrorLo = (TH2D*) hSF->Clone(); hSFErrorLo->SetDirectory(0); hSFErrorLo->Clear(); hSFErrorLo->Reset();
    TH2D *hSFErrorHi = (TH2D*) hSF->Clone(); hSFErrorHi->SetDirectory(0); hSFErrorHi->Clear(); hSFErrorHi->Reset();
    hSFErrorLo->SetName(TString("scalefactors_"+selectionName+"_statErrorLow" ).Data());
    hSFErrorHi->SetName(TString("scalefactors_"+selectionName+"_statErrorHigh").Data());
    hSFErrorLo->SetTitle(TString("-1#sigma quantile stat. unc. for "+selectionName+" scale factors").Data());
    hSFErrorHi->SetTitle(TString("+1#sigma quantile stat. unc. for "+selectionName+" scale factors").Data());
    
    unsigned long int time_now = static_cast<unsigned long int>(time(NULL));
    unsigned int randomToySeed=(time_now-731178000); // random seed based on Dylan's age in seconds
    TRandom3 toymaker(randomToySeed);

    for(int i = 1; i <= hEffMC->GetNbinsX(); i++) { for(int j = 1; j <= hEffMC->GetNbinsY(); j++) {
      unsigned int nbin = hEffMC->GetBin(i,j);
      double sf=hEffData->GetBinContent(nbin) / hEffMC->GetBinContent(nbin);
      double sf_error_hi, sf_error_lo;
      bool badSf=false; if(sf!=sf || sf<=0 || sf>2) { sf=1; badSf=true;}
      hSF->SetBinContent(nbin, sf);
      if(!badSf) {
        double sfEstErrorHi = sf * sqrt( pow(hErrorHiData->GetBinContent(nbin) / hEffData->GetBinContent(nbin), 2) + pow(hErrorLoMC->GetBinContent(nbin) / hEffMC->GetBinContent(nbin),2));
        double sfEstErrorLo = sf * sqrt( pow(hErrorLoData->GetBinContent(nbin) / hEffData->GetBinContent(nbin), 2) + pow(hErrorHiMC->GetBinContent(nbin) / hEffMC->GetBinContent(nbin),2));
        TH1F *toySFs = new TH1F("toySFs","toySFs", 100, TMath::Max(0.,sf-5*sfEstErrorLo), TMath::Min(2.,sf+5*sfEstErrorHi));
        //printf("toySFs %.3f to %.3f\n", toySFs->GetBinLowEdge(1), toySFs->GetBinLowEdge(101));
        for(int iToy=0; iToy<1000; iToy++) {
          // Toys are Gaussian about 0 with sigma 1
          double toyData=toymaker.Gaus(0,1);
          double toyMC=toymaker.Gaus(0,1);
          // Use the upper or lower error bar on the signal size based on the sign of the toys
          double toyEffData = hEffData->GetBinContent(nbin) + toyData * (toyData>0? hErrorHiData->GetBinContent(nbin) : hErrorLoData->GetBinContent(nbin));
          double toyEffMC   = hEffMC  ->GetBinContent(nbin) + toyMC   * (toyMC  >0? hErrorHiMC  ->GetBinContent(nbin) : hErrorLoMC  ->GetBinContent(nbin));
          double toySF      = toyEffData / toyEffMC;
          //printf("toyEff #%d: %.3f\n", iToy, toyEff);
          toySFs->Fill(toySF);
        }
        double quantileProbs[3]={0.159,0.5,0.841};
        double theQuantiles[3];
        toySFs->GetQuantiles(3, theQuantiles, quantileProbs);
        delete toySFs;
        sf_error_lo=TMath::Max(0.,sf - TMath::Max(0.,theQuantiles[0]));
        sf_error_hi=TMath::Max(0.,TMath::Min(2.,theQuantiles[2]) - sf);
      } else {
        sf_error_hi = 1;
        sf_error_lo = 1;
      }
      hSFErrorHi->SetBinContent(nbin, sf_error_hi);
      hSFErrorLo->SetBinContent(nbin, sf_error_lo);
      // Choose the max for the s.f. and eff. histogram errors
      // Analyzers who naively use this uncertainty value will just get worse sensitivity :^)
      hSF->SetBinError(nbin, TMath::Max(
        hSFErrorHi->GetBinContent(nbin),
        hSFErrorLo->GetBinContent(nbin)
      ));
      hEffData->SetBinError(nbin, TMath::Max(
        hErrorLoData->GetBinContent(nbin),
        hErrorHiData->GetBinContent(nbin)
      ));
      hEffMC->SetBinError(nbin, TMath::Max(
        hErrorLoMC->GetBinContent(nbin),
        hErrorHiMC->GetBinContent(nbin)
      ));
    }} 

    // Start drawing stuff 
    gStyle->SetOptStat(0);
    gStyle->SetPaintTextFormat("4.3f");
    TPaletteAxis *palette_axis; TCanvas *canvas[5]; TH2D* theHist[5]; int thePalette[5];
    thePalette[0]=kBlackBody; canvas[0]=new TCanvas("c_eff_data"   ,"c_eff_data"   ,width,height); theHist[0]=hEffData  ; hEffData  ->SetMinimum(  0); hEffData   ->SetMaximum(  1);
    thePalette[1]=kBlackBody; canvas[1]=new TCanvas("c_eff_mc"     ,"c_eff_mc"     ,width,height); theHist[1]=hEffMC    ; hEffMC    ->SetMinimum(  0); hEffMC     ->SetMaximum(  1);
    thePalette[2]=kBlackBody; canvas[2]=new TCanvas("c_sf"         ,"c_sf"         ,width,height); theHist[2]=hSF       ; hSF       ->SetMinimum(0.6); hSF         ->SetMaximum(1.4);
    thePalette[3]=kCool;      canvas[3]=new TCanvas("c_sf_error_hi","c_sf_error_hi",width,height); theHist[3]=hSFErrorHi; hSFErrorHi->SetMinimum(  0); hSFErrorHi->SetMaximum(0.1);
    thePalette[4]=kCool;      canvas[4]=new TCanvas("c_sf_error_lo","c_sf_error_lo",width,height); theHist[4]=hSFErrorLo; hSFErrorLo->SetMinimum(  0); hSFErrorLo->SetMaximum(0.1);

    TLine *lines[1000]; int nLines=0;
    bool xRange=(xmin!=xmax); bool yRange=(ymin!=ymax);
    if(xRange) theHist[0]->GetXaxis()->SetRangeUser(xmin,xmax); else theHist[0]->GetXaxis()->SetRangeUser(theHist[0]->GetXaxis()->GetBinLowEdge(1), theHist[0]->GetXaxis()->GetBinLowEdge(theHist[0]->GetXaxis()->GetNbins()+1));
    if(yRange) theHist[0]->GetYaxis()->SetRangeUser(ymin,ymax); else theHist[0]->GetYaxis()->SetRangeUser(theHist[0]->GetYaxis()->GetBinLowEdge(1), theHist[0]->GetYaxis()->GetBinLowEdge(theHist[0]->GetYaxis()->GetNbins()+1));
    for(unsigned i=(xRange ? theHist[0]->GetXaxis()->FindBin(xmin+.001):1); i<(unsigned)(xRange ? theHist[0]->GetXaxis()->FindBin(xmax-.001) : theHist[0]->GetXaxis()->GetNbins()); i++) {
      lines[nLines]=new TLine(
        theHist[0]->GetXaxis()->GetBinUpEdge(i),   
        yRange? theHist[0]->GetYaxis()->GetBinUpEdge(theHist[0]->GetYaxis()->FindBin(ymin+.001)-1) : theHist[0]->GetYaxis()->GetBinUpEdge(0),
        theHist[0]->GetXaxis()->GetBinUpEdge(i),
        yRange? theHist[0]->GetYaxis()->GetBinUpEdge(theHist[0]->GetYaxis()->FindBin(ymax-.001)) : theHist[0]->GetYaxis()->GetBinUpEdge(theHist[0]->GetYaxis()->GetNbins())
      ); lines[nLines]->SetLineStyle(kDashed); nLines++;
    }
    for(unsigned j=(yRange ? theHist[0]->GetYaxis()->FindBin(ymin+.001):1); j<(unsigned)(yRange ? theHist[0]->GetYaxis()->FindBin(ymax-.001) : theHist[0]->GetYaxis()->GetNbins()); j++) {
      lines[nLines]=new TLine(
        xRange? theHist[0]->GetXaxis()->GetBinUpEdge(theHist[0]->GetXaxis()->FindBin(xmin+.001)-1) : theHist[0]->GetXaxis()->GetBinUpEdge(0),
        theHist[0]->GetYaxis()->GetBinUpEdge(j),   
        xRange? theHist[0]->GetXaxis()->GetBinUpEdge(theHist[0]->GetXaxis()->FindBin(xmax-.001)) : theHist[0]->GetXaxis()->GetBinUpEdge(theHist[0]->GetXaxis()->GetNbins()),
        theHist[0]->GetYaxis()->GetBinUpEdge(j)
      ); lines[nLines]->SetLineStyle(kDashed); nLines++;
    }
   
    for(unsigned i=0;i<5;i++) {
      canvas[i]->cd(); gStyle->SetPalette(thePalette[i]);
      //if(flavor=="Muon") gStyle->SetPaintTextFormat("4.2f"); else gStyle->SetPaintTextFormat("4.2f");
      //if(flavor=="Muon") theHist[i]->Draw("TEXT45 COLZ"); else theHist[i]->Draw("TEXT45 COLZ");
      gStyle->SetPaintTextFormat("4.2f");
      theHist[i]->Draw("TEXT45 COLZ");
      //canvas[i]->SetLogy();
      canvas[i]->Update();
      theHist[i]->GetXaxis()->SetTitle("#eta");
      theHist[i]->GetXaxis()->SetTitleOffset(1.1);
      theHist[i]->GetXaxis()->SetTitleSize(0.04);
      theHist[i]->GetXaxis()->SetLabelSize(0.04);
      theHist[i]->GetYaxis()->SetTitle("p_{T} [GeV]");
      theHist[i]->GetYaxis()->SetTitleOffset(1.1);
      theHist[i]->GetYaxis()->SetTitleSize(0.04);
      theHist[i]->GetYaxis()->SetLabelSize(0.04);
      if(xRange) theHist[i]->GetXaxis()->SetRangeUser(xmin,xmax); else theHist[i]->GetXaxis()->SetRangeUser(theHist[i]->GetXaxis()->GetBinLowEdge(1), theHist[i]->GetXaxis()->GetBinLowEdge(theHist[i]->GetXaxis()->GetNbins()+1));
      if(yRange) theHist[i]->GetYaxis()->SetRangeUser(ymin,ymax); else theHist[i]->GetYaxis()->SetRangeUser(theHist[i]->GetYaxis()->GetBinLowEdge(1), theHist[i]->GetYaxis()->GetBinLowEdge(theHist[i]->GetYaxis()->GetNbins()+1));
      theHist[i]->GetYaxis()->SetMoreLogLabels();
      theHist[i]->SetMarkerSize(1.1);
      palette_axis = (TPaletteAxis*) theHist[i]->GetListOfFunctions()->FindObject("palette"); 
      palette_axis->SetLabelSize(0.02);
      canvas[i]->Update();
      for(int nl=0;nl<nLines;nl++) lines[nl]->Draw("SAME");
      canvas[i]->Print(Form("%s%s_%.1f-%.1f_%.1f-%.1f.png",outputDir.Data(),theHist[i]->GetName(),xmin,xmax,ymin,ymax));
      canvas[i]->Print(Form("%s%s_%.1f-%.1f_%.1f-%.1f.pdf",outputDir.Data(),theHist[i]->GetName(),xmin,xmax,ymin,ymax));
    }

    for(unsigned i=0;i<5;i++) {
      theHist[i]->GetXaxis()->SetRangeUser(theHist[i]->GetXaxis()->GetBinLowEdge(1), theHist[i]->GetXaxis()->GetBinLowEdge(theHist[i]->GetXaxis()->GetNbins()+1));
      theHist[i]->GetYaxis()->SetRangeUser(theHist[i]->GetYaxis()->GetBinLowEdge(1), theHist[i]->GetYaxis()->GetBinLowEdge(theHist[i]->GetYaxis()->GetNbins()+1));
    }
    // Write to file
    outputFile->cd();
    hEffData->Write();
    hEffMC->Write();
    hSF->Write();
    hSFErrorLo->Write();
    hSFErrorHi->Write();
    
    
    f_data->Close();
    f_mc->Close();
  }
  outputFile->Close();
  printf("Saved efficiencies, scale factors, and statistical errors in %s\n", outputFileName.Data());
}

void doAll() {
  scale_factors("2018-02-09/trackToMediumMuonTightIso/tnp_BtoF_sigTemplate_bkgDasPlusExp/eff.root"      , "2018-02-09/trackToMediumMuonTightIso/tnp_BtoF_mc_NLOMC/eff.root"        , "sfMediumMuonBtoFAltBkg.root" , "Medium2016MuonBtoFAltBkg" , "Medium2016MuonBCDEF", "2018-02-09/trackToMediumMuonTightIso/", 2200,800, 0,0, 10,120);  
  //scale_factors("2018-02-09/trackToMediumMuonTightIso/tnp_BtoF_sigTemplate_bkgTemplate_LOMC/eff.root"   , "2018-02-09/trackToMediumMuonTightIso/tnp_BtoF_mc_LOMC/eff.root"         , "sfMediumMuonBtoFAltGen.root" , "Medium2016MuonBtoFAltGen" , "Medium2016MuonBCDEF", "2018-02-09/trackToMediumMuonTightIso/", 2200,800, 0,0, 10,120);  
  scale_factors("2018-02-09/trackToMediumMuonTightIso/tnp_BtoF_sigTemplate_bkgTemplate_LOMC/eff.root"   , "2018-02-09/trackToMediumMuonTightIso/tnp_BtoF_mc_NLOMC/eff.root"         , "sfMediumMuonBtoFAltGen.root" , "Medium2016MuonBtoFAltGen" , "Medium2016MuonBCDEF", "2018-02-09/trackToMediumMuonTightIso/", 2200,800, 0,0, 10,120);  
  scale_factors("2018-02-09/trackToMediumMuonTightIso/tnp_BtoF_sigTemplate_bkgTemplate_altTag/eff.root" , "2018-02-09/trackToMediumMuonTightIso/tnp_BtoF_mc_NLOMC_altTag/eff.root" , "sfMediumMuonBtoFAltTag.root" , "Medium2016MuonBtoFAltTag" , "Medium2016MuonBCDEF", "2018-02-09/trackToMediumMuonTightIso/", 2200,800, 0,0, 10,120);  
  scale_factors("2018-02-09/trackToMediumMuonTightIso/tnp_BtoF_sigTemplate_bkgTemplate_resFunc/eff.root", "2018-02-09/trackToMediumMuonTightIso/tnp_BtoF_mc_NLOMC/eff.root"        , "sfMediumMuonBtoFResFunc.root", "Medium2016MuonBtoFResFunc", "Medium2016MuonBCDEF", "2018-02-09/trackToMediumMuonTightIso/", 2200,800, 0,0, 10,120);  
  scale_factors("2018-02-09/trackToMediumMuonTightIso/tnp_GtoH_sigTemplate_bkgDasPlusExp/eff.root"      , "2018-02-09/trackToMediumMuonTightIso/tnp_GtoH_mc_NLOMC/eff.root"        , "sfMediumMuonGtoHAltBkg.root" , "Medium2016MuonGtoHAltBkg" , "Medium2016MuonGH", "2018-02-09/trackToMediumMuonTightIso/", 2200,800, 0,0, 10,120);  
  scale_factors("2018-02-09/trackToMediumMuonTightIso/tnp_BtoF_sigTemplate_bkgTemplate_powhegPythia/eff.root", "2018-02-09/trackToMediumMuonTightIso/tnp_BtoF_mc_NLOMC_powhegPythia/eff.root", "sfMediumMuonBtoFpowhegPythia.root", "Medium2016MuonBtoFpowhegPythia", "Medium2016MuonBCDEF", "2018-02-09/trackToMediumMuonTightIso/", 2200,800, 0,0, 10,120);  
  scale_factors("2018-02-09/trackToMediumMuonTightIso/tnp_BtoF_sigTemplate_bkgTemplate_powhegPhotos/eff.root", "2018-02-09/trackToMediumMuonTightIso/tnp_BtoF_mc_NLOMC_powhegPhotos/eff.root", "sfMediumMuonBtoFpowhegPhotos.root", "Medium2016MuonBtoFpowhegPhotos", "Medium2016MuonBCDEF", "2018-02-09/trackToMediumMuonTightIso/", 2200,800, 0,0, 10,120);  
  //scale_factors("2018-02-09/trackToMediumMuonTightIso/tnp_BtoF_sigTemplate_bkgTemplate_fsrDown/eff.root", "2018-02-09/trackToMediumMuonTightIso/tnp_BtoF_mc_NLOMC_fsrDown/eff.root", "sfMediumMuonBtoFfsrDown.root", "Medium2016MuonBtoFfsrDown", "Medium2016MuonBCDEF", "2018-02-09/trackToMediumMuonTightIso/", 2200,800, 0,0, 10,120);  
  //scale_factors("2018-02-09/trackToMediumMuonTightIso/tnp_BtoF_sigTemplate_bkgTemplate_fsrUp/eff.root"  , "2018-02-09/trackToMediumMuonTightIso/tnp_BtoF_mc_NLOMC_fsrUp/eff.root"  , "sfMediumMuonBtoFfsrUp.root"  , "Medium2016MuonBtoFfsrUp"  , "Medium2016MuonBCDEF", "2018-02-09/trackToMediumMuonTightIso/", 2200,800, 0,0, 10,120);  
  scale_factors("2018-02-09/trackToMediumMuonTightIso/tnp_BtoF_sigTemplate_bkgTemplate/eff.root"        , "2018-02-09/trackToMediumMuonTightIso/tnp_BtoF_mc_NLOMC/eff.root"        , "sfMediumMuonBtoFNominal.root", "Medium2016MuonBtoFNominal", "Medium2016MuonBCDEF", "2018-02-09/trackToMediumMuonTightIso/", 2200,800, 0,0, 10,120);  
  //scale_factors("2018-02-09/trackToMediumMuonTightIso/tnp_GtoH_sigTemplate_bkgTemplate_LOMC/eff.root"   , "2018-02-09/trackToMediumMuonTightIso/tnp_GtoH_mc_LOMC/eff.root"         , "sfMediumMuonGtoHAltGen.root" , "Medium2016MuonGtoHAltGen" , "Medium2016MuonGH", "2018-02-09/trackToMediumMuonTightIso/", 2200,800, 0,0, 10,120);  
  scale_factors("2018-02-09/trackToMediumMuonTightIso/tnp_GtoH_sigTemplate_bkgTemplate_LOMC/eff.root"   , "2018-02-09/trackToMediumMuonTightIso/tnp_GtoH_mc_NLOMC/eff.root"         , "sfMediumMuonGtoHAltGen.root" , "Medium2016MuonGtoHAltGen" , "Medium2016MuonGH", "2018-02-09/trackToMediumMuonTightIso/", 2200,800, 0,0, 10,120);  
  scale_factors("2018-02-09/trackToMediumMuonTightIso/tnp_GtoH_sigTemplate_bkgTemplate_altTag/eff.root" , "2018-02-09/trackToMediumMuonTightIso/tnp_GtoH_mc_NLOMC_altTag/eff.root" , "sfMediumMuonGtoHAltTag.root" , "Medium2016MuonGtoHAltTag" , "Medium2016MuonGH", "2018-02-09/trackToMediumMuonTightIso/", 2200,800, 0,0, 10,120);  
  scale_factors("2018-02-09/trackToMediumMuonTightIso/tnp_GtoH_sigTemplate_bkgTemplate_resFunc/eff.root", "2018-02-09/trackToMediumMuonTightIso/tnp_GtoH_mc_NLOMC/eff.root"        , "sfMediumMuonGtoHResFunc.root", "Medium2016MuonGtoHResFunc", "Medium2016MuonGH", "2018-02-09/trackToMediumMuonTightIso/", 2200,800, 0,0, 10,120);  
  scale_factors("2018-02-09/trackToMediumMuonTightIso/tnp_GtoH_sigTemplate_bkgTemplate_powhegPythia/eff.root", "2018-02-09/trackToMediumMuonTightIso/tnp_GtoH_mc_NLOMC_powhegPythia/eff.root", "sfMediumMuonGtoHpowhegPythia.root", "Medium2016MuonGtoHpowhegPythia", "Medium2016MuonGH", "2018-02-09/trackToMediumMuonTightIso/", 2200,800, 0,0, 10,120);  
  scale_factors("2018-02-09/trackToMediumMuonTightIso/tnp_GtoH_sigTemplate_bkgTemplate_powhegPhotos/eff.root", "2018-02-09/trackToMediumMuonTightIso/tnp_GtoH_mc_NLOMC_powhegPhotos/eff.root", "sfMediumMuonGtoHpowhegPhotos.root", "Medium2016MuonGtoHpowhegPhotos", "Medium2016MuonGH", "2018-02-09/trackToMediumMuonTightIso/", 2200,800, 0,0, 10,120);  
  //scale_factors("2018-02-09/trackToMediumMuonTightIso/tnp_GtoH_sigTemplate_bkgTemplate_fsrDown/eff.root", "2018-02-09/trackToMediumMuonTightIso/tnp_GtoH_mc_NLOMC_fsrDown/eff.root", "sfMediumMuonGtoHfsrDown.root", "Medium2016MuonGtoHfsrDown", "Medium2016MuonGH", "2018-02-09/trackToMediumMuonTightIso/", 2200,800, 0,0, 10,120);  
  //scale_factors("2018-02-09/trackToMediumMuonTightIso/tnp_GtoH_sigTemplate_bkgTemplate_fsrUp/eff.root"  , "2018-02-09/trackToMediumMuonTightIso/tnp_GtoH_mc_NLOMC_fsrUp/eff.root"  , "sfMediumMuonGtoHfsrUp.root"  , "Medium2016MuonGtoHfsrUp"  , "Medium2016MuonGH", "2018-02-09/trackToMediumMuonTightIso/", 2200,800, 0,0, 10,120);  
  scale_factors("2018-02-09/trackToMediumMuonTightIso/tnp_GtoH_sigTemplate_bkgTemplate/eff.root"        , "2018-02-09/trackToMediumMuonTightIso/tnp_GtoH_mc_NLOMC/eff.root"        , "sfMediumMuonGtoHNominal.root", "Medium2016MuonGtoHNominal", "Medium2016MuonGH", "2018-02-09/trackToMediumMuonTightIso/", 2200,800, 0,0, 10,120);  
}
void doAllEle() {
  scale_factors("2018-02-27/gsfElectronToMedium/tnp_BtoF_sigTemplate_bkgDasPlusExp/eff.root"      , "2018-02-27/gsfElectronToMedium/tnp_BtoF_mc_NLOMC/eff.root"        , "sfMediumElectronBtoFAltBkg.root" , "MediumElectronBtoFAltBkg" , "MediumElectronBCDEF", "2018-02-27/gsfElectronToMedium/", 2200,800, 0,0, 10,100);  
  //scale_factors("2018-02-27/gsfElectronToMedium/tnp_BtoF_sigTemplate_bkgTemplate_LOMC/eff.root"   , "2018-02-27/gsfElectronToMedium/tnp_BtoF_mc_LOMC/eff.root"         , "sfMediumElectronBtoFAltGen.root" , "MediumElectronBtoFAltGen" , "MediumElectronBCDEF", "2018-02-27/gsfElectronToMedium/", 2200,800, 0,0, 10,100);  
  scale_factors("2018-02-27/gsfElectronToMedium/tnp_BtoF_sigTemplate_bkgTemplate_LOMC/eff.root"   , "2018-02-27/gsfElectronToMedium/tnp_BtoF_mc_NLOMC/eff.root"         , "sfMediumElectronBtoFAltGen.root" , "MediumElectronBtoFAltGen" , "MediumElectronBCDEF", "2018-02-27/gsfElectronToMedium/", 2200,800, 0,0, 10,100);  
  scale_factors("2018-02-27/gsfElectronToMedium/tnp_BtoF_sigTemplate_bkgTemplate_altTag/eff.root" , "2018-02-27/gsfElectronToMedium/tnp_BtoF_mc_NLOMC_altTag/eff.root" , "sfMediumElectronBtoFAltTag.root" , "MediumElectronBtoFAltTag" , "MediumElectronBCDEF", "2018-02-27/gsfElectronToMedium/", 2200,800, 0,0, 10,100);  
  scale_factors("2018-02-27/gsfElectronToMedium/tnp_BtoF_sigTemplate_bkgTemplate_resFunc/eff.root", "2018-02-27/gsfElectronToMedium/tnp_BtoF_mc_NLOMC/eff.root"        , "sfMediumElectronBtoFResFunc.root", "MediumElectronBtoFResFunc", "MediumElectronBCDEF", "2018-02-27/gsfElectronToMedium/", 2200,800, 0,0, 10,100);  
  scale_factors("2018-02-27/gsfElectronToMedium/tnp_BtoF_sigTemplate_bkgTemplate_powhegPythia/eff.root", "2018-02-27/gsfElectronToMedium/tnp_BtoF_mc_NLOMC_powhegPythia/eff.root", "sfMediumElectronBtoFpowhegPythia.root", "MediumElectronBtoFpowhegPythia", "MediumElectronBCDEF", "2018-02-27/gsfElectronToMedium/", 2200,800, 0,0, 10,100);  
  scale_factors("2018-02-27/gsfElectronToMedium/tnp_BtoF_sigTemplate_bkgTemplate_powhegPhotos/eff.root", "2018-02-27/gsfElectronToMedium/tnp_BtoF_mc_NLOMC_powhegPhotos/eff.root", "sfMediumElectronBtoFpowhegPhotos.root", "MediumElectronBtoFpowhegPhotos", "MediumElectronBCDEF", "2018-02-27/gsfElectronToMedium/", 2200,800, 0,0, 10,100);  
  //scale_factors("2018-02-27/gsfElectronToMedium/tnp_BtoF_sigTemplate_bkgTemplate_fsrDown/eff.root", "2018-02-27/gsfElectronToMedium/tnp_BtoF_mc_NLOMC_fsrDown/eff.root", "sfMediumElectronBtoFfsrDown.root", "MediumElectronBtoFfsrDown", "MediumElectronBCDEF", "2018-02-27/gsfElectronToMedium/", 2200,800, 0,0, 10,100);  
  //scale_factors("2018-02-27/gsfElectronToMedium/tnp_BtoF_sigTemplate_bkgTemplate_fsrUp/eff.root"  , "2018-02-27/gsfElectronToMedium/tnp_BtoF_mc_NLOMC_fsrUp/eff.root"  , "sfMediumElectronBtoFfsrUp.root"  , "MediumElectronBtoFfsrUp"  , "MediumElectronBCDEF", "2018-02-27/gsfElectronToMedium/", 2200,800, 0,0, 10,100);  
  scale_factors("2018-02-27/gsfElectronToMedium/tnp_BtoF_sigTemplate_bkgTemplate/eff.root"        , "2018-02-27/gsfElectronToMedium/tnp_BtoF_mc_NLOMC/eff.root"        , "sfMediumElectronBtoFNominal.root", "MediumElectronBtoFNominal", "MediumElectronBCDEF", "2018-02-27/gsfElectronToMedium/", 2200,800, 0,0, 10,100);  
  scale_factors("2018-02-27/gsfElectronToMedium/tnp_GtoH_sigTemplate_bkgDasPlusExp/eff.root"      , "2018-02-27/gsfElectronToMedium/tnp_GtoH_mc_NLOMC/eff.root"        , "sfMediumElectronGtoHAltBkg.root" , "MediumElectronGtoHAltBkg" , "MediumElectronGH", "2018-02-27/gsfElectronToMedium/", 2200,800, 0,0, 10,100);  
  //scale_factors("2018-02-27/gsfElectronToMedium/tnp_GtoH_sigTemplate_bkgTemplate_LOMC/eff.root"   , "2018-02-27/gsfElectronToMedium/tnp_GtoH_mc_LOMC/eff.root"         , "sfMediumElectronGtoHAltGen.root" , "MediumElectronGtoHAltGen" , "MediumElectronGH", "2018-02-27/gsfElectronToMedium/", 2200,800, 0,0, 10,100);  
  scale_factors("2018-02-27/gsfElectronToMedium/tnp_GtoH_sigTemplate_bkgTemplate_LOMC/eff.root"   , "2018-02-27/gsfElectronToMedium/tnp_GtoH_mc_NLOMC/eff.root"         , "sfMediumElectronGtoHAltGen.root" , "MediumElectronGtoHAltGen" , "MediumElectronGH", "2018-02-27/gsfElectronToMedium/", 2200,800, 0,0, 10,100);  
  scale_factors("2018-02-27/gsfElectronToMedium/tnp_GtoH_sigTemplate_bkgTemplate_altTag/eff.root" , "2018-02-27/gsfElectronToMedium/tnp_GtoH_mc_NLOMC_altTag/eff.root" , "sfMediumElectronGtoHAltTag.root" , "MediumElectronGtoHAltTag" , "MediumElectronGH", "2018-02-27/gsfElectronToMedium/", 2200,800, 0,0, 10,100);  
  scale_factors("2018-02-27/gsfElectronToMedium/tnp_GtoH_sigTemplate_bkgTemplate_resFunc/eff.root", "2018-02-27/gsfElectronToMedium/tnp_GtoH_mc_NLOMC/eff.root"        , "sfMediumElectronGtoHResFunc.root", "MediumElectronGtoHResFunc", "MediumElectronGH", "2018-02-27/gsfElectronToMedium/", 2200,800, 0,0, 10,100);  
  scale_factors("2018-02-27/gsfElectronToMedium/tnp_GtoH_sigTemplate_bkgTemplate_powhegPythia/eff.root", "2018-02-27/gsfElectronToMedium/tnp_GtoH_mc_NLOMC_powhegPythia/eff.root", "sfMediumElectronGtoHpowhegPythia.root", "MediumElectronGtoHpowhegPythia", "MediumElectronGH", "2018-02-27/gsfElectronToMedium/", 2200,800, 0,0, 10,100);  
  scale_factors("2018-02-27/gsfElectronToMedium/tnp_GtoH_sigTemplate_bkgTemplate_powhegPhotos/eff.root", "2018-02-27/gsfElectronToMedium/tnp_GtoH_mc_NLOMC_powhegPhotos/eff.root", "sfMediumElectronGtoHpowhegPhotos.root", "MediumElectronGtoHpowhegPhotos", "MediumElectronGH", "2018-02-27/gsfElectronToMedium/", 2200,800, 0,0, 10,100);  
  //scale_factors("2018-02-27/gsfElectronToMedium/tnp_GtoH_sigTemplate_bkgTemplate_fsrDown/eff.root", "2018-02-27/gsfElectronToMedium/tnp_GtoH_mc_NLOMC_fsrDown/eff.root", "sfMediumElectronGtoHfsrDown.root", "MediumElectronGtoHfsrDown", "MediumElectronGH", "2018-02-27/gsfElectronToMedium/", 2200,800, 0,0, 10,100);  
  //scale_factors("2018-02-27/gsfElectronToMedium/tnp_GtoH_sigTemplate_bkgTemplate_fsrUp/eff.root"  , "2018-02-27/gsfElectronToMedium/tnp_GtoH_mc_NLOMC_fsrUp/eff.root"  , "sfMediumElectronGtoHfsrUp.root"  , "MediumElectronGtoHfsrUp"  , "MediumElectronGH", "2018-02-27/gsfElectronToMedium/", 2200,800, 0,0, 10,100);  
  scale_factors("2018-02-27/gsfElectronToMedium/tnp_GtoH_sigTemplate_bkgTemplate/eff.root"        , "2018-02-27/gsfElectronToMedium/tnp_GtoH_mc_NLOMC/eff.root"        , "sfMediumElectronGtoHNominal.root", "MediumElectronGtoHNominal", "MediumElectronGH", "2018-02-27/gsfElectronToMedium/", 2200,800, 0,0, 10,100);  
}

