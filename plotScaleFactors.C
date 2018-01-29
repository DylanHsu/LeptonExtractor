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
void scale_factors(string plots_dir, string root_dir, string basename_config, int width=2400, int height=1200, double xmin=0, double xmax=0, double ymin=0, double ymax=0){ 
  // This function calculates the efficiencies based on output from MIT TNP saved in subdirectories of plots_dir
  // The base names of the subdirectories are recorded in the config file
  // The efficiencies are plotted in plots_dir
  // The efficiencies and scale factors are recorded in root_dir in a 
  // rootfile whose filename is taken from basename_config


  // Pad directories with a slash at the end if it's not there
  if( plots_dir[plots_dir.size()-1]  != '/' ) plots_dir = plots_dir + "/";
  if( root_dir[root_dir.size()-1]  != '/' )  root_dir  = root_dir + "/";
  system( Form("mkdir -p %s", root_dir.c_str()) );

  // Open the output rootfile
  string output_rootfile_name;
  size_t ext_location = basename_config.find_last_of(".");
  if(ext_location != string::npos) output_rootfile_name = "scalefactors_"+basename_config.substr(0, ext_location)+".root";
  else output_rootfile_name = "scalefactors_"+basename_config+".root";
  TFile *output_rootfile = TFile::Open((root_dir+output_rootfile_name).c_str(), "RECREATE");
  assert(output_rootfile->IsOpen() && !output_rootfile->IsZombie());

  // Now read config file. Each line has
  // <data basename> <mc basename> <selection> <flavor>
  ifstream config_stream;
  config_stream.open(basename_config.c_str());
  assert(config_stream.is_open());
  string line, data_basename, mc_basename, output_basename, selection, flavor;
  while(getline(config_stream, line)) {
    stringstream ss(line);
    ss >> data_basename >> mc_basename >> selection >> flavor;

    if(flavor == "Electron")  output_basename = selection+"_Electron";
    else if(flavor == "Muon")   output_basename = selection+"_Muon";
    else                    output_basename = selection+"_"+flavor;
    
    TFile *f_data = TFile::Open( ( plots_dir + data_basename   + "/eff.root").c_str(), "READ");
    TH2D *h_eff_data      = (TH2D*) f_data->Get("hEffEtaPt");  h_eff_data     ->SetDirectory(0); 
    TH2D *h_error_lo_data = (TH2D*) f_data->Get("hErrlEtaPt"); h_error_lo_data->SetDirectory(0); 
    TH2D *h_error_hi_data = (TH2D*) f_data->Get("hErrhEtaPt"); h_error_hi_data->SetDirectory(0); 
    h_eff_data->SetName(("eff_data_"+output_basename).c_str());
    h_eff_data->SetTitle(("Efficiency for "+selection+" "+flavor+" selection (Data)").c_str());
    
    TFile *f_mc   = TFile::Open( ( plots_dir + mc_basename + "/eff.root").c_str(), "READ");
    TH2D *h_eff_mc      = (TH2D*) f_mc->Get("hEffEtaPt");  h_eff_mc     ->SetDirectory(0); 
    TH2D *h_error_lo_mc = (TH2D*) f_mc->Get("hErrlEtaPt"); h_error_lo_mc->SetDirectory(0); 
    TH2D *h_error_hi_mc = (TH2D*) f_mc->Get("hErrhEtaPt"); h_error_hi_mc->SetDirectory(0); 
    h_eff_mc->SetName(("eff_mc_"+output_basename).c_str());
    h_eff_mc->SetTitle(("Efficiency for "+selection+" "+flavor+" selection (MC)").c_str());
    
    // Divide Data/MC to get the scale factors
    TH2D *h_sf = (TH2D*) h_eff_data->Clone(); h_sf->SetDirectory(0); h_sf->Clear(); h_sf->Reset(); 
    h_sf->SetName(("scalefactors_"+output_basename).c_str()); 
    h_sf->SetTitle(("Scale factors for "+selection+" "+flavor+" selection").c_str());

    // Propagate asymmetric statistical errors
    TH2D *h_sf_error_lo = (TH2D*) h_sf->Clone(); h_sf_error_lo->SetDirectory(0); h_sf_error_lo->Clear(); h_sf_error_lo->Reset();
    TH2D *h_sf_error_hi = (TH2D*) h_sf->Clone(); h_sf_error_hi->SetDirectory(0); h_sf_error_hi->Clear(); h_sf_error_hi->Reset();
    h_sf_error_lo->SetName(("scalefactors_"+output_basename+"_stat_error_lo").c_str());
    h_sf_error_hi->SetName(("scalefactors_"+output_basename+"_stat_error_hi").c_str());
    h_sf_error_lo->SetTitle(("-1#sigma quantile stat. unc. for "+selection+" "+flavor+" scale factors").c_str());
    h_sf_error_hi->SetTitle(("+1#sigma quantile stat. unc. for "+selection+" "+flavor+" scale factors").c_str());

    for(int i = 1; i <= h_eff_mc->GetNbinsX(); i++) { for(int j = 1; j <= h_eff_mc->GetNbinsY(); j++) {
      unsigned int nbin = h_eff_mc->GetBin(i,j);
      double sf=h_eff_data->GetBinContent(nbin) / h_eff_mc->GetBinContent(nbin);
      double sf_error_hi, sf_error_lo;
      bool badSf=false; if(sf!=sf || sf<=0 || sf>2) { sf=1; badSf=true;}
      h_sf->SetBinContent(nbin, sf);
      if(!badSf) {
        sf_error_hi = sf * sqrt( pow(h_error_hi_data->GetBinContent(nbin) / h_eff_data->GetBinContent(nbin), 2) + pow(h_error_lo_mc->GetBinContent(nbin) / h_eff_mc->GetBinContent(nbin),2));
        sf_error_lo = sf * sqrt( pow(h_error_lo_data->GetBinContent(nbin) / h_eff_data->GetBinContent(nbin), 2) + pow(h_error_hi_mc->GetBinContent(nbin) / h_eff_mc->GetBinContent(nbin),2));
      } else {
        sf_error_hi = 1;
        sf_error_lo = 1;
      }
      h_sf_error_hi->SetBinContent(nbin, sf_error_hi);
      h_sf_error_lo->SetBinContent(nbin, sf_error_lo);
      // Choose the max for the s.f. and eff. histogram errors
      // Analyzers who naively use this uncertainty value will just get worse sensitivity :^)
      h_sf->SetBinError(nbin, TMath::Max(
        h_sf_error_hi->GetBinContent(nbin),
        h_sf_error_lo->GetBinContent(nbin)
      ));
      h_eff_data->SetBinError(nbin, TMath::Max(
        h_error_lo_data->GetBinContent(nbin),
        h_error_hi_data->GetBinContent(nbin)
      ));
      h_eff_mc->SetBinError(nbin, TMath::Max(
        h_error_lo_mc->GetBinContent(nbin),
        h_error_hi_mc->GetBinContent(nbin)
      ));
    }} 

    // Start drawing stuff 
    gStyle->SetOptStat(0);
    gStyle->SetPaintTextFormat("4.3f");
    TPaletteAxis *palette_axis; TCanvas *canvas[5]; TH2D* theHist[5]; int thePalette[5];
    thePalette[0]=kBlackBody; canvas[0]=new TCanvas("c_eff_data"   ,"c_eff_data"   ,width,height); theHist[0]=h_eff_data   ; h_eff_data   ->SetMinimum(  0); h_eff_data   ->SetMaximum(  1);
    thePalette[1]=kBlackBody; canvas[1]=new TCanvas("c_eff_mc"     ,"c_eff_mc"     ,width,height); theHist[1]=h_eff_mc     ; h_eff_mc     ->SetMinimum(  0); h_eff_mc     ->SetMaximum(  1);
    thePalette[2]=kBlackBody; canvas[2]=new TCanvas("c_sf"         ,"c_sf"         ,width,height); theHist[2]=h_sf         ; h_sf         ->SetMinimum(0.8); h_sf         ->SetMaximum(1.2);
    thePalette[3]=kCool;      canvas[3]=new TCanvas("c_sf_error_hi","c_sf_error_hi",width,height); theHist[3]=h_sf_error_hi; h_sf_error_hi->SetMinimum(  0); h_sf_error_hi->SetMaximum(0.1);
    thePalette[4]=kCool;      canvas[4]=new TCanvas("c_sf_error_lo","c_sf_error_lo",width,height); theHist[4]=h_sf_error_lo; h_sf_error_lo->SetMinimum(  0); h_sf_error_lo->SetMaximum(0.1);

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
      if(flavor=="Muon") gStyle->SetPaintTextFormat("4.2f"); else gStyle->SetPaintTextFormat("4.2f");
      if(flavor=="Muon") theHist[i]->Draw("TEXT45 COLZ"); else theHist[i]->Draw("TEXT45 COLZ");
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
      canvas[i]->Print(Form("%s%s_%.1f-%.1f_%.1f-%.1f.png",plots_dir.c_str(),theHist[i]->GetName(),xmin,xmax,ymin,ymax));
      canvas[i]->Print(Form("%s%s_%.1f-%.1f_%.1f-%.1f.pdf",plots_dir.c_str(),theHist[i]->GetName(),xmin,xmax,ymin,ymax));
    }

    for(unsigned i=0;i<5;i++) {
      theHist[i]->GetXaxis()->SetRangeUser(theHist[i]->GetXaxis()->GetBinLowEdge(1), theHist[i]->GetXaxis()->GetBinLowEdge(theHist[i]->GetXaxis()->GetNbins()+1));
      theHist[i]->GetYaxis()->SetRangeUser(theHist[i]->GetYaxis()->GetBinLowEdge(1), theHist[i]->GetYaxis()->GetBinLowEdge(theHist[i]->GetYaxis()->GetNbins()+1));
    }
    // Write to file
    output_rootfile->cd();
    h_eff_data->Write();
    h_eff_mc->Write();
    h_sf->Write();
    h_sf_error_lo->Write();
    h_sf_error_hi->Write();
    
    
    f_data->Close();
    f_mc->Close();
  }
  output_rootfile->Close();
  printf("Saved efficiencies, scale factors, and statistical errors in %s\n", (root_dir+output_rootfile_name).c_str());
}

void doAll() {
  system("rm 2017-10-19/*/*/*_Medium_Muon*.p*");
  system("rm 2017-10-19/*/*/*_Medium_Electron*.p*");
  scale_factors("2017-10-19/electrons/bkg_lpi_emu","2017-10-19/electrons/bkg_lpi_emu","ele_medium_2016.cfg",2400,1200,0,0,10,40);
  scale_factors("2017-10-19/electrons/bkg_lpi_emu","2017-10-19/electrons/bkg_lpi_emu","ele_medium_2016.cfg",2400,1200,0,0,40,100);
  scale_factors("2017-10-19/electrons/sig_fsrUp","2017-10-19/electrons/sig_fsrUp","ele_medium_2016.cfg",2400,1200,0,0,10,40);
  scale_factors("2017-10-19/electrons/sig_fsrUp","2017-10-19/electrons/sig_fsrUp","ele_medium_2016.cfg",2400,1200,0,0,40,100);
  scale_factors("2017-10-19/electrons/sig_fsrDown","2017-10-19/electrons/sig_fsrDown","ele_medium_2016.cfg",2400,1200,0,0,10,40);
  scale_factors("2017-10-19/electrons/sig_fsrDown","2017-10-19/electrons/sig_fsrDown","ele_medium_2016.cfg",2400,1200,0,0,40,100);
  scale_factors("2017-10-19/electrons/sig_genmc","2017-10-19/electrons/sig_genmc","ele_medium_2016.cfg",2400,1200,0,0,10,40);
  scale_factors("2017-10-19/electrons/sig_genmc","2017-10-19/electrons/sig_genmc","ele_medium_2016.cfg",2400,1200,0,0,40,100);
  scale_factors("2017-10-19/electrons/alt_tag","2017-10-19/electrons/alt_tag","ele_medium_2016.cfg",2400,1200,0,0,10,40);
  scale_factors("2017-10-19/electrons/alt_tag","2017-10-19/electrons/alt_tag","ele_medium_2016.cfg",2400,1200,0,0,40,100);
  scale_factors("2017-10-19/electrons/nominal","2017-10-19/electrons/nominal","ele_medium_2016.cfg",2400,1200,0,0,10,40);
  scale_factors("2017-10-19/electrons/nominal","2017-10-19/electrons/nominal","ele_medium_2016.cfg",2400,1200,0,0,40,100);
  scale_factors("2017-10-19/electrons/mc_lo","2017-10-19/electrons/mc_lo","ele_medium_2016_lo.cfg",2400,1200,0,0,10,40);
  scale_factors("2017-10-19/electrons/mc_lo","2017-10-19/electrons/mc_lo","ele_medium_2016_lo.cfg",2400,1200,0,0,40,100);
  
  scale_factors("2017-10-19/muons/bkg_lpi","2017-10-19/muons/bkg_lpi","mu_medium_2016.cfg",2400,1200,0,0,10,40);
  scale_factors("2017-10-19/muons/bkg_lpi","2017-10-19/muons/bkg_lpi","mu_medium_2016.cfg",2400,1200,0,0,40,50);
  scale_factors("2017-10-19/muons/bkg_lpi","2017-10-19/muons/bkg_lpi","mu_medium_2016.cfg",2400,1200,0,0,50,100);
  scale_factors("2017-10-19/muons/sig_fsrUp","2017-10-19/muons/sig_fsrUp","mu_medium_2016.cfg",2400,1200,0,0,10,40);
  scale_factors("2017-10-19/muons/sig_fsrUp","2017-10-19/muons/sig_fsrUp","mu_medium_2016.cfg",2400,1200,0,0,40,50);
  scale_factors("2017-10-19/muons/sig_fsrUp","2017-10-19/muons/sig_fsrUp","mu_medium_2016.cfg",2400,1200,0,0,50,100);
  scale_factors("2017-10-19/muons/sig_fsrDown","2017-10-19/muons/sig_fsrDown","mu_medium_2016.cfg",2400,1200,0,0,10,40);
  scale_factors("2017-10-19/muons/sig_fsrDown","2017-10-19/muons/sig_fsrDown","mu_medium_2016.cfg",2400,1200,0,0,40,50);
  scale_factors("2017-10-19/muons/sig_fsrDown","2017-10-19/muons/sig_fsrDown","mu_medium_2016.cfg",2400,1200,0,0,50,100);
  scale_factors("2017-10-19/muons/sig_genmc","2017-10-19/muons/sig_genmc","mu_medium_2016.cfg",2400,1200,0,0,10,40);
  scale_factors("2017-10-19/muons/sig_genmc","2017-10-19/muons/sig_genmc","mu_medium_2016.cfg",2400,1200,0,0,40,50);
  scale_factors("2017-10-19/muons/sig_genmc","2017-10-19/muons/sig_genmc","mu_medium_2016.cfg",2400,1200,0,0,50,100);
  scale_factors("2017-10-19/muons/alt_tag","2017-10-19/muons/alt_tag","mu_medium_2016.cfg",2400,1200,0,0,10,40);
  scale_factors("2017-10-19/muons/alt_tag","2017-10-19/muons/alt_tag","mu_medium_2016.cfg",2400,1200,0,0,40,50);
  scale_factors("2017-10-19/muons/alt_tag","2017-10-19/muons/alt_tag","mu_medium_2016.cfg",2400,1200,0,0,50,100);
  scale_factors("2017-10-19/muons/nominal","2017-10-19/muons/nominal","mu_medium_2016.cfg",2400,1200,0,0,10,40);
  scale_factors("2017-10-19/muons/nominal","2017-10-19/muons/nominal","mu_medium_2016.cfg",2400,1200,0,0,40,50);
  scale_factors("2017-10-19/muons/nominal","2017-10-19/muons/nominal","mu_medium_2016.cfg",2400,1200,0,0,50,100);
  scale_factors("2017-10-19/muons/mc_lo","2017-10-19/muons/mc_lo","mu_medium_2016_lo.cfg",2400,1200,0,0,10,40);
  scale_factors("2017-10-19/muons/mc_lo","2017-10-19/muons/mc_lo","mu_medium_2016_lo.cfg",2400,1200,0,0,40,50);
  scale_factors("2017-10-19/muons/mc_lo","2017-10-19/muons/mc_lo","mu_medium_2016_lo.cfg",2400,1200,0,0,50,100);
  system("eog 2017-10-19/*/*/scalefactors*_Medium_*on_0*.png &");
}
void doAll2() {
  system("rm 2018-01-08_splitEraPileup/*/*/*_Medium_Muon*.p*");
  system("rm 2018-01-08_splitEraPileup/*/*/*_Medium_Electron*.p*");
  scale_factors("2018-01-08_splitEraPileup/electrons/bkg_lpi_emu_BtoF","2018-01-08_splitEraPileup/electrons/bkg_lpi_emu_BtoF","ele_medium_2016_BtoF.cfg",2400,1200,0,0,10,40);
  scale_factors("2018-01-08_splitEraPileup/electrons/bkg_lpi_emu_BtoF","2018-01-08_splitEraPileup/electrons/bkg_lpi_emu_BtoF","ele_medium_2016_BtoF.cfg",2400,1200,0,0,40,100);
  scale_factors("2018-01-08_splitEraPileup/electrons/bkg_lpi_emu_GtoH","2018-01-08_splitEraPileup/electrons/bkg_lpi_emu_GtoH","ele_medium_2016_GtoH.cfg",2400,1200,0,0,10,40);
  scale_factors("2018-01-08_splitEraPileup/electrons/bkg_lpi_emu_GtoH","2018-01-08_splitEraPileup/electrons/bkg_lpi_emu_GtoH","ele_medium_2016_GtoH.cfg",2400,1200,0,0,40,100);
  scale_factors("2018-01-08_splitEraPileup/muons/bkg_lpi_BtoF","2018-01-08_splitEraPileup/muons/bkg_lpi_BtoF","mu_medium_2016_BtoF.cfg",2400,1200,0,0,10,40);
  scale_factors("2018-01-08_splitEraPileup/muons/bkg_lpi_BtoF","2018-01-08_splitEraPileup/muons/bkg_lpi_BtoF","mu_medium_2016_BtoF.cfg",2400,1200,0,0,40,50);
  scale_factors("2018-01-08_splitEraPileup/muons/bkg_lpi_BtoF","2018-01-08_splitEraPileup/muons/bkg_lpi_BtoF","mu_medium_2016_BtoF.cfg",2400,1200,0,0,50,100);
  scale_factors("2018-01-08_splitEraPileup/muons/bkg_lpi_GtoH","2018-01-08_splitEraPileup/muons/bkg_lpi_GtoH","mu_medium_2016_GtoH.cfg",2400,1200,0,0,10,40);
  scale_factors("2018-01-08_splitEraPileup/muons/bkg_lpi_GtoH","2018-01-08_splitEraPileup/muons/bkg_lpi_GtoH","mu_medium_2016_GtoH.cfg",2400,1200,0,0,40,50);
  scale_factors("2018-01-08_splitEraPileup/muons/bkg_lpi_GtoH","2018-01-08_splitEraPileup/muons/bkg_lpi_GtoH","mu_medium_2016_GtoH.cfg",2400,1200,0,0,50,100);
  //system("eog 2017-10-19/*/*/scalefactors*_Medium_*on_0*.png &");
}
