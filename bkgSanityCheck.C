#include <vector>
#include <fstream>
#include <TString.h>
#include <TFile.h>
#include <TH1D.h>
#include <TTree.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <THStack.h>
#include <TLorentzVector.h>

void bkgSanityCheck() {
 TString folder="/data/t3home000/dhsu/root/lepTrack2/";
 double lumi=35900.;// inverse picobarns
 unsigned nTypes=5;
 vector<int> processColor_; vector<TString> processName_;
 processColor_.push_back(    kBlack ); processName_.push_back("data");
 processColor_.push_back(       901 ); processName_.push_back("Zll");
 processColor_.push_back(  kAzure-9 ); processName_.push_back("Ztautau");
 processColor_.push_back(   kPink+7 ); processName_.push_back("TTTo2L2Nu");
 processColor_.push_back( kViolet+8 ); processName_.push_back("WJetsToLNu");

 vector<TString> fileName_; vector<double> xs_; vector<unsigned> processType_;
 processName_.push_back("data");
 processName_.push_back("Zll");
 processName_.push_back("Ztautau");
 processName_.push_back("TTTo2L2Nu");
 processName_.push_back("WJetsToLNu");
 // 0: data, 1: Zmm/Zee, 2: Z->TauTau, 3: TTbar, 4: W+jets
 processType_.push_back(0); xs_.push_back(      -1 ); fileName_.push_back("SingleElectron2016_BaselineToMedium_ScaleSmearCorrections_electronTnP.root");
 processType_.push_back(0); xs_.push_back(      -1 ); fileName_.push_back("SingleMuon2016_BaselineToMedium_RochesterCorrections_muonTnP.root");
 processType_.push_back(1); xs_.push_back(  6025.2 ); fileName_.push_back("DYJetsToLL_BaselineToMedium_ScaleSmearCorrections_PtInclusiveNLO_electronTnP.root");
 processType_.push_back(1); xs_.push_back(   5451. ); fileName_.push_back("DYJetsToLL_BaselineToMedium_ScaleSmearCorrections_Zpt-0To50-NLO_electronTnP.root");
 processType_.push_back(1); xs_.push_back(   355.3 ); fileName_.push_back("DYJetsToLL_BaselineToMedium_ScaleSmearCorrections_Pt-50To100-NLO_electronTnP.root");
 processType_.push_back(1); xs_.push_back(   81.02 ); fileName_.push_back("DYJetsToLL_BaselineToMedium_ScaleSmearCorrections_Pt-100To250-NLO_electronTnP.root");
 processType_.push_back(1); xs_.push_back(   2.991 ); fileName_.push_back("DYJetsToLL_BaselineToMedium_ScaleSmearCorrections_Pt-250To400-NLO_electronTnP.root");
 processType_.push_back(1); xs_.push_back(  0.3882 ); fileName_.push_back("DYJetsToLL_BaselineToMedium_ScaleSmearCorrections_Pt-400To650-NLO_electronTnP.root");
 processType_.push_back(1); xs_.push_back( 0.03737 ); fileName_.push_back("DYJetsToLL_BaselineToMedium_ScaleSmearCorrections_Pt-650ToInf-NLO_electronTnP.root");
 processType_.push_back(2); xs_.push_back(   1867. ); fileName_.push_back("DYJetsToTauTau_ForcedMuEleDecay_BaselineToMedium_electronTnP.root");
 processType_.push_back(3); xs_.push_back(    76.7 ); fileName_.push_back("TTTo2L2Nu_BaselineToMedium_electronTnP.root");
 processType_.push_back(4); xs_.push_back(  57280. ); fileName_.push_back("WJetsToLNu_Wpt-0To50_BaselineToMedium_electronTnP.root");
 processType_.push_back(4); xs_.push_back(   3258. ); fileName_.push_back("WJetsToLNu_Wpt-50To100_BaselineToMedium_electronTnP.root");
 processType_.push_back(4); xs_.push_back(   682.2 ); fileName_.push_back("WJetsToLNu_Pt-100To250_BaselineToMedium_electronTnP.root");
 processType_.push_back(4); xs_.push_back(   24.10 ); fileName_.push_back("WJetsToLNu_Pt-250To400_BaselineToMedium_electronTnP.root");
 processType_.push_back(4); xs_.push_back(   3.054 ); fileName_.push_back("WJetsToLNu_Pt-400To600_BaselineToMedium_electronTnP.root");
 processType_.push_back(4); xs_.push_back(  0.4590 ); fileName_.push_back("WJetsToLNu_Pt-600ToInf_BaselineToMedium_electronTnP.root");

 TH1D *h_mass[5], *h_massForward[5];
 for(unsigned n=0; n<nTypes; n++) {
  h_mass[n] = new TH1D(Form("h_mass_%s",processName_[n].Data()),Form("h_mass_%s",processName_[n].Data()),30,60,120); h_mass[n]->Sumw2();
  h_massForward[n] = new TH1D(Form("h_massForward_%s",processName_[n].Data()),Form("h_massForward_%s",processName_[n].Data()),30,60,120); h_massForward[n]->Sumw2();
  if(n==0) { h_mass[n]->SetLineColor(kBlack); h_mass[n]->SetMarkerStyle(20); }
  else     { h_mass[n]->SetLineColor(processColor_[n]); h_mass[n]->SetFillColor(processColor_[n]); h_mass[n]->SetFillStyle(1001);}
 }
 TH1D *h_epiData = new TH1D("h_epiData","h_epiData",30,60,120); h_epiData->Sumw2();
 TH1D *h_emuData = new TH1D("h_emuData","h_emuData",30,60,120); h_emuData->Sumw2();
 TH1D *h_epiDataForward = new TH1D("h_epiDataForward","h_epiDataForward",30,60,120); h_epiDataForward->Sumw2();
 TH1D *h_emuDataForward = new TH1D("h_emuDataForward","h_emuDataForward",30,60,120); h_emuDataForward->Sumw2();
 
 TH1D *mcSum=(TH1D*)h_mass[0]->Clone("mcSum");
 TFile *inputFile=0; TTree *events=0; TH1D *sum_weights;
 unsigned long bytesToRead=0; vector<unsigned long> fileBytes_; unsigned long bytesRead=0;
 for(unsigned iFile=0; iFile<fileName_.size(); iFile++) {
   std::ifstream in(Form("%s%s",folder.Data(),fileName_[iFile].Data()), std::ios::binary | std::ios::ate);
   fileBytes_.push_back((unsigned long)in.tellg()); bytesToRead+=((unsigned long)in.tellg());
 }
 unsigned int runNum, lumiSec, evtNum;   // event ID
 unsigned int npv;                       // number of primary vertices
 unsigned int pass;                      // whether probe passes requirements
 float        npu;                       // mean number of expected pileup
 float        scale1fb;                  // event weight per 1/fb
 float        mass;                      // tag-probe mass
 int          qtag, qprobe;              // tag, probe charge
 int          tagPid, probePid;              // tag, probe PID
 float        met;
 int          njets;
 TLorentzVector *tag=0, *probe=0;        // tag, probe 4-vector
 TLorentzVector *genp4_tag=0, *genp4_probe=0;        // tag, probe 4-vector
 std::time_t t = std::time(0);
 unsigned long int time0 = static_cast<unsigned long int>(time(NULL));
 for(unsigned iFile=0; iFile<fileName_.size(); iFile++) {
  printf("Opening \"%s\"...\n",fileName_[iFile].Data());
  inputFile=TFile::Open(Form("%s%s",folder.Data(),fileName_[iFile].Data()),"read"); assert(inputFile && inputFile->IsOpen());
  events=(TTree*)inputFile->Get("events"); assert(events);
  
  //events->SetBranchAddress("runNum",   &runNum);
  //events->SetBranchAddress("lumiSec",  &lumiSec);
  //events->SetBranchAddress("evtNum",   &evtNum);
  events->SetBranchAddress("npv",      &npv);
  events->SetBranchAddress("pass",     &pass);
  //events->SetBranchAddress("npu",      &npu);
  events->SetBranchAddress("scale1fb", &scale1fb);
  events->SetBranchAddress("mass",     &mass);
  events->SetBranchAddress("qtag",     &qtag);
  events->SetBranchAddress("qprobe",   &qprobe);
  events->SetBranchAddress("tagPid",     &tagPid);
  events->SetBranchAddress("probePid",   &probePid);
  events->SetBranchAddress("tag",      &tag);
  events->SetBranchAddress("probe",    &probe);
  bool isData=(processType_[iFile]==0);
  if(!isData) sum_weights=(TH1D*)inputFile->Get("sum_weights"); assert(sum_weights);
  double norm= !isData ? lumi*xs_[iFile]/sum_weights->GetBinContent(1) : 1;
  if(processType_[iFile]==1) norm/=2.;
  for(unsigned int ientry=0; ientry<(unsigned int)events->GetEntries(); ientry++) {
   events->GetEntry(ientry); if(ientry%1000000==0) printf("\tReading entry %d/%d\n",ientry,(unsigned int)events->GetEntries());
   if(pass==1) continue;
   if(probe->Pt()<25 || probe->Pt()>=50) continue;
   bool isForward=(probe->Eta()>1.4442 || probe->Eta()<-1.4442);
   if((tagPid==11||tagPid==-11) && probePid+tagPid==0) {
    if(isForward) h_massForward[processType_[iFile]]->Fill(mass, isData?norm:scale1fb*norm);
    else          h_mass       [processType_[iFile]]->Fill(mass, isData?norm:scale1fb*norm);
   } else if(isData && (tagPid==11||tagPid==-11) && (probePid==211 || probePid==-211) && qtag+qprobe!=0) {
    if(isForward) h_epiDataForward->Fill(mass, isData?norm:scale1fb*norm);
    else          h_epiData       ->Fill(mass, isData?norm:scale1fb*norm);
   } else if(isData && (tagPid==11||tagPid==-11) && (probePid==13||probePid==-13)) { // bug in ntuples
    if(isForward) h_emuDataForward->Fill(mass, isData?norm:scale1fb*norm);
    else          h_emuData       ->Fill(mass, isData?norm:scale1fb*norm);
   }
  }
  //if(isData) events->Draw(
  // Form("mass>>+mass_%s",processName_[processType_[iFile]].Data()),
  // "TMath::Abs(tagPid)==11 && TMath::Abs(probePid)==211 && qtag+qprobe!=0 && pass==0 && probe.Pt()>25 && probe.Pt()<50 && TMath::Abs(probe.Eta())<1.4442",
  // "goff e");
  //else events->Draw(
  // Form("mass>>+mass_%s",processName_[processType_[iFile]].Data()),
  // Form("%f*scale1fb*(TMath::Abs(tagPid)==11 && (probePid+tagPid==0) && pass==0 && probe.Pt()>25 && probe.Pt()<50 && TMath::Abs(probe.Eta())<1.4442)",norm),
  // "goff e");
  inputFile->Close(); bytesRead+=fileBytes_[iFile];
  unsigned long int time1 = static_cast<unsigned long int>(time(NULL)); time1=time1-time0;
  printf("\tDone, Closing file (%.1f%% complete, ETA %d s) \n",TMath::Max(.0001,100.*((double)bytesRead/(double)bytesToRead)), (int)((bytesToRead-bytesRead)*(double)time1/(double)bytesRead));
 }
 // scaling
 TH1D *mcError=0;
 {
  double Ndata=h_mass[0]->Integral(1,30);
  mcSum->Add(h_mass[4]); mcSum->Add(h_mass[3]); mcSum->Add(h_mass[2]);
  double totalBkg=mcSum->Integral(1,30);
  mcSum->Add(h_mass[1]);
  double totalMC=mcSum->Integral(1,30);
  double mcSF = Ndata/totalMC; totalBkg*=mcSF;
  h_mass[1]->Scale(mcSF); h_mass[2]->Scale(mcSF); h_mass[3]->Scale(mcSF); h_mass[4]->Scale(mcSF);
  double epiSF = .6*totalBkg/h_epiData->Integral(1,30);
  double emuSF = .4*totalBkg/h_emuData->Integral(1,30);
  h_epiData->Scale(epiSF); h_emuData->Scale(emuSF);
  mcError=(TH1D*)mcSum->Clone("mcError"); mcError->Scale(mcSF);
 }
 gStyle->SetOptStat(0);
 THStack *hs_mass=new THStack("hs_mass","Z(ee)"); 
 hs_mass->Add(h_mass[4]); hs_mass->Add(h_mass[3]); hs_mass->Add(h_mass[2]); //hs_mass->Add(h_mass[1]);
 //h_mass[0]->Scale(mcSum->Integral()/h_mass[0]->Integral());
 TCanvas *canvas=new TCanvas("c_mass","c_mass");
 hs_mass->Draw("hist"); 
 hs_mass->GetXaxis()->SetTitle("Invariant mass [GeV]");
 hs_mass->GetYaxis()->SetTitle("Events / 2 GeV");
 hs_mass->SetMinimum(0);
 //hs_mass->SetMaximum(1.5*h_mass[0]->GetMaximum());
 mcError->SetMarkerColor(901); mcError->SetMarkerStyle(1); mcError->SetFillStyle(3254); mcError->SetFillColor(kBlack);
 //mcError->Draw("e2 same");
 h_epiData->SetLineColor(kOrange+2); h_epiData->SetLineWidth(2); h_epiData->SetLineStyle(2);
 h_emuData->SetLineColor(kRed+3   ); h_emuData->SetLineWidth(2); h_emuData->SetLineStyle(9);
 TH1D *h_dataTemplate = (TH1D*) h_epiData->Clone("h_dataTemplate"); h_dataTemplate->Add(h_emuData);
 h_dataTemplate->SetLineColor(12); h_dataTemplate->SetLineStyle(3);
 hs_mass->SetMaximum(1.5*h_dataTemplate->GetMaximum());
 h_epiData->Draw("hist same"); h_emuData->Draw("hist same"); h_dataTemplate->Draw("hist same");
 //h_mass[0]->Draw("p e0 same");
}
