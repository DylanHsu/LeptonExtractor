#include <vector>
#include <fstream>
#include <TString.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TTree.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <THStack.h>
#include <TLorentzVector.h>

void systEffectOnXs(string flavor="electrons") {
 TString folder="/data/t3home000/dhsu/root/lepTrack2/";
 double lumi=35900.;// inverse picobarns

 vector<TString> fileName_; vector<double> xs_; vector<unsigned> flavor_;
 // 0: data, 1: Zmm/Zee, 2: Z->TauTau, 3: TTbar, 4: W+jets
 if(flavor=="electrons") {
  flavor_.push_back(11); xs_.push_back(  6025.2/2. ); fileName_.push_back("DYJetsToLL_BaselineToMedium_ScaleSmearCorrections_PtInclusiveNLO_electronTnP.root");
  flavor_.push_back(11); xs_.push_back(   5451./2. ); fileName_.push_back("DYJetsToLL_BaselineToMedium_ScaleSmearCorrections_Zpt-0To50-NLO_electronTnP.root");
  flavor_.push_back(11); xs_.push_back(   355.3/2. ); fileName_.push_back("DYJetsToLL_BaselineToMedium_ScaleSmearCorrections_Pt-50To100-NLO_electronTnP.root");
  flavor_.push_back(11); xs_.push_back(   81.02/2. ); fileName_.push_back("DYJetsToLL_BaselineToMedium_ScaleSmearCorrections_Pt-100To250-NLO_electronTnP.root");
  flavor_.push_back(11); xs_.push_back(   2.991/2. ); fileName_.push_back("DYJetsToLL_BaselineToMedium_ScaleSmearCorrections_Pt-250To400-NLO_electronTnP.root");
  flavor_.push_back(11); xs_.push_back(  0.3882/2. ); fileName_.push_back("DYJetsToLL_BaselineToMedium_ScaleSmearCorrections_Pt-400To650-NLO_electronTnP.root");
  flavor_.push_back(11); xs_.push_back( 0.03737/2. ); fileName_.push_back("DYJetsToLL_BaselineToMedium_ScaleSmearCorrections_Pt-650ToInf-NLO_electronTnP.root");
 } else {
  flavor_.push_back(13); xs_.push_back(  6025.2/2. ); fileName_.push_back("DYJetsToLL_BaselineToMedium_RochesterCorrections_PtInclusiveNLO_muonTnP.root");
  flavor_.push_back(13); xs_.push_back(   5451./2. ); fileName_.push_back("DYJetsToLL_BaselineToMedium_RochesterCorrections_Zpt-0To50-NLO_muonTnP.root");
  flavor_.push_back(13); xs_.push_back(   355.3/2. ); fileName_.push_back("DYJetsToLL_BaselineToMedium_RochesterCorrections_Pt-50To100-NLO_muonTnP.root");
  flavor_.push_back(13); xs_.push_back(   81.02/2. ); fileName_.push_back("DYJetsToLL_BaselineToMedium_RochesterCorrections_Pt-100To250-NLO_muonTnP.root");
  flavor_.push_back(13); xs_.push_back(   2.991/2. ); fileName_.push_back("DYJetsToLL_BaselineToMedium_RochesterCorrections_Pt-250To400-NLO_muonTnP.root");
  flavor_.push_back(13); xs_.push_back(  0.3882/2. ); fileName_.push_back("DYJetsToLL_BaselineToMedium_RochesterCorrections_Pt-400To650-NLO_muonTnP.root");
  flavor_.push_back(13); xs_.push_back( 0.03737/2. ); fileName_.push_back("DYJetsToLL_BaselineToMedium_RochesterCorrections_Pt-650ToInf-NLO_muonTnP.root");
 }
 TFile *sfFile=0;
 if(flavor=="electrons") { sfFile=TFile::Open("2017-10-19/electrons/bkg_lpi_emu/scalefactors_ele_medium_2016.root","READ"); assert(sfFile); }
 else                    { sfFile=TFile::Open("2017-10-19/muons/bkg_lpi/scalefactors_mu_medium_2016.root","READ"); assert(sfFile); }
 vector<TString> systSource_; vector<double> weightedYieldUp_; vector<double> weightedYieldDown_;
 systSource_.push_back("stat_error_hi"     ); weightedYieldUp_.push_back(0); weightedYieldDown_.push_back(0); 
 systSource_.push_back("stat_error_lo"     ); weightedYieldUp_.push_back(0); weightedYieldDown_.push_back(0); 
 systSource_.push_back("signalFsrTNP"      ); weightedYieldUp_.push_back(0); weightedYieldDown_.push_back(0); 
 systSource_.push_back("signalResTNP"      ); weightedYieldUp_.push_back(0); weightedYieldDown_.push_back(0); 
 systSource_.push_back("bkgModelTNP"       ); weightedYieldUp_.push_back(0); weightedYieldDown_.push_back(0); 
 systSource_.push_back("tagBiasTNP"        ); weightedYieldUp_.push_back(0); weightedYieldDown_.push_back(0); 
 systSource_.push_back("generatorChoiceTNP"); weightedYieldUp_.push_back(0); weightedYieldDown_.push_back(0); 
 double nominalYield=0; double sumWeightsSquared=0;
 TH2D *h_sf=0;
 if(flavor=="electrons") h_sf=(TH2D*)sfFile->Get("scalefactors_Medium_Electron");
 else                    h_sf=(TH2D*)sfFile->Get("scalefactors_Medium_Muon");
 assert(h_sf); h_sf->SetDirectory(0);

 TH2D *h_syst_[16];
 for(unsigned iSource=0; iSource<systSource_.size(); iSource++) {
   if(flavor=="electrons") h_syst_[iSource]=(TH2D*)sfFile->Get(Form("scalefactors_Medium_Electron_%s",systSource_[iSource].Data())); 
   else                    h_syst_[iSource]=(TH2D*)sfFile->Get(Form("scalefactors_Medium_Muon_%s",systSource_[iSource].Data()));
   assert(h_syst_[iSource]); h_syst_[iSource]->SetDirectory(0);
 }
 sfFile->Close();

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
  sum_weights=(TH1D*)inputFile->Get("sum_weights"); assert(sum_weights);
  double norm= lumi*xs_[iFile]/sum_weights->GetBinContent(1);
  for(unsigned int ientry=0; ientry<(unsigned int)events->GetEntries(); ientry++) {
   events->GetEntry(ientry); if(ientry%1000000==0) printf("\tReading entry %d/%d\n",ientry,(unsigned int)events->GetEntries());
   if(pass==0) continue;
   if(probe->Pt()<25 || tag->Pt()<25) continue;
   if(flavor=="electrons" && TMath::Abs(tagPid)!=11) continue;
   if(flavor=="muons"     && TMath::Abs(tagPid)!=13) continue;
   if(probePid+tagPid!=0) continue;
   if(mass<76.1876||mass>106.1876) continue;
   double tagAbsEta=5, probeAbsEta=5;
   if(tag->Pt()>0) tagAbsEta=TMath::Abs(tag->Eta()); 
   if(probe->Pt()>0) probeAbsEta=TMath::Abs(probe->Eta());
   if(tagAbsEta>2.4 || probeAbsEta>2.4) continue;
   if(flavor=="electrons" && (tagAbsEta>1.4442&&tagAbsEta<1.566)) continue;
   if(flavor=="electrons" && (probeAbsEta>1.4442&&probeAbsEta<1.566)) continue;
   double weight=scale1fb*norm;
   double sfTag=h_sf->GetBinContent(h_sf->FindBin(tag->Eta(), tag->Pt()));
   double sfProbe=h_sf->GetBinContent(h_sf->FindBin(probe->Eta(), probe->Pt()));
   nominalYield += weight*sfTag*sfProbe;
   sumWeightsSquared += weight*weight*sfTag*sfTag*sfProbe*sfProbe;
   for(unsigned iSource=0; iSource<systSource_.size(); iSource++) {
     double sfUncTag=h_syst_[iSource]->GetBinContent(h_syst_[iSource]->FindBin(tag->Eta(), tag->Pt()));
     double sfUncProbe=h_syst_[iSource]->GetBinContent(h_syst_[iSource]->FindBin(probe->Eta(), probe->Pt()));
     weightedYieldUp_[iSource] += weight*(sfTag+sfUncTag)*(sfProbe+sfUncProbe);
     weightedYieldDown_[iSource] += weight*(sfTag-sfUncTag)*(sfProbe-sfUncProbe);
   } 
  }
  inputFile->Close(); bytesRead+=fileBytes_[iFile];
  unsigned long int time1 = static_cast<unsigned long int>(time(NULL)); time1=time1-time0;
  printf("\tDone, Closing file (%.1f%% complete, ETA %d s) \n",TMath::Max(.0001,100.*((double)bytesRead/(double)bytesToRead)), (int)((bytesToRead-bytesRead)*(double)time1/(double)bytesRead));
 }
 printf("Nominal yield: %f +/- %f\n", nominalYield, sqrt(sumWeightsSquared));
 for(unsigned iSource=0; iSource<systSource_.size(); iSource++) {
  printf("Yield from syst source %s : %f up, %f down (+%.2f%%,-%.2f%%)\n",
   systSource_[iSource].Data(), 
   weightedYieldUp_[iSource], 
   weightedYieldDown_[iSource],
   100.*(weightedYieldUp_[iSource]/nominalYield-1.), 
   100.*(1.-weightedYieldDown_[iSource]/nominalYield)
  );
 } 
}
