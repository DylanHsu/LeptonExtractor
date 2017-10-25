#include <TH1D.h>
#include <TTree.h>
#include <TFile.h>
#include <TString.h>
#include <TLorentzVector.h>
#include <cassert>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <unistd.h>

//const Float_t elePtBins[] = {10,15,20,25, 30, 35, 40, 50, 70, 100, 1000};
//const Float_t eleEtaBins[] = {-2.5,-2.0,-1.566,-1.4442,-1,-0.5,0,0.5,1,1.4442,1.566,2.0,2.5};
//const Int_t nElePtBins=10;
//const Int_t nEleEtaBins=12;
//const Float_t muPtBins[] = {10,15,20,25,30,35,40,50,70,100,1000};
//const Float_t muEtaBins[] = {-2.4,-2.1,-1.8,-1.6,-1.3,-1.1,-0.8,-0.6,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.6,0.8,1.1,1.3,1.6,1.8,2.1,2.4};
//const Int_t nMuPtBins=10;
//const Int_t nMuEtaBins=24;

// Selection: Zll | em | e+e+ | m+m+ | pm | pe
enum massSelType { kZll, kEM, kZeeSameSign, kZmmSameSign, kMuPi, kEPi, kGenZll, knSelTypes};
void binSplitter(
  TString inputFileName, 
  TString outputFileName, 
  std::string binFile, 
  double scalefactor=1, 
  bool isMC=true,
  massSelType selection = kZll,
  bool alt_tag=false
) {
  if(inputFileName==outputFileName) { printf("check output file name!\n"); return; }
  // parse bin file
  std::vector<double> fPtBinEdgesv, fEtaBinEdgesv, fPhiBinEdgesv, fNPVBinEdgesv, fJetsBinEdgesv, fMETBinEdgesv;
  // flags for |eta| and |phi| binning
  bool fDoAbsEta, fDoAbsPhi;
  // flags for binnings to compute efficiencies for
  bool fDoPt, fDoEta, fDoPhi, fDoEtaPt, fDoEtaPhi, fDoNPV, fDoJets, fDoMET;
  
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


  TFile *inputFile=TFile::Open(inputFileName,"READ");
  TTree *inputTree=(TTree*)inputFile->Get("events");
  // set up input tree
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
  
  //inputTree->SetBranchAddress("runNum",   &runNum);
  //inputTree->SetBranchAddress("lumiSec",  &lumiSec);
  //inputTree->SetBranchAddress("evtNum",   &evtNum);
  inputTree->SetBranchAddress("npv",      &npv);
  inputTree->SetBranchAddress("pass",     &pass);
  //inputTree->SetBranchAddress("npu",      &npu);
  inputTree->SetBranchAddress("scale1fb", &scale1fb);
  inputTree->SetBranchAddress("mass",     &mass);
  inputTree->SetBranchAddress("qtag",     &qtag);
  inputTree->SetBranchAddress("qprobe",   &qprobe);
  inputTree->SetBranchAddress("tagPid",     &tagPid);
  inputTree->SetBranchAddress("probePid",   &probePid);
  inputTree->SetBranchAddress("tag",      &tag);
  inputTree->SetBranchAddress("probe",    &probe);
  if(selection==kGenZll) inputTree->SetBranchAddress("genTag",      &genp4_tag);
  if(selection==kGenZll) inputTree->SetBranchAddress("genProbe",    &genp4_probe);
  //inputTree->SetBranchAddress("met",      &met);
  //inputTree->SetBranchAddress("njets",    &njets);
  
  TFile *outputFile = TFile::Open(outputFileName,"RECREATE");
  TH1D *histosPass[2048], *histosFail[2048];
  assert(outputFile && outputFile->IsOpen());
  int iHisto=0;
  for(unsigned iPt=0;iPt<fPtBinEdgesv.size()-1; iPt++) { for(unsigned iEta=0; iEta<fEtaBinEdgesv.size()-1; iEta++) {
    char histName[256];
    sprintf(histName,"pass_ptBin%d_etaBin%d", iPt, iEta);
    histosPass[iHisto] = new TH1D(histName,histName,50,40,140);
    sprintf(histName,"fail_ptBin%d_etaBin%d", iPt, iEta);
    histosFail[iHisto] = new TH1D(histName,histName,50,40,140);
    iHisto++;
  }}
  int iHistoMax=iHisto-1;
  TH1F *sum_weights = (TH1F*)inputFile->Get("sum_weights");
  assert(!isMC || sum_weights);
  
  for(unsigned int ientry=0; ientry<(unsigned int)inputTree->GetEntries(); ientry++) {
    inputTree->GetEntry(ientry);
    if(ientry%1000000==0) printf("reading entry %d/%d\n",ientry,(unsigned int)inputTree->GetEntries());
    if(mass<40.||mass>=140.) continue;
    bool passSelection=false;
    if     (selection==         kZll ) passSelection=(tagPid+probePid==0 && (tagPid==13||tagPid==-13||tagPid==11||tagPid==-11));
    else if(selection==      kGenZll ) passSelection=(tagPid+probePid==0 && (tagPid==13||tagPid==-13||tagPid==11||tagPid==-11));
    else if(selection==          kEM ) passSelection=((TMath::Abs(tagPid)==13 && TMath::Abs(probePid)==11)||(TMath::Abs(tagPid)==11 && TMath::Abs(probePid==13)))&&(tag->Pt()>=30.);
    else if(selection== kZeeSameSign ) passSelection=(tagPid==11&&probePid==11)||(tagPid==-11&&probePid==-11);
    else if(selection== kZmmSameSign ) passSelection=(tagPid==13&&probePid==13)||(tagPid==-13&&probePid==-13);
    else if(selection==        kMuPi ) passSelection=(tagPid==13&&probePid==-211)||(tagPid==-13&&probePid==211);
    else if(selection==         kEPi ) passSelection=(tagPid==11&&probePid==-211)||(tagPid==-11&&probePid==211);
    else                               passSelection=false;
    if(alt_tag) passSelection &= (
      ((tagPid==13||tagPid==-13) && tag->Pt()>30 && (tag->Eta()<2.1||tag->Eta()>-2.1)) ||
      ((tagPid==11||tagPid==-11) && tag->Pt()>35 && (tag->Eta()<2.1||tag->Eta()>-2.1))
    );
    if(!passSelection) continue;
    iHisto=0;
    bool ignoreEtaHighEnergy=(selection==kMuPi || selection==kEPi || selection==kEM);
    bool ignorePassFlag=(selection==kMuPi || selection==kEPi || selection==kEM);
    for(unsigned iPt=0;iPt<fPtBinEdgesv.size()-1; iPt++) { for(unsigned iEta=0; iEta<fEtaBinEdgesv.size()-1; iEta++) {
      double pt=probe->Pt(); 
      if(pt >= fPtBinEdgesv[iPt] && pt < fPtBinEdgesv[iPt+1]) { 
      double eta = fDoAbsEta? fabs(probe->Eta()) : probe->Eta();
      bool passEta=(eta >= fEtaBinEdgesv[iEta] && eta < fEtaBinEdgesv[iEta+1]) || (pt >= 50 && ignoreEtaHighEnergy);
        if(passEta) {
          double weight;
          if(isMC) weight = scale1fb / sum_weights->GetBinContent(1) * scalefactor;
          else     weight = 1;
          if(selection==kGenZll) { TLorentzVector genDilep = (*genp4_tag) + (*genp4_probe); mass=genDilep.M(); }
          if(pass  || ignorePassFlag) histosPass[iHisto]->Fill(mass,weight);
          if(!pass || ignorePassFlag) histosFail[iHisto]->Fill(mass,weight);
        }
      }
      iHisto++;
    }}
  }
  for(iHisto=0; iHisto <= iHistoMax; iHisto++) histosPass[iHisto]->Write();
  for(iHisto=0; iHisto <= iHistoMax; iHisto++) histosFail[iHisto]->Write();
  outputFile->Close();
}

void generateJobArgs(string outDir, std::string binFile, string flavor="electrons") {
  // parse bin file
  std::vector<double> fPtBinEdgesv, fEtaBinEdgesv, fPhiBinEdgesv, fNPVBinEdgesv, fJetsBinEdgesv, fMETBinEdgesv;
  // flags for |eta| and |phi| binning
  bool fDoAbsEta, fDoAbsPhi;
  // flags for binnings to compute efficiencies for
  bool fDoPt, fDoEta, fDoPhi, fDoEtaPt, fDoEtaPhi, fDoNPV, fDoJets, fDoMET;
  
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
  std::ofstream ofs; ofs.open(Form("%s/jobArgs.txt",outDir.c_str())); assert(ofs.is_open());
  for(unsigned iPt=0;iPt<fPtBinEdgesv.size()-1; iPt++) { for(unsigned iEta=0; iEta<fEtaBinEdgesv.size()-1; iEta++) {
    TString titleStringPass,titleStringFail;
    if(flavor=="muons") {
      titleStringPass=Form("Z#rightarrow#mu#mu, passing probes p_{T}[%d,%d] #eta[%.4f,%.4f]",(int)fPtBinEdgesv[iPt],(int)fPtBinEdgesv[iPt+1],fEtaBinEdgesv[iEta],fEtaBinEdgesv[iEta+1]);
      titleStringFail=Form("Z#rightarrow#mu#mu, failing probes p_{T}[%d,%d] #eta[%.4f,%.4f]",(int)fPtBinEdgesv[iPt],(int)fPtBinEdgesv[iPt+1],fEtaBinEdgesv[iEta],fEtaBinEdgesv[iEta+1]);
    } else {
      titleStringPass=Form("Z#rightarrowee, passing probes p_{T}[%d,%d] #eta[%.4f,%.4f]",(int)fPtBinEdgesv[iPt],(int)fPtBinEdgesv[iPt+1],fEtaBinEdgesv[iEta],fEtaBinEdgesv[iEta+1]);
      titleStringFail=Form("Z#rightarrowee, failing probes p_{T}[%d,%d] #eta[%.4f,%.4f]",(int)fPtBinEdgesv[iPt],(int)fPtBinEdgesv[iPt+1],fEtaBinEdgesv[iEta],fEtaBinEdgesv[iEta+1]);
    }
    titleStringPass.ReplaceAll(" ","~");
    titleStringFail.ReplaceAll(" ","~");
    ofs
      << outDir.c_str() 
      << Form(" pass_ptBin%d_etaBin%d",iPt,iEta)
      << ((flavor=="muons")? Form(" %s/SingleMuon2016_BaselineToMedium_RochesterCorrections_muonTnP_binnedHistos.root",outDir.c_str())                      : Form(" %s/SingleElectron2016_BaselineToMedium_ScaleSmearCorrections_electronTnP_binnedHistos.root",outDir.c_str()))
      << ((flavor=="muons")? Form(" %s/DYJetsToLL_BaselineToMedium_RochesterCorrections_PtBinnedPlusInclusiveNLO_muonTnP_binnedHistos.root",outDir.c_str()) : Form(" %s/DYJetsToLL_BaselineToMedium_ScaleSmearCorrections_PtBinnedPlusInclusiveNLO_electronTnP_binnedHistos.root",outDir.c_str()))
      << " fitterShape::kTemplateConvGaus fitterShape::kBkgDasPlusExp " 
      << titleStringPass.Data() 
      << std::endl;
    ofs
      << outDir.c_str() 
      << Form(" fail_ptBin%d_etaBin%d",iPt,iEta)
      << ((flavor=="muons")? Form(" %s/SingleMuon2016_BaselineToMedium_RochesterCorrections_muonTnP_binnedHistos.root",outDir.c_str())                      : Form(" %s/SingleElectron2016_BaselineToMedium_ScaleSmearCorrections_electronTnP_binnedHistos.root",outDir.c_str()))
      << ((flavor=="muons")? Form(" %s/DYJetsToLL_BaselineToMedium_RochesterCorrections_PtBinnedPlusInclusiveNLO_muonTnP_binnedHistos.root",outDir.c_str()) : Form(" %s/DYJetsToLL_BaselineToMedium_ScaleSmearCorrections_PtBinnedPlusInclusiveNLO_electronTnP_binnedHistos.root",outDir.c_str()))
      << " fitterShape::kTemplateConvGaus fitterShape::kBkgDasPlusExp " 
      << titleStringFail.Data() 
      << std::endl;
  }}


}
