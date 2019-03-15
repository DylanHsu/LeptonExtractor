#include <TH1D.h>
#include <TTree.h>
#include <TFile.h>
#include <TString.h>
#include <TSystem.h>
#include <TLorentzVector.h>
#include <cassert>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <unistd.h>
#include "MuonPogFitterTree.C"

//const int nMassBins=200;
const int year=2016;
const int nMassBins=35;
const float xmin=60;
const float xmax=130;
//const bool doBtoF=true;
//const bool doGtoH=false;
enum massSelType {
  kZll,           //0
  kEM,            //1
  kZeeSameSign,   //2
  kZmmSameSign,   //3
  kMuPi,          //4
  kEPi,           //5
  kZMuMuFromTrack,//6
  kZMuPiFromTrack,//7
  kGenZll,        //8
  kSingleMuTrig,  //9
  knSelTypes      
};
TString dirPath = TString(gSystem->Getenv("CMSSW_BASE")) + "/src/";

void binSplitterMuonPog( // Muon POG AOD Bin Splitter
  TString inputFileName, 
  TString outputFileName, 
  std::string binFile, 
  double scalefactor=1, 
  bool isMC=true,
  massSelType selection = kZll,
  bool alt_tag=false,
  int muonTestSel=5, 
  // 0-global or tracker
  // 4-loose ID&iso
  // 5-medium2016ID+tight iso
  // 6-tightID+tight iso
  // 100- single muon trigger given tight muon ID/iso
  TString pileupProfile="",
  TString reweightBandFileName=""

) {
  if(inputFileName==outputFileName) { printf("check output file name!\n"); return; }
  if(!(muonTestSel==0 || muonTestSel==4 || muonTestSel==5 || muonTestSel==6
  || muonTestSel==7 || muonTestSel==100)) { printf("invalid muonTestSel (0|4|5|6)\n"); return; }
  if(selection!=kZll && selection!=kGenZll && selection!=kMuPi && selection!=kZmmSameSign && selection!=kZMuMuFromTrack && selection!=kZMuPiFromTrack) { 
    printf("selection argument not supported by AOD module\n");
    return; 
  }
  // Load pileup weights
  bool doBtoF=false;
  bool doGtoH=false;
  TH1D *puWeights=0; {
    if(pileupProfile=="BtoF")      doBtoF=true;
    else if(pileupProfile=="GtoH") doGtoH=true;
    TFile *puFile = 0;
    if(year==2016) {
      if(doBtoF)      puFile = TFile::Open(Form("%sLeptonExtractor/puWeights_2016_bf.root", dirPath.Data()), "READ");
      else if(doGtoH) puFile = TFile::Open(Form("%sLeptonExtractor/puWeights_2016_gh.root", dirPath.Data()), "READ");
      else            puFile = TFile::Open(Form("%sLeptonExtractor/puWeights_80x_37ifb.root", dirPath.Data()), "READ");
    } else {
      puFile = TFile::Open(Form("%sLeptonExtractor/puWeights_90x.root", dirPath.Data()), "READ");
    }
    puWeights = (TH1D*)puFile->Get("puWeights");
    assert(puWeights);
    puWeights->SetDirectory(0);
    puFile->Close();  
  }
  // Load FSR reweight band file, but we don't load all the histograms into memory
  TFile *reweightBandFile=0; bool useReweightBand=false;
  if(reweightBandFileName!="" && isMC) {
    reweightBandFile = TFile::Open(reweightBandFileName.Data());
    assert(reweightBandFile);
    useReweightBand=true;
  }

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


  TFile *inputFile=TFile::Open(inputFileName,"READ"); assert(inputFile);
  
  TTree *inputTree=0;
  if(selection==kMuPi || selection==kZMuPiFromTrack || selection==kZmmSameSign)
    inputTree=(TTree*)inputFile->Get("tpTreeSameSign/fitter_tree");
  else
    inputTree=(TTree*)inputFile->Get("tpTree/fitter_tree");

  vector<TString> mpftBranches = {
    "run", "lumi", "event",
    "mass", "pair_probeMultiplicity", "tag_nVertices",
    "pt", "eta", "abseta", "Loose", "Medium2016", "Medium", "Tight2012", "combRelIsoPF04dBeta",
    "tag_combRelIsoPF04dBeta", "tag_pt", "tag_eta", "tag_abseta", "tag_IsoTkMu24", "tag_IsoMu24",
    "TM","Glb","PF"
  };
  MuonPogFitterTree mpft(inputTree, isMC);
  mpft.fChain->SetBranchStatus("*",0);
  if(isMC) { 
    mpft.fChain->SetBranchStatus("mc*", 1); 
    mpft.fChain->SetBranchStatus("pair_genWeight", 1); 
  }
  // need pileup reweighting
  for(unsigned i=0; i<mpftBranches.size(); i++) mpft.fChain->SetBranchStatus(mpftBranches[i],1);

  TFile *outputFile = TFile::Open(outputFileName,"RECREATE");
  TH1D *histosPass[2048], *histosFail[2048];
  assert(outputFile && outputFile->IsOpen());
  int iHisto=0;
  for(unsigned iPt=0;iPt<fPtBinEdgesv.size()-1; iPt++) { for(unsigned iEta=0; iEta<fEtaBinEdgesv.size()-1; iEta++) {
    char histName[256];
    sprintf(histName,"pass_ptBin%d_etaBin%d", iPt, iEta);
    histosPass[iHisto] = new TH1D(histName,histName,nMassBins,xmin,xmax);
    histosPass[iHisto]->Sumw2();
    sprintf(histName,"fail_ptBin%d_etaBin%d", iPt, iEta);
    histosFail[iHisto] = new TH1D(histName,histName,nMassBins,xmin,xmax);
    histosFail[iHisto]->Sumw2();
    iHisto++;
  }}
  int iHistoMax=iHisto-1;
  
  for(unsigned int ientry=0; ientry<(unsigned int)inputTree->GetEntries(); ientry++) {
    if(ientry%1000000==0) printf("reading entry %d/%d\n",ientry,(unsigned int)inputTree->GetEntries());

    mpft.fChain->GetBranch("run")->GetEntry(ientry);
    if(!isMC) {
      if(doBtoF && mpft.run>=278803) continue; // before HIP problem fixed
      if(doGtoH && mpft.run<278803) continue; // after HIP problem fixed
    }
    mpft.fChain->GetBranch("mass")->GetEntry(ientry);
    float mass=mpft.mass;
    if(mass<xmin || mass>=xmax) continue;
    mpft.fChain->GetEntry(ientry);

    bool passNominalTag = mpft.tag_pt>26 && mpft.tag_combRelIsoPF04dBeta < 0.15 && mpft.tag_combRelIsoPF04dBeta > -0.5 && (mpft.tag_IsoMu24==1 || mpft.tag_IsoTkMu24==1);
    bool passAltTag = mpft.tag_pt>31 && mpft.tag_combRelIsoPF04dBeta < 0.2 && mpft.tag_combRelIsoPF04dBeta > -0.5 && (mpft.tag_IsoMu24==1 || mpft.tag_IsoTkMu24==1);
    bool probeAcceptance = mpft.abseta<=2.4 && mpft.pt >= 10;
    bool probeIsMuon = mpft.TM!=0 && mpft.Glb!=0;
    bool probeIsTightMuon = mpft.Tight2012 && mpft.combRelIsoPF04dBeta < 0.15 && mpft.combRelIsoPF04dBeta > -0.5;
    bool probePassTestSel = 
      (muonTestSel==0 && probeIsMuon) ||
      (muonTestSel==4 && mpft.Loose && mpft.combRelIsoPF04dBeta < 0.2 && mpft.combRelIsoPF04dBeta > -0.5) ||
      (muonTestSel==5 && mpft.Medium2016 && mpft.combRelIsoPF04dBeta < 0.15 && mpft.combRelIsoPF04dBeta > -0.5) ||
      (muonTestSel==6 && mpft.Tight2012 && mpft.combRelIsoPF04dBeta < 0.15 && mpft.combRelIsoPF04dBeta > -0.5) || 
      (muonTestSel==7 && mpft.Tight2012 && mpft.combRelIsoPF04dBeta < 0.06 && mpft.combRelIsoPF04dBeta > -0.5) ||
      (muonTestSel==100 && probeIsTightMuon && (mpft.IsoMu24==1 || mpft.IsoTkMu24==1));

    bool passMult = mpft.pair_probeMultiplicity==1;
    bool passTag = (!alt_tag && passNominalTag) || (alt_tag && passAltTag);

    bool passSelection=false;
    if     (selection==           kZll ) passSelection=passTag && probeIsMuon && probeAcceptance;
    else if(selection==        kGenZll ) passSelection=passTag && probeIsMuon && probeAcceptance;
    else if(selection==   kZmmSameSign ) passSelection=passTag && probeIsMuon && probeAcceptance;
    else if(selection==          kMuPi ) passSelection=passTag && probeIsMuon && probeAcceptance;
    else if(selection== kZMuPiFromTrack) passSelection=passTag && passMult && probeAcceptance;
    else if(selection== kZMuMuFromTrack) passSelection=passTag && passMult && probeAcceptance;
    else if(selection==   kSingleMuTrig) passSelection=passTag && probeIsTightMuon;
    else                                 passSelection=false;
    //if(!passSelection) printf("failing selection, passTag=%d, probeIsMuon=%d, probeAcceptance=%d\n",passTag,probeIsMuon,probeAcceptance);
    if(!passSelection) continue;
    
    if((selection == kZMuMuFromTrack || selection==kZMuPiFromTrack) && isMC && mpft.mcTrue==0) continue;
    
    iHisto=0;
    bool ignoreEtaHighEnergy=(selection==kMuPi);
    bool ignorePassFlag=(selection==kMuPi || selection==kZMuPiFromTrack);

    for(unsigned iPt=0;iPt<fPtBinEdgesv.size()-1; iPt++) { for(unsigned iEta=0; iEta<fEtaBinEdgesv.size()-1; iEta++) {
      double pt=mpft.pt;
      if(pt >= fPtBinEdgesv[iPt] && pt < fPtBinEdgesv[iPt+1]) { 
      double eta = fDoAbsEta? mpft.abseta:mpft.eta;
      bool passEta=(eta >= fEtaBinEdgesv[iEta] && eta < fEtaBinEdgesv[iEta+1]) || (pt >= 50 && ignoreEtaHighEnergy);
        if(passEta) {
          double weight;
          if(isMC) weight = mpft.pair_genWeight * puWeights->GetBinContent(mpft.tag_nVertices);
          else     weight = 1;
          //if(isMC) printf("weight %f, pair_genWeight %f, tag_nVertices %f\n", weight, mpft.pair_genWeight,mpft.tag_nVertices);
          if(selection==kGenZll) mass=mpft.mcMass;
          if(useReweightBand) { // apply mass reweighting
            TH1D *theReweightBand = (TH1D*)reweightBandFile->Get(Form("reweightBand_ptBin%d_etaBin%d",iPt,iEta));
            assert(theReweightBand);
            int m = theReweightBand->FindBin(mpft.mcMass);
            if(m>=1 && m<=theReweightBand->GetNbinsX())
              weight *= theReweightBand->GetBinContent(m);
          }
          if     (probePassTestSel  || ignorePassFlag) histosPass[iHisto]->Fill(mass,weight);
          if     (!probePassTestSel || ignorePassFlag) histosFail[iHisto]->Fill(mass,weight);
        }
      }
      iHisto++;
    }}
  }
  for(iHisto=0; iHisto <= iHistoMax; iHisto++) histosPass[iHisto]->Write();
  for(iHisto=0; iHisto <= iHistoMax; iHisto++) histosFail[iHisto]->Write();
  outputFile->Close();
}
