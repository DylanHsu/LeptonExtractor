#include <map>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/lexical_cast.hpp>
#include "PandaTree/Objects/interface/Event.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include <sstream>
#include "RoccoR.cc"
#include <TRandom3.h>

const float massForArbitration=91.1876;
const int year=2017;

struct leptonVariables {
  int tag_id            = 6; 
  int tag_iso           = 6; 
  int probe_id          = 0; 
  int probe_iso         = 0; 
  int passing_probe_id  = 6; 
  int passing_probe_iso = 6; 
  double tag_pt_min = 30;
  double tag_eta_max = 2.1;
  long double xsec;
} gettagprobe;
// Above are the cut variables.
// Baseline cuts determined by probe_* (0 = default cut). 
// Tag cuts determined by tag_* (6 = Tight, 5 = Medium, 4 = Loose).
// Probe cuts determined by passing_probe_* (6 = Tight, 5 = Medium, 4 = Loose).

TString dirPath = TString(gSystem->Getenv("CMSSW_BASE")) + "/src/";

typedef std::map<UInt_t,std::vector<std::pair <UInt_t, UInt_t> > > MapType;
//string jsonFile = "certs/Cert_294927-300575_13TeV_PromptReco_Collisions17_JSON.txt";

void animator(int iEntry, int nEntries);
bool selector(panda::Electron const&, int, int);
bool selector(panda::Muon const&, int, int);
bool passJetId(panda::Jet const&, Float_t fMVACut[4][4]);
void InitializeJetIdCuts(Float_t fMVACut[4][4]);
vector<int> passProbe(panda::Electron const&,leptonVariables gettagprobe, bool);
vector<int> passProbe(panda::Muon const&,leptonVariables gettagprobe, bool);

void make_tnp_skim(
                   string input_file_name,
                   string output_basename,
                   bool do_electrons = true,  //perform pair matching on electrons 
                   bool do_muons = true,      //perform pair matching on muons
                   bool real_data = false,     //state whether data or Monte Carlo
                   bool verbose = false,      //DEBUG
                   bool truth_matching = true,//used for Monte Carlo truth matching
                   int max_entries = -1,      //DEBUG (max entries in fileset to be run)
                     int electron_trigger=3,
                   int muon_trigger=6,
                   double truth_matching_dR = 0.3){//max Delta-R for truth matching 
  string jsonFile;
  if(year==2016) jsonFile = Form("%sLeptonExtractor/certs/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt",dirPath.Data());
  else if(year==2017) jsonFile = Form("%sLeptonExtractor/certs/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON_v1.txt",dirPath.Data());
  gSystem->Load("libPandaTreeObjects.so"); 
  
  //OUTPUT FILE:

  //declare output variables
  unsigned out_runNum, // event ID
    out_lumiSec,
    out_evtNum,
    out_npv, // number of primary vertices
    pass; // whether probe passes requirements
  float         npu=1;                     // mean number of expected pileup
  float         scale1fb=1;                // event weight per 1/fb
  float         mass;                      // tag-probe mass
  char          qtag, qprobe;              // tag, probe charge
  float         truth_tag, truth_probe;    // tag, probe truth
  float         met;                       // missing ET
  unsigned char njets;                     // number of jets
  short         tagPid ,probePid;          // particle ID
  unsigned char probeMultiplicity;
  bool          bestMass;
  bool          probeFromHardestPV;
  TLorentzVector *p4_tag=0, *p4_probe=0;        // tag, probe 4-vector 
  TLorentzVector *genp4_tag=0, *genp4_probe=0;        // tag, probe 4-vector at generator level

  //Output filename setup
  string output_dir = "";
  string electron_filename = output_dir +output_basename +"_electronTnP.root";
  string muon_filename     = output_dir +output_basename +"_muonTnP.root";


  //Output file initialization
  TFile *electron_outfile = nullptr;
  TTree *electron_pair_tree = nullptr;
  TFile *muon_outfile = nullptr;
  TTree *muon_pair_tree = nullptr;

  
  //Grow the branches of the tree
  if(do_electrons) {
    electron_outfile = TFile::Open(electron_filename.c_str(),"RECREATE");
    electron_pair_tree = new TTree("events", "Electron skim for TnP script");
    electron_pair_tree->Branch("runNum",   &out_runNum,   "runNum/i"   );  
    electron_pair_tree->Branch("lumiSec",  &out_lumiSec,  "lumiSec/i"  );  
    electron_pair_tree->Branch("evtNum",   &out_evtNum,   "evtNum/i"   );  
    electron_pair_tree->Branch("npv",      &out_npv,      "npv/i"      );  
    electron_pair_tree->Branch("pass",     &pass,     "pass/i"     );  
    electron_pair_tree->Branch("npu",      &npu,      "npu/F"      );  
    electron_pair_tree->Branch("scale1fb", &scale1fb, "scale1fb/F" );
    electron_pair_tree->Branch("mass",     &mass,     "mass/F"     );  
    electron_pair_tree->Branch("qtag",     &qtag,     "qtag/B"     );  
    electron_pair_tree->Branch("qprobe",   &qprobe,   "qprobe/B"   );  
    electron_pair_tree->Branch("njets",    &njets,    "njets/b"   );  
    electron_pair_tree->Branch("met",      &met,      "met/F"   );  
    electron_pair_tree->Branch("tag",   "TLorentzVector", &p4_tag   );  
    electron_pair_tree->Branch("probe", "TLorentzVector", &p4_probe );          
    electron_pair_tree->Branch("tagPid",      &tagPid,      "tagPid/S"   );  
    electron_pair_tree->Branch("probePid",      &probePid,      "probePid/S"   );  
    electron_pair_tree->Branch("probeMultiplicity",    &probeMultiplicity,    "probeMultiplicity/b"   );  
    electron_pair_tree->Branch("bestMass",             &bestMass         ,    "bestMass/O"            );  
    electron_pair_tree->Branch("probeFromHardestPV",   &probeFromHardestPV,   "probeFromHardestPV/O"  );  
    if(!real_data) {
      electron_pair_tree->Branch("truth_tag",     &truth_tag,     "truth_tag/F"     );  
      electron_pair_tree->Branch("truth_probe",   &truth_probe,   "truth_probe/F"   );
      electron_pair_tree->Branch("genTag",   "TLorentzVector", &genp4_tag   );  
      electron_pair_tree->Branch("genProbe", "TLorentzVector", &genp4_probe );          
    }
  }
  

  if(do_muons) {
    muon_outfile = TFile::Open(muon_filename.c_str(),"RECREATE");
    muon_pair_tree = new TTree("events", "Muon skim for TnP script");
    muon_pair_tree->Branch("runNum",   &out_runNum,   "runNum/i"   );  
    muon_pair_tree->Branch("lumiSec",  &out_lumiSec,  "lumiSec/i"  );  
    muon_pair_tree->Branch("evtNum",   &out_evtNum,   "evtNum/i"   );  
    muon_pair_tree->Branch("npv",      &out_npv,      "npv/i"      );  
    muon_pair_tree->Branch("pass",     &pass,     "pass/i"     );  
    muon_pair_tree->Branch("npu",      &npu,      "npu/F"      );  
    muon_pair_tree->Branch("scale1fb", &scale1fb, "scale1fb/F" );
    muon_pair_tree->Branch("mass",     &mass,     "mass/F"     );  
    muon_pair_tree->Branch("qtag",     &qtag,     "qtag/B"     );  
    muon_pair_tree->Branch("qprobe",   &qprobe,   "qprobe/B"   );  
    muon_pair_tree->Branch("tag",   "TLorentzVector", &p4_tag   );  
    muon_pair_tree->Branch("probe", "TLorentzVector", &p4_probe );      
    muon_pair_tree->Branch("njets",    &njets,    "njets/b"   );  
    muon_pair_tree->Branch("met",      &met,      "met/F"   );  
    muon_pair_tree->Branch("tagPid",   &tagPid,   "tagPid/S"   );  
    muon_pair_tree->Branch("probePid", &probePid, "probePid/S"   );  
    muon_pair_tree->Branch("probeMultiplicity",    &probeMultiplicity,    "probeMultiplicity/b"   );  
    muon_pair_tree->Branch("bestMass",             &bestMass         ,    "bestMass/O"            );  
    muon_pair_tree->Branch("probeFromHardestPV",   &probeFromHardestPV,   "probeFromHardestPV/O"  );  
    if(!real_data) {
      muon_pair_tree->Branch("genTag",   "TLorentzVector", &genp4_tag   );  
      muon_pair_tree->Branch("genProbe", "TLorentzVector", &genp4_probe );      
      muon_pair_tree->Branch("truth_tag",     &truth_tag,     "truth_tag/F"     );  
      muon_pair_tree->Branch("truth_probe",   &truth_probe,   "truth_probe/F"   );  
    }
  }
  p4_tag=new TLorentzVector; p4_probe=new TLorentzVector;
  genp4_tag=new TLorentzVector; genp4_probe=new TLorentzVector;

  //////////////////////////////////////////////////////////////////////////////////////////
  
  //Load pileup corrections
  //TFile *puFile = TFile::Open(Form("%sLeptonExtractor/puWeights_80x_37ifb.root", dirPath.Data()), "READ");
  //TFile *puFile = TFile::Open(Form("%sLeptonExtractor/puWeights_2016_bf.root", dirPath.Data()), "READ");
  //TFile *puFile = TFile::Open(Form("%sLeptonExtractor/puWeights_2016_gh.root", dirPath.Data()), "READ");
  //TH1D *puWeights = (TH1D*)puFile->Get("puWeights"); puWeights->SetDirectory(0);
  //puFile->Close();  

  //INPUT FILE:
  //Setting up Tree files and various input variables  
  printf("Trying to open file %s...\n", input_file_name.c_str());
  TFile *input_file;
  int retries=0;
  while(true) {
    input_file = TFile::Open(input_file_name.c_str(),"read");
    if(input_file && input_file->IsOpen()) break;
    retries++;
    if(retries>100) { throw std::runtime_error("Error opening input file"); return; }
  }
  //TFile*  input_file=TFile::Open(input_file_name.c_str(),"READ");
  //if(!input_file || !input_file->IsOpen()) {
  //  printf("\"I couldn't open it\"~Alex Jones\n");
  //  assert(0); return;
  //}
  TTree* tree = (TTree*)input_file->Get("events"); // get the tree object from the file
  panda::Event event; // create an Event object
  event.setStatus(*tree, {"!*"});
  event.setAddress(*tree, {"runNumber", "lumiNumber", "eventNumber", "muons", "electrons", "npv", "npvTrue", "genParticles","isData","pfMet","chsAK4Jets","pfCandidates","vertices","weight","tracks","vertices", "triggers","triggerObjects"}); 
  
  event.run.setLoadTrigger(true);
  TString ele32Filter1="hltEle32L1DoubleEGWPTightGsfTrackIsoFilter";
  TString ele32Filter2="hltEGL1SingleEGOrFilter";

  Long64_t nEntries = tree->GetEntries();
  Long64_t sum_mc_weights=0;  
  Float_t fMVACut[4][4];
  InitializeJetIdCuts(fMVACut);
  int min_runNum=99999999, max_runNum=0;
  
  if(!real_data) {
    TH1D* all_tree=(TH1D*)input_file->FindObjectAny("hSumW");
    sum_mc_weights = all_tree->GetBinContent(1);
  
  }
  
  //cout << input_file_name;
  
  //////////////////////////////////////////////////////////////////////////////////////////
  
  //Read json file into boost property tree
  MapType fMap;
  boost::property_tree::ptree jsonTree;
  boost::property_tree::read_json(jsonFile.c_str(),jsonTree);
  
  //Loop through boost property tree and fill the MapType structure with the list of good lumi
  //ranges for each run
  for (boost::property_tree::ptree::const_iterator it = jsonTree.begin(); it!=jsonTree.end(); ++it) {
    UInt_t runNumber = boost::lexical_cast<UInt_t>(it->first);
    MapType::mapped_type &lumiPairList = fMap[runNumber];
    boost::property_tree::ptree lumiPairListTree = it->second;
    for (boost::property_tree::ptree::const_iterator jt = lumiPairListTree.begin(); jt!=lumiPairListTree.end(); ++jt) {
      boost::property_tree::ptree lumiPairTree = jt->second;
      if (lumiPairTree.size()==2) {
        UInt_t firstLumi = boost::lexical_cast<UInt_t>(lumiPairTree.begin()->second.data());
        UInt_t lastLumi = boost::lexical_cast<UInt_t>((++lumiPairTree.begin())->second.data());
        lumiPairList.push_back(std::pair<UInt_t,UInt_t>(firstLumi,lastLumi));
      }
    }
  }

  //Load rochester corrections
  // Initialize the random seed based on Dylan's age in seconds
  std::time_t t = std::time(0);
  unsigned long int time_now = static_cast<unsigned long int>(time(NULL));
  TRandom3 rng(time_now-731178000);
  RoccoR rochesterCorrection(Form("%sLeptonExtractor/rcdata.2016.v3",dirPath.Data()));
  
  //ACTIVE FUNCTION BEGINS:

  //Beginning of loop over entries:    

  if(do_electrons)
    printf("Electrons are go!\n");
  if(do_muons)
    printf("Muons are go!\n");
  printf("Parsing through events in fileset...\n");
  for (Long64_t iEntry = 0; iEntry != nEntries && iEntry != max_entries; ++iEntry) {
    event.getEntry(*tree, iEntry);
    if (!verbose)
      animator(iEntry,nEntries);

    if (verbose)
      printf("Begin event number #:%lld/%lld\n", iEntry,nEntries);
    
    // Check data certification
    bool certifiedEvent=false;
    std::pair<unsigned, unsigned> runLumi(event.runNumber, event.lumiNumber);      
    MapType::const_iterator it = fMap.find(runLumi.first);
    if (it!=fMap.end()) {
      //check lumis
      const MapType::mapped_type &lumiPairList = it->second;
      for (MapType::mapped_type::const_iterator jt = lumiPairList.begin(); jt<lumiPairList.end(); ++jt) {
        if (runLumi.second >= jt->first && runLumi.second <= jt->second) {
          //found lumi in accepted range
          certifiedEvent=true;
        }
      }
    }
    if(real_data && !certifiedEvent) { if(verbose) printf("failed golden json\n"); continue; }
 
    //Filling event-specific branches
    out_runNum=event.runNumber;
    out_lumiSec=event.lumiNumber;
    out_evtNum=event.eventNumber;
    out_npv=event.npv;
    npu=event.npvTrue;
    double eMC_weights = event.weight;
    //    printf("%d", eMC_weights);
    


    
    //Looping over electrons in event to identify which pass tag criteria
    unsigned goodIsTight = 0;
    if(event.electrons.size() != 32) for (unsigned iE = 0; iE != event.electrons.size(); ++iE) {
      auto& electron = event.electrons[iE];
      if (selector(electron, gettagprobe.tag_id, gettagprobe.tag_iso)){
        goodIsTight++;
      }
    }
    
    //Looping over muons in event to identify which pass tag criteria
    if(event.muons.size() != 32) for (unsigned iM = 0; iM != event.muons.size(); ++iM) {
        //      cout<<iM<<endl;
      auto& muon = event.muons[iM];
      //      cout<<event.eventNumber<<endl;
      //cout<<muon.pt()<<endl;

      if (selector(muon, gettagprobe.tag_id, gettagprobe.tag_iso)){
        goodIsTight++;
      }
    }

    if(goodIsTight==0) continue;
    //Reaching this point == pass 


    //Checking delta-R and compiling lists of Leptons
    vector<int> idJet;
    for (unsigned nJ = 0; nJ != event.chsAK4Jets.size(); ++nJ){
      auto& jet = event.chsAK4Jets[nJ];

      //Passing criteria for Jet Pt
      if(jet.pt()<30) continue;
      
      //Determine if Jet is actually electron: if not, add to list of Jets
      bool isElectron = false; bool isMuon = false;
      for(unsigned iE = 0; iE != event.electrons.size(); ++iE){
        auto& electron  = event.electrons[iE];
        if (electron.loose){
          if(jet.dR(electron) < 0.16) { isElectron = true; break; }
        }
      }
      if (isElectron == true) continue;
        
      //Determine if Jet is actually muon: if not, add to list of Jets
      for(unsigned iM = 0; iM != event.muons.size(); ++iM){
        auto& muon  = event.muons[iM];
        if(muon.loose){
          if(jet.dR(muon) < 0.16) { isMuon = true; break; }
        }
      }
      if (isMuon == true) continue;
      idJet.push_back(nJ);
    }
    //return the number of Jets:
    njets = idJet.size();
    met = event.pfMet.pt;
      

    //PROBE:
    
    //Initialize vectors to hold lepton information:
    std::vector<TLorentzVector> p4_ele_tag_, p4_ele_passing_probe_, p4_ele_failing_probe_, p4_mu_tag_, p4_mu_passing_probe_, p4_mu_failing_probe_;
    std::vector<TLorentzVector> genp4_ele_tag_, genp4_ele_passing_probe_, genp4_ele_failing_probe_, genp4_mu_tag_, genp4_mu_passing_probe_, genp4_mu_failing_probe_;
    std::vector<int> q_ele_tag_, q_ele_passing_probe_, q_ele_failing_probe_, q_mu_tag_, q_mu_passing_probe_, q_mu_failing_probe_;
    std::vector<double>truth_ele_tag_,truth_ele_failing_probe_,truth_ele_passing_probe_,truth_mu_tag_,truth_mu_passing_probe_,truth_mu_failing_probe_;
    std::vector<float> vz_mu_tag_, vz_ele_tag_;

    // Loop over Electrons
    if (event.electrons.size()>0) { for (unsigned iE = 0; iE != event.electrons.size(); ++iE) {
      auto& electron = event.electrons[iE];
      if(electron.smearedPt < 10.) continue;
      bool electron_trigger_matched = false;
      double truth = 0;
      int charge = electron.charge;
      bool pass_tag_trigger = false;
      float vz = -99;
      if(electron.vertex.isValid()) vz = electron.vertex.get()->z;
      double dR;
      TLorentzVector genP4;
      if (!real_data && truth_matching){
        for (unsigned iG = 0; iG != event.genParticles.size(); ++iG){
          auto& genParticle = event.genParticles[iG];
          if (genParticle.finalState != 1)
            continue;
          TVector3 genParticle3;
          genParticle3.SetPtEtaPhi(genParticle.pt(),genParticle.eta(),genParticle.phi());
          dR = genParticle3.DeltaR(electron.p4().Vect());
          if (TMath::Abs(genParticle.pdgid) == 11 && dR < truth_matching_dR) {
            truth = TMath::Max(0.001, dR);
            genP4.SetPtEtaPhiM(genParticle.pt(),genParticle.eta(),genParticle.phi(),511e-6);
            break;
          }
        }
        if(truth==0) continue;
      }

      if(event.isData){
        //            pass_tag_trigger = (electron.triggerMatch[panda::Electron::fHLT_Ele27_eta2p1_WPLoose_Gsf]);
        if(year==2016) pass_tag_trigger = (electron.triggerMatch[panda::Electron::fEl27Tight]);
        else if(year==2017) {
          //pass_tag_trigger = (electron.triggerMatch[panda::Electron::fEl35Tight]); 
          pass_tag_trigger=false;
          panda::HLTObjectStore::HLTObjectVector filter1Objects, filter2Objects;
          filter1Objects=event.triggerObjects.filterObjects(ele32Filter1.Data());
          filter2Objects=event.triggerObjects.filterObjects(ele32Filter2.Data());
          if(filter1Objects.size()>0 && filter2Objects.size()>0) {
            TLorentzVector filter1ObjectP4, filter2ObjectP4;
            filter1ObjectP4.SetPtEtaPhiM(filter1Objects[0]->pt(), filter1Objects[0]->eta(), filter1Objects[0]->phi(), filter1Objects[0]->m());
            filter2ObjectP4.SetPtEtaPhiM(filter2Objects[0]->pt(), filter2Objects[0]->eta(), filter2Objects[0]->phi(), filter2Objects[0]->m());

            if(filter1ObjectP4.DeltaR(electron.p4())<0.1 && filter2ObjectP4.DeltaR(electron.p4())<0.1)
              pass_tag_trigger=true;
          }
        }
      } else{
        pass_tag_trigger = true;
      }
      TLorentzVector lepP4; lepP4.SetPtEtaPhiM(electron.smearedPt, electron.eta(), electron.phi(), 511e-6);
      vector<int> passing_ids = passProbe(electron, gettagprobe, pass_tag_trigger);
      for (unsigned long i = 0; i < passing_ids.size(); i++){
        switch(passing_ids[i]){
        case 1:
          //Fails Probe
          p4_ele_failing_probe_.push_back(lepP4);
          genp4_ele_failing_probe_.push_back(genP4);
          q_ele_failing_probe_.push_back(charge);
          truth_ele_failing_probe_.push_back(truth);
          break;
        case 2:
          //Passes Probe
          p4_ele_passing_probe_.push_back(lepP4);
          genp4_ele_passing_probe_.push_back(genP4);
          q_ele_passing_probe_.push_back(charge);
          truth_ele_passing_probe_.push_back(truth);
          break;
        case 3:
          //Passes Tag
          p4_ele_tag_.push_back(lepP4);
          genp4_ele_tag_.push_back(genP4);
          q_ele_tag_.push_back(charge);
          truth_ele_tag_.push_back(truth);
          vz_ele_tag_.push_back(vz);
          break;
        default:
          break;
        }
      }
    }}
    if (verbose){
      printf("Number of vectors in electron tag list:%lu\n", p4_ele_tag_.size());
      printf("Number of vectors in electron pass probe list:%lu\n", p4_ele_passing_probe_.size());
    }

    // Loop over Muons
    if (event.muons.size()>0) { for (unsigned iM = 0; iM != event.muons.size(); ++iM) {
      auto& muon = event.muons[iM];
      if(muon.pt()<1.) continue;
      bool muon_trigger_matched = false;
      int charge = muon.charge;
      double truth;
      bool pass_tag_trigger = false;
      double ptCorrection=1;
      float vz = -99;
      if(muon.vertex.isValid()) vz = muon.vertex.get()->z;
      double dR;
      TLorentzVector genP4;
      if (!real_data && truth_matching){
        for (unsigned iG = 0; iG != event.genParticles.size(); ++iG){
          auto& genParticle = event.genParticles[iG];
          if (genParticle.finalState != 1)
            continue;
           TVector3 genParticle3;
          genParticle3.SetPtEtaPhi(genParticle.pt(),genParticle.eta(),genParticle.phi());
          dR = genParticle3.DeltaR(muon.p4().Vect());
          if (TMath::Abs(genParticle.pdgid)==13 && dR < truth_matching_dR) {
            truth = TMath::Max(0.001, dR);
            double random1=rng.Rndm();
            genP4.SetPtEtaPhiM(genParticle.pt(),genParticle.eta(),genParticle.phi(),.106);
            ptCorrection=rochesterCorrection.kScaleFromGenMC(charge, muon.pt(), muon.eta(), muon.phi(), muon.trkLayersWithMmt, genParticle.pt(), random1, 0, 0);
            break;
          }
        }
        if(truth==0) continue;
      } else if(!real_data) {
        double random1=rng.Rndm();
        double random2=rng.Rndm();
        ptCorrection=rochesterCorrection.kScaleAndSmearMC(charge, muon.pt(), muon.eta(), muon.phi(), muon.trkLayersWithMmt, random1, random2, 0, 0);
        
      } else {
        ptCorrection=rochesterCorrection.kScaleDT(charge, muon.pt(), muon.eta(), muon.phi(), 0, 0);
      }
      if (event.isData){
        //            pass_tag_trigger = (muon.triggerMatch[panda::Muon::fIsoMu27] || muon.triggerMatch[panda::Muon::fIsoTkMu27]);
        pass_tag_trigger = (muon.triggerMatch[panda::Muon::fIsoMu24] || muon.triggerMatch[panda::Muon::fIsoTkMu24]);
      } else pass_tag_trigger = true;
      double correctedPt = muon.pt()*ptCorrection;
      if(muon.pt()>10. && muon.pt()*ptCorrection < 10. && verbose)
        printf("Rejecting muon due to corrections. (pt,eta,phi)=%f,%f,%f ; trackLayersWithMeasurement=%d\n", muon.pt(), muon.eta(), muon.phi(), muon.trkLayersWithMmt);
      if(correctedPt<10) continue;
      TLorentzVector lepP4; lepP4.SetPtEtaPhiM(correctedPt,muon.eta(),muon.phi(),0.106);
      if(lepP4.Pt()<10) continue;
    
      int i = 0;
      vector<int> passing_ids = passProbe(muon, gettagprobe, pass_tag_trigger);
      for (unsigned long i = 0; i < passing_ids.size(); i++){
        switch(passing_ids[i]){
        case 1: 
          //Fails Probe
          p4_mu_failing_probe_.push_back(lepP4);
          genp4_mu_failing_probe_.push_back(genP4);
          q_mu_failing_probe_.push_back(charge);
          truth_mu_failing_probe_.push_back(truth);
          break;
        case 2:
          //Passes Probe
          p4_mu_passing_probe_.push_back(lepP4);
          genp4_mu_passing_probe_.push_back(genP4);
          q_mu_passing_probe_.push_back(charge);
          truth_mu_passing_probe_.push_back(truth);
          break;
        case 3:
          //Passes Tag
          p4_mu_tag_.push_back(lepP4);
          genp4_mu_tag_.push_back(genP4);
          q_mu_tag_.push_back(charge);
          truth_mu_tag_.push_back(truth);
          vz_mu_tag_.push_back(vz);
          break;
        }
      }
    }}
    if (verbose){
      printf("Number of vectors in muon tag list:%lu\n", p4_mu_tag_.size());
      printf("Number of vectors in muon pass probe list:%lu\n", p4_mu_passing_probe_.size());
    }
    
    if (!real_data) scale1fb = event.weight; // *puWeights->GetBinContent(event.npvTrue);
    
    //ELECTRON PAIR ASSOCIATION:
    //associating electron pairs and filling tree
    unsigned nPassingElePairs=0, nFailingElePairs=0; if(do_electrons) {
      for(unsigned iTag=0; iTag < p4_ele_tag_.size(); iTag++) {
        // passing probes for electrons
        for(unsigned iProbe=0; iProbe < p4_ele_passing_probe_.size(); iProbe++) {
          pass = 1;
          *p4_tag     = p4_ele_tag_[iTag];
          *p4_probe   = p4_ele_passing_probe_[iProbe];
          *genp4_tag     = genp4_ele_tag_[iTag];
          *genp4_probe   = genp4_ele_passing_probe_[iProbe];
          qtag    = q_ele_tag_[iTag];
          qprobe  = q_ele_passing_probe_[iProbe];
          truth_tag    = truth_ele_tag_[iTag];
          truth_probe  = truth_ele_passing_probe_[iProbe];
          tagPid = -11*q_ele_tag_[iTag];
          probePid = -11*q_ele_passing_probe_[iProbe];
          // make sure they aren't the same particle with tiny delta R veto
          if(( p4_tag->DeltaR(*p4_probe) < .0001) ) continue;
          TLorentzVector systemP4 = (*p4_tag) + (*p4_probe);
          mass = systemP4.M();
          if(mass<40.||mass>140.) continue;
          electron_pair_tree->Fill();
          nPassingElePairs++;
        }  
        // failing probes for electrons and filling tree
        for(unsigned iProbe=0; iProbe < p4_ele_failing_probe_.size(); iProbe++) {
          pass = 0;
          *p4_tag     = p4_ele_tag_[iTag];
          *p4_probe   = p4_ele_failing_probe_[iProbe];
          *genp4_tag     = genp4_ele_tag_[iTag];
          *genp4_probe   = genp4_ele_failing_probe_[iProbe];
          qtag    = q_ele_tag_[iTag];
          qprobe  = q_ele_failing_probe_[iProbe];
          truth_tag    = truth_ele_tag_[iTag];
          truth_probe  = truth_ele_failing_probe_[iProbe];
          tagPid = -11*q_ele_tag_[iTag];
          probePid = -11*q_ele_failing_probe_[iProbe];
          // make sure they aren't the same particle with tiny delta R veto
          if(( p4_tag->DeltaR(*p4_probe) < .0001) ) continue;
          TLorentzVector systemP4 = (*p4_tag) + (*p4_probe);
          mass = systemP4.M();          
          if(mass<40.||mass>140.) continue;
          electron_pair_tree->Fill();
          nFailingElePairs++;
        } 
      }
    }
    //MUON PAIR ASSOCIATION:
    // associating muon pairs and filling tree
    unsigned nPassingMuPairs=0, nFailingMuPairs=0; if (do_muons) {
      for (unsigned iTag=0; iTag < p4_mu_tag_.size(); iTag++) {
        // passing probes for muons
        for (unsigned iProbe=0; iProbe < p4_mu_passing_probe_.size(); iProbe++) {
          pass = 1;
          *p4_tag     = p4_mu_tag_[iTag];
          *p4_probe   = p4_mu_passing_probe_[iProbe];
          *genp4_tag     = genp4_mu_tag_[iTag];
          *genp4_probe   = genp4_mu_passing_probe_[iProbe];
          qtag    = q_mu_tag_[iTag];
          qprobe  = q_mu_passing_probe_[iProbe];
          tagPid = -13*q_mu_tag_[iTag];
          probePid = -13*q_mu_passing_probe_[iProbe];
          truth_tag    = truth_mu_tag_[iTag];
          truth_probe  = truth_mu_passing_probe_[iProbe];
          // make sure they aren't the same particle with tiny delta R veto
          if(( p4_tag->DeltaR(*p4_probe) < .0001) ) continue;
          TLorentzVector systemP4 = (*p4_tag) + (*p4_probe);
          mass = systemP4.M();
          if(mass<40.||mass>140.) continue;
          if (verbose) printf("\t\tmade a PASSING muon pair! pTs %f, %f; system mass %f, total charge %d e\n", p4_tag->Pt(), p4_probe->Pt(), mass, qtag+qprobe);
          muon_pair_tree->Fill();
          nPassingMuPairs++;
        }
        
        // failing probes for muons and filling tree
        for(unsigned iProbe=0; iProbe < p4_mu_failing_probe_.size(); iProbe++) {
          pass = 0;
          *p4_tag     = p4_mu_tag_[iTag];  
          *p4_probe   = p4_mu_failing_probe_[iProbe];
          *genp4_tag     = genp4_mu_tag_[iTag];  
          *genp4_probe   = genp4_mu_failing_probe_[iProbe];
          qtag    = q_mu_tag_[iTag];
          qprobe  = q_mu_failing_probe_[iProbe];
          truth_tag    = truth_mu_tag_[iTag];
          truth_probe  = truth_mu_failing_probe_[iProbe];
          tagPid = -13*q_mu_tag_[iTag];
          probePid = -13*q_mu_failing_probe_[iProbe];
          // make sure they aren't the same particle with tiny delta R veto
          if(( p4_tag->DeltaR(*p4_probe) < .0001) ) continue;
          TLorentzVector systemP4 = (*p4_tag) + (*p4_probe);
          mass = systemP4.M();
          if(mass<40.||mass>140.) continue;
          if (verbose) printf("\t\tmade a FAILING muon pair! pTs %f, %f; system mass %f, total charge %d e\n", p4_tag->Pt(), p4_probe->Pt(), mass, qtag+qprobe);
          muon_pair_tree->Fill();
          nFailingMuPairs++;
        } 
      }
    }
    
    // Make e-mu pairs
    if(do_electrons) { // electron tag muon probe pairs
      for(unsigned iTag=0; iTag < p4_ele_tag_.size(); iTag++) {
        // passing probes for muons
        for (unsigned iProbe=0; iProbe < p4_mu_passing_probe_.size(); iProbe++) {
          pass = 1;
          *p4_tag     = p4_ele_tag_[iTag];
          *p4_probe   = p4_mu_passing_probe_[iProbe];
          qtag    = q_ele_tag_[iTag];
          qprobe  = q_mu_passing_probe_[iProbe];
          tagPid = -11*q_ele_tag_[iTag];
          probePid = -13*q_mu_passing_probe_[iProbe];
          truth_tag    = truth_ele_tag_[iTag];
          truth_probe  = truth_mu_passing_probe_[iProbe];
          TLorentzVector systemP4 = (*p4_tag) + (*p4_probe);
          mass = systemP4.M();
          if (verbose) printf("\t\tmade a PASSING e-mu pair! pTs %f, %f; system mass %f, total charge %d e\n", p4_tag->Pt(), p4_probe->Pt(), mass, qtag+qprobe);
          electron_pair_tree->Fill();
        }
        // failing probes for muons
        for(unsigned iProbe=0; iProbe < p4_mu_failing_probe_.size(); iProbe++) {
          pass = 0;
          *p4_tag     = p4_ele_tag_[iTag];  
          *p4_probe   = p4_mu_failing_probe_[iProbe];
          qtag    = q_ele_tag_[iTag];
          qprobe  = q_mu_failing_probe_[iProbe];
          truth_tag    = truth_ele_tag_[iTag];
          truth_probe  = truth_mu_failing_probe_[iProbe];
          tagPid = -11*q_ele_tag_[iTag];
          probePid = -13*q_mu_failing_probe_[iProbe];
          TLorentzVector systemP4 = (*p4_tag) + (*p4_probe);
          mass = systemP4.M();
          if (verbose) printf("\t\tmade a FAILING e-mu pair! pTs %f, %f; system mass %f, total charge %d e\n", p4_tag->Pt(), p4_probe->Pt(), mass, qtag+qprobe);
          electron_pair_tree->Fill();
        } 

      }
    }
    if(do_muons) { // muon tag electron probe pairs
      for(unsigned iTag=0; iTag < p4_mu_tag_.size(); iTag++) {
        for (unsigned iProbe=0; iProbe < p4_ele_passing_probe_.size(); iProbe++) {
          pass = 1;
          *p4_tag     = p4_mu_tag_[iTag];
          *p4_probe   = p4_ele_passing_probe_[iProbe];
          *genp4_tag     = genp4_mu_tag_[iTag];
          *genp4_probe   = genp4_ele_passing_probe_[iProbe];
          qtag    = q_mu_tag_[iTag];
          qprobe  = q_ele_passing_probe_[iProbe];
          tagPid = -13*q_mu_tag_[iTag];
          probePid = -11*q_ele_passing_probe_[iProbe];
          truth_tag    = truth_mu_tag_[iTag];
          truth_probe  = truth_ele_passing_probe_[iProbe];
          TLorentzVector systemP4 = (*p4_tag) + (*p4_probe);
          mass = systemP4.M();
          if(mass<40.||mass>140.) continue;
          if (verbose)
            printf("\t\tmade a PASSING mu-e pair! pTs %f, %f; system mass %f, total charge %d e\n", p4_tag->Pt(), p4_probe->Pt(), mass, qtag+qprobe);
          muon_pair_tree->Fill();
        }
        for(unsigned iProbe=0; iProbe < p4_ele_failing_probe_.size(); iProbe++) {
          pass = 0;
          *p4_tag     = p4_mu_tag_[iTag];  
          *p4_probe   = p4_ele_failing_probe_[iProbe];
          qtag    = q_mu_tag_[iTag];
          qprobe  = q_ele_failing_probe_[iProbe];
          truth_tag    = truth_mu_tag_[iTag];
          truth_probe  = truth_ele_failing_probe_[iProbe];
          tagPid = -13*q_mu_tag_[iTag];
          probePid = -11*q_ele_failing_probe_[iProbe];
          TLorentzVector systemP4 = (*p4_tag) + (*p4_probe);
          mass = systemP4.M();
          if(mass<40.||mass>140.) continue;
          if (verbose)
            printf("\t\tmade a FAILING mu-e pair! pTs %f, %f; system mass %f, total charge %d e\n", p4_tag->Pt(), p4_probe->Pt(), mass, qtag+qprobe);
          muon_pair_tree->Fill();
        } 

      }
    }
    
    //////////////////////////
    // Make Lepton-CH pairs //
    //////////////////////////
    unsigned nZCand=nPassingMuPairs+nFailingMuPairs+nPassingElePairs+nFailingElePairs;
    
    // First, assemble the collection of CH probes
    vector<panda::PFCand*> chProbes;
    vector<bool> chEleProbePass, chEleProbeReco, chMuProbePass, chMuProbeReco;
    vector<bool> chFromTheHardestPV;
    vector<float> chProbeTruths;
    vector<TLorentzVector> chProbeGenP4s;
    UShort_t iPF=0; unsigned pfRangeMax=2084;
    chEleProbePass.reserve(event.pfCandidates.size());
    chEleProbeReco.reserve(event.pfCandidates.size());
    chMuProbePass.reserve(event.pfCandidates.size());
    chMuProbeReco.reserve(event.pfCandidates.size());
    chFromTheHardestPV.reserve(event.pfCandidates.size());
    chProbeTruths.reserve(event.pfCandidates.size());
    chProbeGenP4s.reserve(event.pfCandidates.size());
    for(auto& the_vertex : event.vertices) if(the_vertex.pfRangeMax<pfRangeMax && the_vertex.pfRangeMax>0) pfRangeMax=the_vertex.pfRangeMax;
    for (auto& pfCand : event.pfCandidates) {
      bool fromTheHardestPV = (iPF < pfRangeMax);
      //if(iPF >= pfRangeMax) break; // only pf candidates from the primary vertex
      if(pfCand.pt()<3.) continue; 
      if(pfCand.eta()<-2.5 || pfCand.eta()>2.5) continue;

      // pf candidate species
      bool isChargedTrack=true;
      switch(pfCand.ptype) {
        case panda::PFCand::hp:
        case panda::PFCand::ep:
        case panda::PFCand::mup:
          probePid=211;
          qprobe=1;
          break;
        case panda::PFCand::hm:
        case panda::PFCand::em:
        case panda::PFCand::mum:
          probePid=-211;
          qprobe=-1;
          break;
        default:
          isChargedTrack=false;
          break;
      }
      if(!isChargedTrack) continue;
      if (!pfCand.track.isValid()) {
        if(verbose) printf("\tSkipping charged PF candidate because it does not have a track reference!\n");
        continue;
      }

      bool isChargedHadron = (pfCand.ptype == panda::PFCand::hp || pfCand.ptype == panda::PFCand::hm);
      *p4_probe=pfCand.p4();

      // Try to match the charged track to a reconstructed electron/muon
      // If it matches, see if the lepton passes the probe test selection.
      bool chMatchedToRecoMuon=false, chMatchedToRecoElectron=false;
      bool chRecoMuonIsPassing=false, chRecoElectronIsPassing=false;
      if (event.muons.size()>0) for (unsigned iM = 0; iM != event.muons.size(); ++iM){
        //if(p4_probe->DeltaR(((panda::Muon)event.muons[iM]).p4())<0.3) {chMatchedToRecoMuon=true; break;}
        if(((panda::Muon)event.muons[iM]).matchedPF.isValid() && ((panda::Muon)event.muons[iM]).matchedPF.get() == &pfCand) {
          chMatchedToRecoMuon=true;
          vector<int> passing_ids = passProbe(event.muons[iM], gettagprobe, false);
          for (unsigned long i = 0; i < passing_ids.size(); i++) if(passing_ids[i]==2) chRecoMuonIsPassing=true;
        }
      } 
      if (event.electrons.size()>0) for (unsigned iE = 0; iE != event.electrons.size(); ++iE){
        //if(p4_probe->DeltaR(((panda::Electron)event.electrons[iE]).p4())<0.3) {chMatchedToRecoElectron=true; break;}
        if(((panda::Electron)event.electrons[iE]).matchedPF.isValid() && ((panda::Electron)event.electrons[iE]).matchedPF.get() == &pfCand) {
          chMatchedToRecoElectron=true;
          vector<int> passing_ids = passProbe(event.electrons[iE], gettagprobe, false);
          for (unsigned long i = 0; i < passing_ids.size(); i++) if(passing_ids[i]==2) chRecoElectronIsPassing=true;
        }
      } 
      
      bool chMatchedToGenMuon=false, chMatchedToGenElectron=false;
      float chTruth=0, dRCHLepton; 
      TLorentzVector genP4;
      if (!real_data && truth_matching) {
        for (unsigned iG = 0; iG != event.genParticles.size(); ++iG){
          auto& genParticle = event.genParticles[iG];
          if (genParticle.finalState != 1) continue;
          TVector3 genParticle3; genParticle3.SetPtEtaPhi(genParticle.pt(),genParticle.eta(),genParticle.phi());
          dRCHLepton = genParticle3.DeltaR(pfCand.p4().Vect());
          if (TMath::Abs(genParticle.pdgid) == 13 && dRCHLepton < truth_matching_dR) {
            chTruth = TMath::Max(0.001f, dRCHLepton);
            genP4.SetPtEtaPhiM(genParticle.pt(),genParticle.eta(),genParticle.phi(),0.106);
            chMatchedToGenMuon=true;
            break;
          } else if (TMath::Abs(genParticle.pdgid) == 11 && dRCHLepton < truth_matching_dR) {
            chTruth = TMath::Max(0.001f, dRCHLepton);
            genP4.SetPtEtaPhiM(genParticle.pt(),genParticle.eta(),genParticle.phi(),511e-6);
            chMatchedToGenElectron=true;
            break;
          }
        }
        if(chTruth==0) continue;
      }

      bool duplicateTrack=false; panda::PFCand* theDupe;
      // Debugging: Check if two charged PF candidates have the same track
      if(chProbes.size()>0) for(unsigned iCH=0; iCH<chProbes.size(); iCH++) {
        if(pfCand.track.get() == chProbes[iCH]->track.get()) { duplicateTrack=true; theDupe=chProbes[iCH]; break; }
      }
      if(duplicateTrack) {
        if(verbose) printf("\tWarning: duplicate track between this charged PF (pT=%.1f, eta=%.2f, type=%d) and existing probe (pT=%.1f, eta=%.2f, type=%d)\n", pfCand.pt(), pfCand.eta(), pfCand.ptype, theDupe->pt(),theDupe->eta(), theDupe->ptype);
        continue;
      }

      chProbes.push_back(&pfCand);
      chProbeGenP4s.push_back(genP4);
      chEleProbeReco.push_back(chMatchedToRecoElectron);
      chMuProbeReco.push_back(chMatchedToRecoMuon);
      chEleProbePass.push_back(chMatchedToRecoElectron && chRecoElectronIsPassing);
      chMuProbePass.push_back(chMatchedToRecoMuon && chRecoMuonIsPassing);
      chProbeTruths.push_back(chTruth);
      chFromTheHardestPV.push_back(fromTheHardestPV);
      iPF++;
    } // Done assembling the collection of CH probes
    
    // Make same sign L-Pi pairs; also perform arbitration for the probes
    // Muon tags
    if(do_muons) for (unsigned iTag=0; iTag < p4_mu_tag_.size(); iTag++) {
      *p4_tag     = p4_mu_tag_[iTag];
      *genp4_tag  = genp4_mu_tag_[iTag];
      qtag    = q_mu_tag_[iTag];
      tagPid = -13*q_mu_tag_[iTag];
      truth_tag    = truth_mu_tag_[iTag];
      float bestProbeMass, bestProbeMassDiff=999; unsigned bestProbeIdx=0; bool foundBestProbe=false;
      vector<float> chPairMass; chPairMass.reserve(chProbes.size());
      vector<float> chPairDVz ; chPairDVz .reserve(chProbes.size());
      
      // We have to loop over the probes completely twice
      // This first loop computes the pair masses and the probe multiplicity
      if(verbose) printf("\tbeginning Muon tag probe multiplicity calculation\n");
      probeMultiplicity=0; 
      for(unsigned iCH=0; iCH < chProbes.size(); iCH++) { 
        *p4_probe=chProbes[iCH]->p4(); 
        *genp4_probe = chProbeGenP4s[iCH];
        TLorentzVector systemP4 = (*p4_tag) + (*p4_probe);
        mass = systemP4.M();
        chPairMass.push_back(mass);
        qprobe=(chProbes[iCH]->ptype==panda::PFCand::hp || chProbes[iCH]->ptype==panda::PFCand::ep || chProbes[iCH]->ptype==panda::PFCand::mup)? 1:-1;
        float vz = 99;
        if(chProbes[iCH]->vertex.isValid()) vz = chProbes[iCH]->vertex.get()->z;
        if(verbose) printf("\tTag (pT=%.1f, q=%d, vz=%.2f) probe (pT=%.1f, q=%d, vz=%.2f) mass=%.2f\n", p4_tag->Pt(), qtag, vz_mu_tag_[iTag], p4_probe->Pt(), qprobe, vz, mass);
        if((qtag+qprobe)!=0) continue;
        //if( p4_tag->DeltaR(*p4_probe) < .05) continue;
        if(mass<60.) continue;
        float dvz=fabs(vz_mu_tag_[iTag] - vz);
        chPairDVz.push_back(dvz);
        if(dvz > 4) continue;
        probeMultiplicity++;
        if(verbose) printf("\tincrementing probeMultiplicity, new value -> %d\n", probeMultiplicity);
        // Check if this is the pair with the mass closest to the Z mass
        float massDiff = fabs(mass-massForArbitration);
        if(!foundBestProbe || massDiff<bestProbeMassDiff) { 
          bestProbeIdx = iCH;
          bestProbeMass=mass;
          bestProbeMassDiff=massDiff;
          foundBestProbe=true;
        }
      }
      
      // Now that we know the probe multiplicity, store the Mu-Pi Pairs
      for(unsigned iCH=0; iCH < chProbes.size(); iCH++) { 
        mass=chPairMass[iCH];
        if(chPairDVz[iCH] > 4) continue;
        if(mass<40.||mass>140.) continue;
        *p4_probe=chProbes[iCH]->p4(); 
        *genp4_probe = chProbeGenP4s[iCH];
        qprobe=(chProbes[iCH]->ptype==panda::PFCand::hp || chProbes[iCH]->ptype==panda::PFCand::ep || chProbes[iCH]->ptype==panda::PFCand::mup)? 1:-1;
        probePid=211*qprobe;
        truth_probe  = chProbeTruths[iCH];
        pass=0;
        if((qtag+qprobe)!=0) { // Handle Same-Sign Mu-Pi Pairs
          if(nZCand!=0 || chMuProbeReco[iCH] || chEleProbeReco[iCH]) continue; // Z veto
          if (verbose) printf("\t\tmade a Same Sign Mu Track pair! pTs %f, %f; system mass %f, total charge %d e, pass=%d\n", p4_tag->Pt(), p4_probe->Pt(), mass, qtag+qprobe, pass);
          muon_pair_tree->Fill();
        } else { // Handle Opposite Sign Mu-Pi Pairs
          bestMass = (bestProbeIdx==iCH);
          if(chEleProbePass[iCH]) pass |= 1<<11;
          if(chMuProbePass[iCH])  pass |= 1<<13;
          if (verbose) printf("\t\tmade an Opposite Sign Mu Track pair! pTs %f, %f; system mass %f, total charge %d e\n", p4_tag->Pt(), p4_probe->Pt(), mass, qtag+qprobe);
          muon_pair_tree->Fill();
        }
      }
      // Now store all opposite sign Mu-Pi pairs, but take note of the probe multiplicity per tag and mark the pair that has the "best mass"
      //if(!foundBestMuProbe) continue;
      //panda::PFCand *bestProbe=chProbes[bestMuProbeIdx];
      //*p4_probe=bestProbe->p4();
      //*genp4_probe = chProbeGenP4s[bestMuProbeIdx];
      //truth_probe = chProbeTruths[bestMuProbeIdx];
      //pass=chMuProbePass[bestMuProbeIdx];
      //mass=bestMuProbeMass;
      //qprobe=(bestProbe->ptype==panda::PFCand::hp || bestProbe->ptype==panda::PFCand::ep || bestProbe->ptype==panda::PFCand::mup)? 1:-1;
      //probePid=211*qprobe;
      //if(qprobe+qtag==0 && ((tagPid==13&&probePid==211)||(tagPid==-13&&probePid==-211))) {
      //  if (verbose) printf("\t\tmade a Mu-Track pair! pTs %f, %f; system mass %f, total charge %d e\n", p4_tag->Pt(), p4_probe->Pt(), mass, qtag+qprobe);
      //  muon_pair_tree->Fill();
      //}
    }
    //Electron tags
    if(do_electrons) for (unsigned iTag=0; iTag < p4_ele_tag_.size(); iTag++) {
      *p4_tag     = p4_ele_tag_[iTag];
      *genp4_tag  = genp4_ele_tag_[iTag];
      qtag    = q_ele_tag_[iTag];
      tagPid = -13*q_ele_tag_[iTag];
      truth_tag    = truth_ele_tag_[iTag];
      
      float bestProbeMass, bestProbeMassDiff=999; unsigned bestProbeIdx=0; bool foundBestProbe=false;
      vector<float> chPairMass; chPairMass.reserve(chProbes.size());
      vector<float> chPairDVz ; chPairDVz .reserve(chProbes.size());
      
      // We have to loop over the probes completely twice
      // This first loop computes the pair masses and the probe multiplicity
      probeMultiplicity=0; 
      for(unsigned iCH=0; iCH < chProbes.size(); iCH++) { 
        *p4_probe=chProbes[iCH]->p4(); 
        *genp4_probe = chProbeGenP4s[iCH];
        TLorentzVector systemP4 = (*p4_tag) + (*p4_probe);
        mass = systemP4.M();
        chPairMass.push_back(mass);
        qprobe=(chProbes[iCH]->ptype==panda::PFCand::hp || chProbes[iCH]->ptype==panda::PFCand::ep || chProbes[iCH]->ptype==panda::PFCand::mup)? 1:-1;
        if((qtag+qprobe)!=0) continue;
        //if( p4_tag->DeltaR(*p4_probe) < .05) continue;
        if(mass<60.) continue;
        float vz = 99;
        if(chProbes[iCH]->vertex.isValid()) vz = chProbes[iCH]->vertex.get()->z;
        float dvz=fabs(vz_ele_tag_[iTag] - vz);
        chPairDVz.push_back(dvz);
        if(dvz > 4) continue;
        probeMultiplicity++;
        // Check if this is the pair with the mass closest to the Z mass
        float massDiff = fabs(mass-massForArbitration);
        if(!foundBestProbe || massDiff<bestProbeMassDiff) { 
          bestProbeIdx = iCH;
          bestProbeMass=mass;
          bestProbeMassDiff=massDiff;
          foundBestProbe=true;
        }
      }
      // Now that we know the probe multiplicity, store the Mu-Pi Pairs
      for(unsigned iCH=0; iCH < chProbes.size(); iCH++) { 
        mass=chPairMass[iCH];
        if(mass<40.||mass>140.) continue;
        if(chPairDVz[iCH] > 4) continue;
        *p4_probe=chProbes[iCH]->p4(); 
        *genp4_probe = chProbeGenP4s[iCH];
        qprobe=(chProbes[iCH]->ptype==panda::PFCand::hp || chProbes[iCH]->ptype==panda::PFCand::ep || chProbes[iCH]->ptype==panda::PFCand::mup)? 1:-1;
        probePid=211*qprobe;
        truth_probe  = chProbeTruths[iCH];
        pass=0;
        if((qtag+qprobe)!=0) { // Handle Same-Sign E-Pi Pairs
          if(nZCand!=0 || chMuProbeReco[iCH] || chEleProbeReco[iCH]) continue; // Z veto
          if (verbose) printf("\t\tmade a Same Sign E Track pair! pTs %f, %f; system mass %f, total charge %d e, pass=%d\n", p4_tag->Pt(), p4_probe->Pt(), mass, qtag+qprobe, pass);
          electron_pair_tree->Fill();
        } else { // Handle Opposite Sign Mu-Pi Pairs
          bestMass = (bestProbeIdx==iCH);
          if(chEleProbePass[iCH]) pass |= 1<<11;
          if(chMuProbePass[iCH])  pass |= 1<<13;
          if (verbose) printf("\t\tmade an Opposte Sign E Track pair! pTs %f, %f; system mass %f, total charge %d e\n", p4_tag->Pt(), p4_probe->Pt(), mass, qtag+qprobe);
          electron_pair_tree->Fill();
        }
      }
      
    }
  }
  input_file->Close();
  delete p4_tag; delete p4_probe;
  delete genp4_tag; delete genp4_probe;

  //Histogram containing total number of entries
  TH1F nEntriesHist("hDTotalMCWeight","hDTotalMCWeight", 1,0,2);
  nEntriesHist.SetBinContent(1, nEntries);
  TH1D sum_MCW("sum_weights","sum_weights", 1, 0, 2);
  sum_MCW.SetBinContent(1,sum_mc_weights);
  //  TH1F Xsec("Xsec","Xsec", 1, 0 ,2);
  //Xsec.SetBinContent(1,gettagprobe.xsec);
  // TH1F*sum_MCW = new TH1F("sum_weights","sum_weights",1,0,10000000000);
  //sum_MCW->Fill(sum_mc_weights);
  //sum_MCW->Draw();

  //save tnp trees
  if(do_electrons) {
    electron_outfile->cd();
    electron_pair_tree->Write(electron_pair_tree->GetName(), TObject::kOverwrite);
    electron_outfile->WriteTObject(&nEntriesHist);
    electron_outfile->WriteTObject(&sum_MCW);
    //electron_outfile->WriteTObject(&Xsec);
    printf("%lld electron pair events!\n", electron_pair_tree->GetEntries());
    electron_outfile->Close();
  }
  if(do_muons) {
    muon_outfile->cd();
    muon_pair_tree->Write(muon_pair_tree->GetName(), TObject::kOverwrite);
    muon_outfile->WriteTObject(&nEntriesHist);
    muon_outfile->WriteTObject(&sum_MCW);
    //muon_outfile->WriteTObject(&Xsec);
    printf("%lld muon pair events!\n", muon_pair_tree->GetEntries());
    muon_outfile->Close();
  }
  
  //Summary of completed events:
  printf("%lld\n", sum_mc_weights);
  printf("Complete. %lld events processed, total MC events %lld \n\n\n", nEntries, sum_mc_weights);
  //  printf("run num [%d,%d] \n", min_runNum, max_runNum); 

  //delete puWeights;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////


//ADDITIONAL FUNCTIONS:

double selectIsoCut(int type, int pdgId, double eta) {
  ///selectIsoCut is used in PROBE section in conjuction with selector
    bool isEB = TMath::Abs(eta) < 1.479;
  if     (TMath::Abs(pdgId) == 15) {
    if (type==0) return 10000;
    return 4.5;
  }
  if     (TMath::Abs(pdgId) == 13) {
    if (type==0) return 10000;
    if (type==1) return .25;
    if (type==4) return .25;
    if (type==6) return .06;
    return 0.15;
  }
  else if(TMath::Abs(pdgId) == 11) {
    //if     (type == "veto")   return (isEB ? 0.1260 : 0.1440);
    //else if(type == "loose")  return (isEB ? 0.0893 : 0.1210);
    //else if(type == "medium") return (isEB ? 0.0766 : 0.0678);
    //else if(type == "tight")  return (isEB ? 0.0354 : 0.0646);
    if     (type == 0)   return 10000;
    else if(type == 1)   return (isEB ? 0.1260 : 0.1440);
    else if(type == 4)   return (isEB ? 0.0893 : 0.1210);
    else if(type == 5)   return (isEB ? 0.0766 : 0.0678);
    else if(type == 6)   return (isEB ? 0.0354 : 0.0646);
    else if(type == 7)   return 0.06;
  }
  //printf("Problem with selectIsoCut! type=%d, pdgId=%d, eta=%f\n", type, pdgId,eta);
  //assert(0);

  return 0.0;
}

bool selector(
              //Uses the standards assigned by CMS to qualify leptons(muons)
              panda::Muon const& muon,
              int id_bit,
              int iso_bit
              ) {
  double iso = muon.combIso() / muon.pt();

  switch (id_bit) {
  case 100:
    return muon.loose && iso < selectIsoCut(iso_bit, 13, muon.eta());
  case 101:
    return (muon.medium || muon.mediumBtoF) && /*       (sel_bits & (0x1 << 9)) != 0 && // Tight IP Cut*/
      iso < selectIsoCut(iso_bit, 13, muon.eta());
  case 102:
    return muon.tight && /*       (sel_bits & (0x1 << 9)) != 0 && // Tight IP Cut*/
      iso < selectIsoCut(iso_bit, 13, muon.eta());
  case 103:
    return (muon.medium || muon.mediumBtoF) && /*       (sel_bits & (0x1 << 8)) != 0 && // Medium IP Cut*/
      iso < selectIsoCut(iso_bit, 13, muon.eta());
  case 4:
    return muon.loose && iso < selectIsoCut(iso_bit, 13, muon.eta());
  case 5:
    return (muon.medium || muon.mediumBtoF) && iso < selectIsoCut(iso_bit, 13, muon.eta());
  case 6:
    return muon.tight && iso < selectIsoCut(iso_bit, 13, muon.eta());
  case 7:
    return muon.tight && iso < selectIsoCut(iso_bit, 13, muon.eta());
  case 0:
    return muon.loose;
  default:
    return true;

  }
}


bool selector(
              //Uses the standards assigned by CMS to qualify leptons(electrons)
              panda::Electron const& electron,
              int id_bit,
              int iso_bit
              ) {
  double iso = electron.combIso() / electron.pt();
  float aeta = fabs(electron.eta());
  float pt = electron.pt();
  switch (id_bit) {
  case 200:
    return electron.tight &&
      //      (sel_bits & (0x1 << 15)) != 0 && // Triple charge requirement
      //      (sel_bits & (0x1 << 16)) != 0 && // No Missing Hits
      iso < selectIsoCut(iso_bit, 11, electron.eta());
  case 4:
    return electron.loose && iso < selectIsoCut(iso_bit, 11, electron.eta());
  case 5:
    return electron.medium && iso < selectIsoCut(iso_bit, 11, electron.eta());
  case 6:
    return electron.tight && iso < selectIsoCut(iso_bit, 11, electron.eta());
  case 7:
    return (electron.mvaWP80 && 
      ((
        aeta < 1.4442 && 
        electron.sieie < 0.012 && 
        electron.hOverE < 0.09 &&
        electron.ecalIso < 0.4*pt && 
        electron.hcalIso < 0.25*pt && 
        electron.trackIso < 0.18*pt &&
        fabs(electron.dEtaInSeed) < 0.0095 && 
        fabs(electron.dPhiIn) < 0.065
      ) || (
        aeta > 1.5660 && 
        electron.sieie < 0.033 && 
        electron.hOverE < 0.09 &&
        electron.ecalIso < 0.45*pt && 
        electron.hcalIso < 0.28*pt &&
        electron.trackIso < 0.18*pt
    )) && iso < selectIsoCut(iso_bit, 11, electron.eta()));
  default:
    return true;
  }
}



void InitializeJetIdCuts(Float_t fMVACut[4][4])
//Creates 4x4 matrix of cut values
{ 
  float cutValues[4][4] = {
    -0.95, -0.96 ,-0.94, -0.95,
    -0.95, -0.96 ,-0.94, -0.95,
    -0.15, -0.26 ,-0.16, -0.16,
    -0.15, -0.26 ,-0.16, -0.16
  };
  
  for(int i=0; i<4; i++){
    for(int j=0; j<4; j++){
      fMVACut[i][j] = cutValues[i][j];
    }
  }
}


bool passJetId(panda::Jet const& jet, Float_t fMVACut[4][4]){
  //Determines the validity of a given jet
    int lPtId = 3;
  if     (jet.pt() < 10.)
    lPtId = 0;    
  else if(jet.pt() < 20.)
    lPtId = 1;
  else if(jet.pt() < 30.)
    lPtId = 2;      
  
  int lEtaId = 3;                                           
  if     (jet.eta() < 2.50)
    lEtaId = 0;
  else if(jet.eta() < 2.75)
    lEtaId = 1;            
  else if(jet.eta() < 3.00)
    lEtaId = 2;
  
  if (jet.puid > fMVACut[lPtId][lEtaId])
    return true;
  
  return false;                         
  
}



vector<int> passProbe(panda::Electron const& electron, leptonVariables gettagprobe, bool pass_tag_trigger){
  vector<int> list_of_ids;
  //Determines whether a lepton(electron) passes the probe
  int i = 0;   
  if ((gettagprobe.passing_probe_iso >=0 && (selector(electron, gettagprobe.probe_id, gettagprobe.probe_iso) && !selector(electron, gettagprobe.passing_probe_id, gettagprobe.passing_probe_iso)))||(gettagprobe.passing_probe_iso < 0 && (selector(electron, gettagprobe.probe_id, gettagprobe.probe_iso)&& !electron.triggerMatch[gettagprobe.passing_probe_id]))){
     i = 1;
     list_of_ids.push_back(i);
   }
   if ((gettagprobe.passing_probe_iso >=0 && (selector(electron, gettagprobe.probe_id, gettagprobe.probe_iso)&& selector(electron, gettagprobe.passing_probe_id, gettagprobe.passing_probe_iso)))||(gettagprobe.passing_probe_iso < 0 && (selector(electron, gettagprobe.probe_id, gettagprobe.probe_iso)&& electron.triggerMatch[gettagprobe.passing_probe_id]))){
     i = 2;
     list_of_ids.push_back(i);
   }
   if ((electron.pt()>=gettagprobe.tag_pt_min && TMath::Abs(electron.eta()) <=gettagprobe.tag_eta_max && selector(electron, gettagprobe.tag_id, gettagprobe.tag_iso) && pass_tag_trigger)){
     i = 3;
     list_of_ids.push_back(i);
   }
   return list_of_ids;
}




vector<int> passProbe(panda::Muon const& muon,leptonVariables gettagprobe, bool pass_tag_trigger){
  //Determines whether a lepton(muon) passes the probe
  vector<int> list_of_ids;
  int i = 0;
  if ((gettagprobe.passing_probe_iso >=0 && (selector(muon, gettagprobe.probe_id, gettagprobe.probe_iso)&& !selector(muon, gettagprobe.passing_probe_id, gettagprobe.passing_probe_iso)))||(gettagprobe.passing_probe_iso < 0 && (selector(muon, gettagprobe.probe_id, gettagprobe.probe_iso)&& !muon.triggerMatch[gettagprobe.passing_probe_id]))){
    i = 1;    
    list_of_ids.push_back(i);
  }
  if ((gettagprobe.passing_probe_iso >=0 && (selector(muon, gettagprobe.probe_id, gettagprobe.probe_iso)&& selector(muon, gettagprobe.passing_probe_id, gettagprobe.passing_probe_iso)))||(gettagprobe.passing_probe_iso < 0 && (selector(muon, gettagprobe.probe_id, gettagprobe.probe_iso)&& muon.triggerMatch[gettagprobe.passing_probe_id]))){
    i = 2;
    list_of_ids.push_back(i);
  }
  if ((muon.pt()>=gettagprobe.tag_pt_min && TMath::Abs(muon.eta()) <=gettagprobe.tag_eta_max && selector(muon, gettagprobe.tag_id, gettagprobe.tag_iso) && pass_tag_trigger)){
    i = 3;
    list_of_ids.push_back(i);
  }
  return list_of_ids;
}






void animator(int iEntry,int nEntries){
  if (iEntry %50 == 0){
    if (iEntry <= .25 * nEntries){
      printf("|<~~> * * * * * * * *|\r");
    }
    else if (iEntry <= .5 * nEntries){
      printf("|      <~~> * * * * *|\r");
    }
    else if (iEntry <= .75 * nEntries){
      printf("|            <~~> * *|\r");
    }
    else if (iEntry <= .99 * nEntries){
      printf("|                <~~>|\r");
    }
  }
}
  
