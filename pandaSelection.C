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

typedef std::map<UInt_t,std::vector<std::pair <UInt_t, UInt_t> > > MapType;
//string jsonFile = "/home/dhsu/CMSSW_8_0_26_patch1/src/MitVBFAnalysis/data/Cert_294927-300575_13TeV_PromptReco_Collisions17_JSON.txt";
string jsonFile = "/home/dhsu/CMSSW_8_0_26_patch1/src/MitVBFAnalysis/data/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt";

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
  
  //OUTPUT FILE:

  //declare output variables
  unsigned int out_runNum, // event ID
    out_lumiSec,
    out_evtNum,
    out_npv, // number of primary vertices
    pass; // whether probe passes requirements
  float        npu=1;                     // mean number of expected pileup
  float        scale1fb=1;                  // event weight per 1/fb
  float        mass;                      // tag-probe mass
  int          qtag, qprobe;              // tag, probe charge
  float        truth_tag, truth_probe;              // tag, probe truth
  float        met;                             // missing ET
  int          njets;                           // number of jets
  int          tagPid,probePid; // particle ID
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
    electron_pair_tree->Branch("qtag",     &qtag,     "qtag/I"     );  
    electron_pair_tree->Branch("qprobe",   &qprobe,   "qprobe/I"   );  
    electron_pair_tree->Branch("njets",    &njets,    "njets/I"   );  
    electron_pair_tree->Branch("met",      &met,      "met/F"   );  
    electron_pair_tree->Branch("tag",   "TLorentzVector", &p4_tag   );  
    electron_pair_tree->Branch("probe", "TLorentzVector", &p4_probe );          
    electron_pair_tree->Branch("tagPid",      &tagPid,      "tagPid/I"   );  
    electron_pair_tree->Branch("probePid",      &probePid,      "probePid/I"   );  
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
    muon_pair_tree->Branch("qtag",     &qtag,     "qtag/I"     );  
    muon_pair_tree->Branch("qprobe",   &qprobe,   "qprobe/I"   );  
    muon_pair_tree->Branch("tag",   "TLorentzVector", &p4_tag   );  
    muon_pair_tree->Branch("probe", "TLorentzVector", &p4_probe );      
    muon_pair_tree->Branch("njets",    &njets,    "njets/I"   );  
    muon_pair_tree->Branch("met",      &met,      "met/F"   );  
    muon_pair_tree->Branch("tagPid",      &tagPid,      "tagPid/I"   );  
    muon_pair_tree->Branch("probePid",      &probePid,      "probePid/I"   );  
    if(!real_data) {
      muon_pair_tree->Branch("genTag",   "TLorentzVector", &genp4_tag   );  
      muon_pair_tree->Branch("genProbe", "TLorentzVector", &genp4_probe );      
      muon_pair_tree->Branch("truth_tag",     &truth_tag,     "truth_tag/F"     );  
      muon_pair_tree->Branch("truth_probe",   &truth_probe,   "truth_probe/F"   );  
    }
  }


  //////////////////////////////////////////////////////////////////////////////////////////
  
  //Load pileup corrections
  TFile *puFile = TFile::Open("/home/dhsu/CMSSW_8_0_26_patch1/src/Clyde_UROP/puWeights_80x_37ifb.root", "READ");
  TH1D *puWeights = (TH1D*)puFile->Get("puWeights"); puWeights->SetDirectory(0);
  puFile->Close();  

  //INPUT FILE:
  //Setting up Tree files and various input variables  
  printf("Trying to open file %s...\n", input_file_name.c_str());
  TFile*  input_file=TFile::Open(input_file_name.c_str(),"READ");
  if(!input_file || !input_file->IsOpen()) {
    printf("\"I couldn't open it\"~Alex Jones\n");
    assert(0); return;
  }
  TTree* tree = (TTree*)input_file->Get("events"); // get the tree object from the file
  panda::Event event; // create an Event object
  event.setStatus(*tree, {"!*"});
  event.setAddress(*tree, {"runNumber", "lumiNumber", "eventNumber", "muons", "electrons", "npv", "npvTrue", "genParticles","isData","pfMet","chsAK4Jets","pfCandidates","vertices","weight"}); 
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
  RoccoR rochesterCorrection("/home/dhsu/CMSSW_8_0_26_patch1/src/Clyde_UROP/rcdata.2016.v3");
  
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
    std::pair<unsigned int, unsigned int> runLumi(event.runNumber, event.lumiNumber);      
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
    unsigned int goodIsTight = 0;
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
      
      //DEBUG: Jet Pt printout:
      if (verbose)
        std::cout << jet.pt() << std::endl;

      //Determine if Jet is actually electron: if not, add to list of Jets
      bool isElectron = false; bool isMuon = false;
      for(unsigned int iE = 0; iE != event.electrons.size(); ++iE){
        auto& electron  = event.electrons[iE];
        if (electron.loose){
          if(jet.dR(electron) < 0.16) { isElectron = true; break; }
        }
      }
      if (isElectron == true) continue;
        
      //Determine if Jet is actually muon: if not, add to list of Jets
      for(unsigned int iM = 0; iM != event.muons.size(); ++iM){
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

    // Loop over Electrons
    if (event.electrons.size()>0){
      for (unsigned iE = 0; iE != event.electrons.size(); ++iE){
        auto& electron = event.electrons[iE];
        if(electron.smearedPt < 10.) continue;
        bool electron_trigger_matched = false;
        double truth = 0;
        int charge = electron.charge;
        bool pass_tag_trigger = false;
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
          pass_tag_trigger = (electron.triggerMatch[panda::Electron::fEl27Tight]);
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
            break;
          default:
            break;
          }
        }
      }        
    }
    if (verbose){
      printf("Number of vectors in electron tag list:%lu\n", p4_ele_tag_.size());
      printf("Number of vectors in electron pass probe list:%lu\n", p4_ele_passing_probe_.size());
    }

    // Loop over Muons
    if (event.muons.size()>0){
      for (unsigned iM = 0; iM != event.muons.size(); ++iM){
        auto& muon = event.muons[iM];
        bool muon_trigger_matched = false;
        int charge = muon.charge;
        double truth;
        bool pass_tag_trigger = false;
        double ptCorrection=1;
        TLorentzVector genP4;
        if (!real_data && truth_matching){
          for (unsigned iG = 0; iG != event.genParticles.size(); ++iG){
            auto& genParticle = event.genParticles[iG];
            double dR;
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
            break;
          }
        }
      }
    }
    if (verbose){
      printf("Number of vectors in muon tag list:%lu\n", p4_mu_tag_.size());
      printf("Number of vectors in muon pass probe list:%lu\n", p4_mu_passing_probe_.size());
    }
    
    //ELECTRON PAIR ASSOCIATION:
    //associating electron pairs and filling tree
    if(do_electrons) {
      for(unsigned int iTag=0; iTag < p4_ele_tag_.size(); iTag++) {
        // passing probes for electrons
        for(unsigned int iProbe=0; iProbe < p4_ele_passing_probe_.size(); iProbe++) {
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
          if(!real_data) scale1fb = event.weight*puWeights->GetBinContent(event.npvTrue);
          // make sure they aren't the same particle with tiny delta R veto
          if(!( p4_tag->DeltaR(*p4_probe) < .0001) ) {
            TLorentzVector systemP4 = (*p4_tag) + (*p4_probe);
            mass = systemP4.M();
            electron_pair_tree->Fill();
          }
        }  
        // failing probes for electrons and filling tree
        for(unsigned int iProbe=0; iProbe < p4_ele_failing_probe_.size(); iProbe++) {
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
          if (!real_data) scale1fb = event.weight*puWeights->GetBinContent(event.npvTrue);
          // make sure they aren't the same particle with tiny delta R veto
          if(!( p4_tag->DeltaR(*p4_probe) < .0001) ) {
            TLorentzVector systemP4 = (*p4_tag) + (*p4_probe);
            mass = systemP4.M();          
            electron_pair_tree->Fill();
          }         
        } 
      }
    }
    //MUON PAIR ASSOCIATION:
    // associating muon pairs and filling tree
    if (do_muons) {
      for (unsigned int iTag=0; iTag < p4_mu_tag_.size(); iTag++) {
        // passing probes for muons
        for (unsigned int iProbe=0; iProbe < p4_mu_passing_probe_.size(); iProbe++) {
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
          if (!real_data) scale1fb = event.weight*puWeights->GetBinContent(event.npvTrue);
          // make sure they aren't the same particle with tiny delta R veto
          if(!( p4_tag->DeltaR(*p4_probe) < .0001) ) {
            TLorentzVector systemP4 = (*p4_tag) + (*p4_probe);
            mass = systemP4.M();
            if (verbose)
              printf("\t\tmade a PASSING muon pair! pTs %f, %f; system mass %f, total charge %d e\n", p4_tag->Pt(), p4_probe->Pt(), mass, qtag+qprobe);
            muon_pair_tree->Fill();
          }
        }
        
        // failing probes for muons and filling tree
        for(unsigned int iProbe=0; iProbe < p4_mu_failing_probe_.size(); iProbe++) {
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
          if (!real_data)
            scale1fb = event.weight;
            //scale1fb  = (eMC_weights*gettagprobe.xsec*1000)/(sum_mc_weights);
          // make sure they aren't the same particle with tiny delta R veto
          if(!( p4_tag->DeltaR(*p4_probe) < .0001) ) {
            TLorentzVector systemP4 = (*p4_tag) + (*p4_probe);
            mass = systemP4.M();
            if (verbose)
              printf("\t\tmade a FAILING muon pair! pTs %f, %f; system mass %f, total charge %d e\n", p4_tag->Pt(), p4_probe->Pt(), mass, qtag+qprobe);
            muon_pair_tree->Fill();
          }
        } 
      }
    }
    if (verbose){
      printf("Electron tag/probe: %lu/%lu\n", q_ele_tag_.size(), q_ele_passing_probe_.size());    
      printf("Muon tag/probe: %lu/%lu\n", q_mu_tag_.size(), q_mu_passing_probe_.size());    
    }
    
    // Make e-mu pairs
    if(do_electrons && real_data) { // electron tag muon probe pairs
      for(unsigned int iTag=0; iTag < p4_ele_tag_.size(); iTag++) {
        // passing probes for muons
        for (unsigned int iProbe=0; iProbe < p4_mu_passing_probe_.size(); iProbe++) {
          pass = 1;
          *p4_tag     = p4_ele_tag_[iTag];
          *p4_probe   = p4_mu_passing_probe_[iProbe];
          qtag    = q_ele_tag_[iTag];
          qprobe  = q_mu_passing_probe_[iProbe];
          tagPid = -11*q_ele_tag_[iTag];
          probePid = -13*q_mu_passing_probe_[iProbe];
          truth_tag    = truth_ele_tag_[iTag];
          truth_probe  = truth_mu_passing_probe_[iProbe];
          if(!real_data) scale1fb = event.weight*puWeights->GetBinContent(event.npvTrue);
          TLorentzVector systemP4 = (*p4_tag) + (*p4_probe);
          mass = systemP4.M();
          if (verbose)
            printf("\t\tmade a PASSING e-mu pair! pTs %f, %f; system mass %f, total charge %d e\n", p4_tag->Pt(), p4_probe->Pt(), mass, qtag+qprobe);
          electron_pair_tree->Fill();
        }
        // failing probes for muons
        for(unsigned int iProbe=0; iProbe < p4_mu_failing_probe_.size(); iProbe++) {
          pass = 0;
          *p4_tag     = p4_ele_tag_[iTag];  
          *p4_probe   = p4_mu_failing_probe_[iProbe];
          qtag    = q_ele_tag_[iTag];
          qprobe  = q_mu_failing_probe_[iProbe];
          truth_tag    = truth_ele_tag_[iTag];
          truth_probe  = truth_mu_failing_probe_[iProbe];
          tagPid = -11*q_ele_tag_[iTag];
          probePid = -13*q_mu_failing_probe_[iProbe];
          if(!real_data) scale1fb = event.weight*puWeights->GetBinContent(event.npvTrue);
          TLorentzVector systemP4 = (*p4_tag) + (*p4_probe);
          mass = systemP4.M();
          if (verbose)
            printf("\t\tmade a FAILING e-mu pair! pTs %f, %f; system mass %f, total charge %d e\n", p4_tag->Pt(), p4_probe->Pt(), mass, qtag+qprobe);
          electron_pair_tree->Fill();
        } 

      }
    }
    if(do_muons && real_data) { // muon tag electron probe pairs
      for(unsigned int iTag=0; iTag < p4_mu_tag_.size(); iTag++) {
        for (unsigned int iProbe=0; iProbe < p4_ele_passing_probe_.size(); iProbe++) {
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
          if(!real_data) scale1fb = event.weight*puWeights->GetBinContent(event.npvTrue);
          TLorentzVector systemP4 = (*p4_tag) + (*p4_probe);
          mass = systemP4.M();
          if (verbose)
            printf("\t\tmade a PASSING mu-e pair! pTs %f, %f; system mass %f, total charge %d e\n", p4_tag->Pt(), p4_probe->Pt(), mass, qtag+qprobe);
          muon_pair_tree->Fill();
        }
        for(unsigned int iProbe=0; iProbe < p4_ele_failing_probe_.size(); iProbe++) {
          pass = 0;
          *p4_tag     = p4_mu_tag_[iTag];  
          *p4_probe   = p4_ele_failing_probe_[iProbe];
          qtag    = q_mu_tag_[iTag];
          qprobe  = q_ele_failing_probe_[iProbe];
          truth_tag    = truth_mu_tag_[iTag];
          truth_probe  = truth_ele_failing_probe_[iProbe];
          tagPid = -11*q_mu_tag_[iTag];
          probePid = -13*q_ele_failing_probe_[iProbe];
          if(!real_data) scale1fb = event.weight*puWeights->GetBinContent(event.npvTrue);
          TLorentzVector systemP4 = (*p4_tag) + (*p4_probe);
          mass = systemP4.M();
          if (verbose)
            printf("\t\tmade a FAILING mu-e pair! pTs %f, %f; system mass %f, total charge %d e\n", p4_tag->Pt(), p4_probe->Pt(), mass, qtag+qprobe);
          muon_pair_tree->Fill();
        } 

      }
    }
    // Make Lepton-CH pairs
    if(real_data){
      UShort_t iPF=0; int pfRangeMax=2084;
      for(auto& the_vertex : event.vertices)
        if(the_vertex.pfRangeMax<pfRangeMax && the_vertex.pfRangeMax>0) pfRangeMax=the_vertex.pfRangeMax;
      for (auto& pfCand : event.pfCandidates) {
        if(pfCand.ptype != panda::PFCand::hp && pfCand.ptype != panda::PFCand::hm) continue; // only pion/kaon/proton candidates
        if(iPF >= pfRangeMax) break; // only pf candidates from the primary vertex
        if(pfCand.pt()<10.) break; // stop after we start getting into the softer candidates
        if(pfCand.eta()<-2.5 || pfCand.eta()>2.5) continue;
        probePid=(pfCand.ptype== panda::PFCand::hp)? 211:-211; // PDG code for pi+/pi-.
        *p4_probe=pfCand.p4();

        // Clean the CH of leptons in DR cone 0.3
        bool chIsBadLepton=false;
        if (event.muons.size()>0) for (unsigned iM = 0; iM != event.muons.size(); ++iM){
          if(p4_probe->DeltaR(((panda::Muon)event.muons[iM]).p4())<0.3) {chIsBadLepton=true; break;}
        } if(chIsBadLepton) continue;
        if (event.electrons.size()>0) for (unsigned iE = 0; iE != event.electrons.size(); ++iE){
          if(p4_probe->DeltaR(((panda::Electron)event.electrons[iE]).p4())<0.3) {chIsBadLepton=true; break;}
        } if(chIsBadLepton) continue;

        // Muon tags
        if(do_muons) for (unsigned int iTag=0; iTag < p4_mu_tag_.size(); iTag++) {
          pass = 0; // all pf candidates fail lepton ID
          *p4_tag     = p4_mu_tag_[iTag];
          qtag    = q_mu_tag_[iTag];
          qprobe  = (pfCand.ptype== panda::PFCand::hp)? 1:-1;
          tagPid = -13*q_mu_tag_[iTag];
          truth_tag    = truth_mu_tag_[iTag];
          truth_probe  = 0; //dummy
          if( p4_tag->DeltaR(*p4_probe) < .1) continue;
          if (!real_data) scale1fb = event.weight*puWeights->GetBinContent(event.npvTrue);
          TLorentzVector systemP4 = (*p4_tag) + (*p4_probe);
          mass = systemP4.M();
          if(mass<40.||mass>140.) continue;
          if (verbose) printf("\t\tmade a Mu-Pi pair! pTs %f, %f; system mass %f, total charge %d e\n", p4_tag->Pt(), p4_probe->Pt(), mass, qtag+qprobe);
          muon_pair_tree->Fill();
        }
        //Electron tags
        if(do_electrons) for (unsigned int iTag=0; iTag < p4_ele_tag_.size(); iTag++) {
          pass = 0; // all pf candidates fail lepton ID
          *p4_tag     = p4_ele_tag_[iTag];
          qtag    = q_ele_tag_[iTag];
          qprobe  = (pfCand.ptype== panda::PFCand::hp)? 1:-1;
          tagPid = -11*q_ele_tag_[iTag];
          truth_tag    = truth_ele_tag_[iTag];
          truth_probe  = 0; //dummy
          if( p4_tag->DeltaR(*p4_probe) < .1) continue;
          if (!real_data) scale1fb = event.weight*puWeights->GetBinContent(event.npvTrue);
          TLorentzVector systemP4 = (*p4_tag) + (*p4_probe);
          mass = systemP4.M();
          if(mass<40.||mass>140.) continue;
          if (verbose) printf("\t\tmade a Mu-Pi pair! pTs %f, %f; system mass %f, total charge %d e\n", p4_tag->Pt(), p4_probe->Pt(), mass, qtag+qprobe);
          electron_pair_tree->Fill();
        }
        iPF++;
      }
    }
  }
  input_file->Close();

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

  delete puWeights;
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
  
