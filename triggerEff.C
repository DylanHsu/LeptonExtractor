#include "triggerEff.h"
#include "PandaTree/Objects/interface/Event.h"

typedef std::map<UInt_t,std::vector<std::pair <UInt_t, UInt_t> > > MapType;
string jsonFile = "PandaAnalysis/data/certs/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON_v1.txt";

using namespace panda;

void triggerEff(
  TString trigType = "OneSMTrig", // "SMTrigSoup","DMRefTrig","DMTrigSoup"
  TString inputFileName="",
  TString outputFileName="",
  bool debug=false
) {
  // Register the trigger tokens for the reference and test triggers
  // Empty filter means we only care if the event passes the trigger
  // TO DO: special hack for Ele32
  vector<unsigned> refTriggerTokens, testTriggerTokens;
  vector<pair<TString,TString>> refTriggersAndFilters,testTriggersAndFilters;
  if (trigType=="OneSMTrig") {
    refTriggersAndFilters = {
      {"HLT_IsoMu27", "hltL3crIsoL1sMu22Or25L1f0L2f10QL3f27QL3trkIsoFiltered0p07"}
    };
    testTriggersAndFilters = {
      {"HLT_IsoMu27", "hltL3crIsoL1sMu22Or25L1f0L2f10QL3f27QL3trkIsoFiltered0p07"}
    };
  } else if (trigType=="SMTrigSoup") {
    refTriggersAndFilters = {
      {"HLT_IsoMu27", "hltL3crIsoL1sMu22Or25L1f0L2f10QL3f27QL3trkIsoFiltered0p07"}
    };
    testTriggersAndFilters = {
      {"HLT_IsoMu24", "hltL3crIsoL1sSingleMu22erL1f0L2f10QL3f24QL3trkIsoFiltered0p07"},
      {"HLT_IsoMu27", "hltL3crIsoL1sMu22Or25L1f0L2f10QL3f27QL3trkIsoFiltered0p07"    },
      {"HLT_IsoMu30", "hltL3crIsoL1sMu22Or25L1f0L2f10QL3f30QL3trkIsoFiltered0p07"    },
      {"HLT_Mu50"   , "hltL3fL1sMu22Or25L1f0L2f10QL3Filtered50Q"                     },
    };
  } else if (trigType=="DMRefTrig") {
    refTriggersAndFilters = {
      {"HLT_IsoMu27", "hltL3crIsoL1sMu22Or25L1f0L2f10QL3f27QL3trkIsoFiltered0p07"}
    };
    testTriggersAndFilters = {
      {"HLT_Mu20"   , "hltL3fL1sMu18L1f0L2f10QL3Filtered20Q"}
    };
  } else if (trigType=="DMTrigSoup") {
    refTriggersAndFilters = {
      //{"HLT_Mu20"   , "hltL3fL1sMu18L1f0L2f10QL3Filtered20Q"},
      //{"HLT_Mu20"   , "hltL1fForIterL3L1fL1sMu18L1Filtered0"},
      //{"HLT_Mu20"   , "hltL2fL1sMu18L1f0L2Filtered10Q"}
      {"HLT_IsoMu24", "hltL3crIsoL1sSingleMu22erL1f0L2f10QL3f24QL3trkIsoFiltered0p07"},
    };
    testTriggersAndFilters = {
      {"HLT_IsoMu24", ""},
      {"HLT_IsoMu27", ""},
      {"HLT_IsoMu30", ""},
      {"HLT_Mu50"   , ""},
      {"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8", ""},
      {"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8"  , ""},
	  {"HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8", ""},
	  {"HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8"  , ""}
    };
  } else if (trigType=="OneSETrig") {
    refTriggersAndFilters = {
      {"HLT_Ele32_WPTight_Gsf",""} // special case
    };
    testTriggersAndFilters = {
      {"HLT_Ele32_WPTight_Gsf", ""}
    };
  } else if (trigType=="SETrigSoup") {
    refTriggersAndFilters = {
      {"HLT_Ele32_WPTight_Gsf",""} // special case
    };
    testTriggersAndFilters = {
      {"HLT_Ele115_CaloIdVT_GsfTrkIdT"   , "hltEle115CaloIdVTGsfTrkIdTGsfDphiFilter"   },  
      {"HLT_Ele27_WPTight_Gsf"           , "hltEle27WPTightGsfTrackIsoFilter"          },
      {"HLT_Ele32_WPTight_Gsf"           , ""                                          },
      {"HLT_Ele35_WPTight_Gsf"           , "hltEle35noerWPTightGsfTrackIsoFilter"      },
      {"HLT_Ele32_WPTight_Gsf_L1DoubleEG", "hltEle32L1DoubleEGWPTightGsfTrackIsoFilter"},
      {"HLT_Photon200"                   , "hltEG200HEFilter"                          },
    };
  } else if (trigType=="DERefTrig") {
    refTriggersAndFilters = {
      {"HLT_Ele32_WPTight_Gsf",""} // special case
    };
    testTriggersAndFilters = {
      {"HLT_Ele27_WPTight_Gsf"           , "hltEle27WPTightGsfTrackIsoFilter"     },
    };
  } else if (trigType=="DETrigSoup") {
    refTriggersAndFilters = {
      {"HLT_Ele27_WPTight_Gsf"           , "hltEle27WPTightGsfTrackIsoFilter"     },
    };
    testTriggersAndFilters = {
      {"HLT_Ele115_CaloIdVT_GsfTrkIdT"             , ""}, 
      {"HLT_Ele27_WPTight_Gsf"                     , ""}, 
      {"HLT_Ele32_WPTight_Gsf"                     , ""}, 
      {"HLT_Ele35_WPTight_Gsf"                     , ""}, 
      {"HLT_Ele32_WPTight_Gsf_L1DoubleEG"          , ""}, 
      {"HLT_Photon200"                             , ""}, 
      {"HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ" , ""}, 
      {"HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL"    , ""}, 
      {"HLT_DiEle27_WPTightCaloOnly_L1DoubleEG"    , ""}, 
      {"HLT_DoubleEle33_CaloIdL_MW"                , ""}, 
      {"HLT_DoubleEle25_CaloIdL_MW"                , ""}, 
      {"HLT_DoublePhoton70"                        , ""}, 
    };
  } else {
    printf("trigType not supported\n");
    return; 
  }

  UInt_t runNumber, lumiSection;
  ULong64_t eventNumber;
  //Read json file into boost property tree
  MapType fMap;
  boost::property_tree::ptree jsonTree;
  boost::property_tree::read_json(jsonFile.c_str(),jsonTree);
  
  //Loop through boost property tree and fill the MapType structure with the list of good lumi
  //ranges for each run
  for (boost::property_tree::ptree::const_iterator it = jsonTree.begin(); it!=jsonTree.end(); ++it) {
    runNumber = boost::lexical_cast<UInt_t>(it->first);
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
  //If running in debug mode, dump run and lumi ranges from MapType structure to verify correct json parsing
  if (debug) {  printf("Iterating over parsed JSON:\n"); for (MapType::const_iterator it = fMap.begin(); it != fMap.end(); ++it) {
    printf("  Run %u:\n",it->first);
    for (MapType::mapped_type::const_iterator jt = it->second.begin(); jt < it->second.end(); ++jt) printf("    Lumis %u - %u\n",jt->first,jt->second);
  }}
  
  // Open the file
  TFile *inputFile=0;
  int retries=0;
  while(true) {
    inputFile = TFile::Open(inputFileName,"read");
    if(inputFile && inputFile->IsOpen()) break;
    usleep(2e5);
    retries++;
    if(retries>100) { throw std::runtime_error("Error opening input file"); return; }
  }
  TTree* tree = (TTree*)inputFile->Get("events"); // get the tree object from the file
  
  // Initialize the Panda event object
  Event event;
  event.setStatus(*tree, {"!*"});
  event.setAddress(*tree, {"chsAK4Jets", "electrons", "muons", "pfMet", "triggers", "triggerObjects", "runNumber", "lumiNumber", "eventNumber"});
  //TString outputFileName =  outputDir + inputFileName(inputFileName.Last('/')+1,inputFileName.Length())
  TFile *outputFile = TFile::Open(outputFileName, "RECREATE", "", ROOT::CompressionSettings(ROOT::kZLIB,9));
  // Declare the Tree and set up the branches.
 
  // Declare the simple variables for the tree branches
  bool passTrigger;
  float tagPt, tagEta, tagPhi, probePt, probeEta, probePhi, mass;
  int tagPdgId, probePdgId;
  
  TTree *effTree = new TTree("effTree", "effTree");
  effTree->Branch("runNumber"          , &runNumber         );
  effTree->Branch("lumiSection"        , &lumiSection       );
  effTree->Branch("eventNumber"        , &eventNumber       );
  effTree->Branch("passTrigger"        , &passTrigger       );
  effTree->Branch("tagPt"              , &tagPt             );
  effTree->Branch("tagEta"             , &tagEta            );
  effTree->Branch("tagPhi"             , &tagPhi            );
  effTree->Branch("tagPdgId"           , &tagPdgId          );
  effTree->Branch("probePt"            , &probePt           );
  effTree->Branch("probeEta"           , &probeEta          );
  effTree->Branch("probePhi"           , &probePhi          );
  effTree->Branch("probePdgId"         , &probePdgId        );
  effTree->Branch("mass"               , &mass              );
  if(debug) printf("Initialized the tree successfully\n");
  
  // Register the triggers
  for(auto const &triggerAndFilter: refTriggersAndFilters)
    refTriggerTokens.push_back( event.registerTrigger( triggerAndFilter.first.Data() ));
  for(auto const &triggerAndFilter: testTriggersAndFilters)
    testTriggerTokens.push_back( event.registerTrigger( triggerAndFilter.first.Data()));
  
  // Loop over all of the events in the TChain
  long iEntry = 0;
  while (event.getEntry(*tree, iEntry++) > 0) {
    if(debug) printf("Run %u, LS %u, evt %lld:\n", event.runNumber, event.lumiNumber, event.eventNumber);
    
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
    if(!certifiedEvent) { if(debug) printf("failed golden json\n"); continue; }

    // Require reference trigger
    bool passRefTriggers=false;
    for(unsigned i=0; i<refTriggerTokens.size() && !passRefTriggers; i++) 
      if(event.triggerFired(refTriggerTokens[i]) || 
         refTriggersAndFilters[i].first=="HLT_Ele32_WPTight_Gsf"
      ) passRefTriggers=true;
    if(!passRefTriggers) {  if(debug) printf("failed reference triggers\n"); continue; }

    // Begin finding leptons
    vector<Lepton*> tightLeps, refMatchedTightLeps; // tight leptons

    // Find tight electrons
    if(
      trigType == "OneSETrig" ||
      trigType == "SETrigSoup"||
      trigType == "DERefTrig" ||
      trigType == "DETrigSoup"
    ) for (auto& ele : event.electrons) {
     float pt = ele.pt(); float eta = ele.eta(); float aeta = fabs(eta);
      if (pt<=10 || aeta>2.5)// || (aeta>1.4442 && aeta<1.566)) // electron acceptance cuts
        continue;
      if(!ele.tight) continue;
      if (!ElectronIP(ele.eta(),ele.dxy,ele.dz)) continue;
      
      tightLeps.push_back(&ele);
      bool isMatched=false;
      for(unsigned i=0; i<refTriggerTokens.size() && !isMatched; i++) {
        if(refTriggersAndFilters[i].first=="HLT_Ele32_WPTight_Gsf") {
          // Special case for Ele32
          if(debug) printf("checking Ele32\n");
          if(checkEle32(&event, &ele)==true)
            isMatched=true;
        } else {
          if(!event.triggerFired(refTriggerTokens[i])) continue;
          if(matchLepToFilter(&event, &ele, refTriggersAndFilters[i].second.Data()))
            isMatched=true;
        }
      }
      if(isMatched) refMatchedTightLeps.push_back(&ele);
    }
    // Find tight muons
    
    if(
      trigType=="OneSMTrig" ||
      trigType=="SMTrigSoup"||
      trigType=="DMRefTrig" ||
      trigType=="DMTrigSoup"
    ) for (auto& mu : event.muons) {
      float pt = mu.pt(); float eta = mu.eta(); float aeta = fabs(eta);
      if (pt<=10 || aeta>2.4) continue; //muon acceptance cuts
      if(!mu.tight) continue;
      if(!MuonIsolation(mu.pt(),mu.eta(),mu.combIso(),kTight)) continue;

      bool isMatched=false;
      for(unsigned i=0; i<refTriggerTokens.size() && !isMatched; i++) {
        if(!event.triggerFired(refTriggerTokens[i])) continue;
        if(matchLepToFilter(&event, &mu, refTriggersAndFilters[i].second.Data()))
          isMatched=true;
      }
      tightLeps.push_back(&mu);
      if(isMatched) refMatchedTightLeps.push_back(&mu);
    }

    if(refMatchedTightLeps.size()==0) { if(debug) printf("failed reference trigger object matching\n"); continue; }
    
    runNumber = event.runNumber;
    lumiSection = event.lumiNumber;
    eventNumber = event.eventNumber;
    
    for(auto& tagLep: refMatchedTightLeps) for(auto& probeLep: tightLeps) {
      if(tagLep==probeLep) continue;
      passTrigger = false;
      unsigned iT=0;
      for(auto const &triggerAndFilter: testTriggersAndFilters) {
        if(triggerAndFilter.first=="HLT_Ele32_WPTight_Gsf") {
          if(checkEle32(&event, probeLep)==true) {
            passTrigger=true; break;
          }
        } else if(triggerAndFilter.second=="") { // Only care about the event passing
          if(event.triggerFired(testTriggerTokens[iT])) {
            passTrigger=true;
            break;
          }
        } else { // Have to match to the trigger
          bool matchedProbe = matchLepToFilter(&event, probeLep, triggerAndFilter.second.Data());
          if(matchedProbe) {
            passTrigger=true;
            break;
          }
        }
        iT++;
      }
      TLorentzVector dilep = tagLep->p4() + probeLep->p4();
      mass = dilep.M();
      if(mass<40) continue;
      if((
        trigType== "OneSMTrig"   ||
        trigType== "SMTrigSoup"  ||
        trigType== "DMRefTrig"   ||
        trigType== "DMTrigSoup"
       ) && (tagLep->charge + probeLep->charge !=0))
        continue;

      Electron* tagIsEle = dynamic_cast<Electron*>(tagLep);
      Electron* probeIsEle = dynamic_cast<Electron*>(probeLep);
      tagPt    = tagLep->pt();
      tagEta   = tagLep->eta();
      tagPhi   = tagLep->phi();
      tagPdgId = tagLep->charge * (tagIsEle? -11:-13);
      probePt    = probeLep->pt();
      probeEta   = probeLep->eta();
      probePhi   = probeLep->phi();
      probePdgId = probeLep->charge * (probeIsEle? -11:-13);
      if(debug) printf("Filling tree: tagPt=%.2f, tagPdgId=%d, probePt=%.2f, probePdgId=%d, pass=%d\n", tagPt, tagPdgId, probePt, probePdgId, passTrigger);
      effTree->Fill();
    }
  }
  outputFile->cd();
  effTree->Write();
  outputFile->Close();
}

bool matchLepToFilter(Event* event, Lepton* lepton, const char* filterName) {
  HLTObjectStore::HLTObjectVector objects = event->triggerObjects.filterObjects(filterName);
  //if(objects.size()==0) printf("Warning: HLTObjectVector is empty for filter \"%s\"\n", filterName);
  for(auto& object : objects)
    if(lepton->p4().DeltaR(object->p4())<0.1) return true;
  return false;
}

bool checkEle32(Event* event, Lepton *lepton) {
  HLTObjectStore::HLTObjectVector filter1Objects, filter2Objects;
  filter1Objects=event->triggerObjects.filterObjects(ele32Filter1);
  filter2Objects=event->triggerObjects.filterObjects(ele32Filter2);
  if(filter1Objects.size()==0 || filter2Objects.size()==0) {
    return false;
  }
  bool matchFilter1=false, matchFilter2=false;
  for(auto& object : filter1Objects)
    if(lepton->p4().DeltaR(object->p4())<0.1) { matchFilter1=true; break; }
  for(auto& object : filter2Objects)
    if(lepton->p4().DeltaR(object->p4())<0.1) { matchFilter2=true; break; }
  return matchFilter1 && matchFilter2;
}