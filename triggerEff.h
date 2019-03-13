#include <map>

// Root libraries
#include <Compression.h>
#include <TChain.h>
#include <TFile.h>
#include <TFileCollection.h>
#include <TString.h>
#include <TTree.h>

// Panda
#include "PandaCore/Tools/interface/Common.h"
#include "PandaCore/Tools/interface/DataTools.h"
#include "PandaCore/Tools/interface/JERReader.h"
#include "PandaTree/Objects/interface/Event.h"
#include "PandaTree/Objects/interface/HLTObjectStore.h"

//Boost
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/lexical_cast.hpp>

using namespace panda;

namespace panda {
  enum IDWorkingPoint {
    kVeto,
    kLoose,
    kMedium,
    kTight,
    nIDWorkingPoints
  };
}

inline bool MuonIsolation(double pt, double eta, double iso, panda::IDWorkingPoint isoType) {
    float maxIso=0;
    maxIso = (isoType == panda::kTight) ? 0.15 : 0.25;
    return (iso < pt*maxIso);
}

inline bool ElectronIP(double eta, double dxy, double dz) {
  double aeta = fabs(eta);
  if (aeta<1.4442) {
    return (dxy < 0.05 && dz < 0.10) ;
  } else {
    return (dxy < 0.10 && dz < 0.20);
  }
}

const char* ele32Filter1="hltEle32L1DoubleEGWPTightGsfTrackIsoFilter";
const char* ele32Filter2="hltEGL1SingleEGOrFilter";

bool matchLepToFilter(Event *event, Lepton* lepton, const char* filterName, bool debug=false);

bool checkEle32(Event *event, Lepton* lepton, bool debug=false);
