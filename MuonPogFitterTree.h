//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Feb  8 18:26:19 2018 by ROOT version 6.06/01
// from TTree fitter_tree/fitter_tree
// found on file: /home/dhsu/muonPOG/CMSSW_8_0_29/src/MuonAnalysis/TagAndProbe/test/zmumu/tnpZ_MC.root
//////////////////////////////////////////////////////////

#ifndef MuonPogFitterTree_h
#define MuonPogFitterTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class MuonPogFitterTree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Float_t         IP;
   Float_t         IPError;
   Float_t         SIP;
   Float_t         abseta;
   Float_t         caloCompatibility;
   Float_t         charge;
   Float_t         chargedHadIso03;
   Float_t         chargedHadIso04;
   Float_t         chargedParticleIso03;
   Float_t         chargedParticleIso04;
   Float_t         chi2LocMom;
   Float_t         chi2LocPos;
   Float_t         combRelIso;
   Float_t         combRelIsoPF03;
   Float_t         combRelIsoPF03dBeta;
   Float_t         combRelIsoPF04;
   Float_t         combRelIsoPF04dBeta;
   Float_t         dB;
   Float_t         ecalIso;
   Float_t         edB;
   Float_t         emEnergy;
   Float_t         emS9Energy;
   Float_t         eta;
   Float_t         glbChi2;
   Float_t         glbPtError;
   Float_t         glbSigmaPtOverPt;
   Float_t         glbTrackProb;
   Float_t         glbValidMuHits;
   Float_t         hadEnergy;
   Float_t         hadS9Energy;
   Float_t         hcalIso;
   Float_t         l1dphi;
   Float_t         l1dr;
   Float_t         l1drByQ;
   Float_t         l1eta;
   Float_t         l1phi;
   Float_t         l1pt;
   Float_t         l1ptByQ;
   Float_t         l1q;
   Float_t         l1qByQ;
   Float_t         l2dr;
   Float_t         l2eta;
   Float_t         l2pt;
   Float_t         l3dr;
   Float_t         l3pt;
   Float_t         neutralHadIso03;
   Float_t         neutralHadIso04;
   Float_t         numberOfMatchedStations;
   Float_t         numberOfMatches;
   Float_t         p;
   Float_t         phi;
   Float_t         photonIso03;
   Float_t         photonIso04;
   Float_t         pt;
   Float_t         puIso03;
   Float_t         puIso04;
   Float_t         relEcalIso;
   Float_t         relHcalIso;
   Float_t         relTkIso;
   Float_t         segmentCompatibility;
   Float_t         staQoverP;
   Float_t         staQoverPerror;
   Float_t         staValidStations;
   Float_t         tkChi2;
   Float_t         tkExpHitIn;
   Float_t         tkExpHitOut;
   Float_t         tkHitFract;
   Float_t         tkIso;
   Float_t         tkKink;
   Float_t         tkPixelLay;
   Float_t         tkPtError;
   Float_t         tkSigmaPtOverPt;
   Float_t         tkTrackerLay;
   Float_t         tkValidHits;
   Float_t         tkValidPixelHits;
   Float_t         JetBTagCSV;
   Float_t         JetNDauCharged;
   Float_t         JetPtRatio;
   Float_t         JetPtRel;
   Float_t         activity_miniIsoCharged;
   Float_t         activity_miniIsoNeutrals;
   Float_t         activity_miniIsoPUCharged;
   Float_t         activity_miniIsoPhotons;
   Float_t         dxyBS;
   Float_t         dxyPVdzmin;
   Float_t         dzPV;
   Float_t         fixedGridRhoFastjetAll;
   Float_t         fixedGridRhoFastjetAllCalo;
   Float_t         fixedGridRhoFastjetCentralCalo;
   Float_t         fixedGridRhoFastjetCentralChargedPileUp;
   Float_t         fixedGridRhoFastjetCentralNeutral;
   Float_t         isoTrk03Abs;
   Float_t         isoTrk03Rel;
   Float_t         miniIsoCharged;
   Float_t         miniIsoNeutrals;
   Float_t         miniIsoPUCharged;
   Float_t         miniIsoPhotons;
   Float_t         mt;
   Float_t         muPFIsoValueCHR04PUPPI;
   Float_t         muPFIsoValueCHR04PUPPINoLep;
   Float_t         muPFIsoValueNHR04PUPPI;
   Float_t         muPFIsoValueNHR04PUPPINoLep;
   Float_t         muPFIsoValuePhR04PUPPI;
   Float_t         muPFIsoValuePhR04PUPPINoLep;
   Float_t         nSplitTk;
   Int_t           Calo;
   Int_t           DiMuonGlb17Glb8RelTrkIsoFiltered0p4;
   Int_t           DoubleIsoMu17Mu8_IsoMu17leg;
   Int_t           DoubleIsoMu17Mu8_IsoMu8leg;
   Int_t           DoubleIsoMu17Mu8_Mu17leg;
   Int_t           DoubleIsoMu17Mu8_Mu8leg;
   Int_t           DoubleIsoMu17Mu8dZ_Mu17leg;
   Int_t           DoubleIsoMu17TkMu8_IsoMu17leg;
   Int_t           DoubleIsoMu17TkMu8_IsoMu8leg;
   Int_t           DoubleIsoMu17TkMu8_Mu17leg;
   Int_t           DoubleIsoMu17TkMu8_TkMu8leg;
   Int_t           DoubleIsoMu17TkMu8dZ_Mu17;
   Int_t           DoubleMu30TkMu11;
   Int_t           DoubleMu30TkMu11_Mu30leg;
   Int_t           DoubleMu30TkMu11_TkMu11leg;
   Int_t           DoubleIsoMu17Mu8dZ_Mass3p8;
   Int_t           DoubleIsoMu17Mu8dZ_Mass8;
   Int_t           Glb;
   Int_t           GlbPT;
   Int_t           HLT_TkMu50;
   Int_t           HWWID;
   Int_t           HighPt;
   Int_t           IsoMu18;
   Int_t           IsoMu20;
   Int_t           IsoMu22;
   Int_t           IsoMu22_eta2p1;
   Int_t           IsoMu24;
   Int_t           IsoMu24_eta2p1;
   Int_t           IsoMu27;
   Int_t           IsoTkMu18;
   Int_t           IsoTkMu20;
   Int_t           IsoTkMu22;
   Int_t           IsoTkMu22_eta2p1;
   Int_t           IsoTkMu24;
   Int_t           IsoTkMu24_eta2p1;
   Int_t           IsoTkMu27;
   Int_t           L1sMu16;
   Int_t           L1sMu18;
   Int_t           L1sMu20;
   Int_t           L1sMu25Eta2p1;
   Int_t           L2fL1sDoubleMu114L1f0L2Filtered10OneMu;
   Int_t           L2fL1sDoubleMu114L1f0OneMuL2Filtered10;
   Int_t           L2fL1sMu18L1f0L2Filtered10Q;
   Int_t           L2pfL1sDoubleMu114L1f0L2PreFiltered0;
   Int_t           L3crIsoL1sMu16L1f0L2f10QL3f18QL3trkIsoFiltered0p09;
   Int_t           L3fL1sDoubleMu114L1f0L2f10L3Filtered17;
   Int_t           L3fL1sDoubleMu114L1f0L2f10OneMuL3Filtered17;
   Int_t           L3fL1sMu16L1f0L2f10QL3Filtered18Q;
   Int_t           L3fL1sMu16L1f0L2f10QL3Filtered18QL3pfecalIsoRhoFilteredEB0p11EE0p08;
   Int_t           L3fL1sMu16L1f0L2f10QL3Filtered18QL3pfhcalIsoRhoFilteredHB0p21HE0p22;
   Int_t           L3fL1sMu16f0TkFiltered18Q;
   Int_t           L3pfL1sDoubleMu114L1f0L2pf0L3PreFiltered8;
   Int_t           Loose;
   Int_t           Medium;
   Int_t           Medium2016;
   Int_t           Mu17;
   Int_t           Mu17_IsoTrkVVL;
   Int_t           Mu20;
   Int_t           Mu23_TrkIsoVVL;
   Int_t           Mu45_eta2p1;
   Int_t           Mu50;
   Int_t           MuIDForOutsideInTk;
   Int_t           PF;
   Int_t           TM;
   Int_t           TMA;
   Int_t           TMLSAT;
   Int_t           TMLST;
   Int_t           TMOSL;
   Int_t           TMOST;
   Int_t           TMOSTQual;
   Int_t           Tight2012;
   Int_t           Track_HP;
   Int_t           VBTF;
   Int_t           VBTF_nL8;
   Int_t           VBTF_nL9;
   Int_t           hltL2fL1sMu10lqL1f0L2Filtered10;
   Int_t           hltL2fL1sMu20L1f0L2Filtered10Q;
   Int_t           hltL2fL1sMu22L1f0L2Filtered10Q;
   Int_t           hltL3crIsoL1sMu18L1f0L2f10QL3f20QL3trkIsoFiltered0p09;
   Int_t           hltL3crIsoL1sMu20L1f0L2f10QL3f22QL3trkIsoFiltered0p09;
   Int_t           hltL3crIsoL1sMu22L1f0L2f10QL3f24QL3trkIsoFiltered0p09;
   Int_t           hltL3fL1sMu10lqL1f0L2f10L3Filtered17;
   Int_t           hltL3fL1sMu18L1f0L2f10QL3Filtered20Q;
   Int_t           hltL3fL1sMu18L1f0L2f10QL3Filtered20QL3pfecalIsoRhoFilteredEB0p11EE0p08;
   Int_t           hltL3fL1sMu18L1f0L2f10QL3Filtered20QL3pfhcalIsoRhoFilteredHB0p21HE0p22;
   Int_t           hltL3fL1sMu18L1f0Tkf20QL3trkIsoFiltered0p09;
   Int_t           hltL3fL1sMu18f0TkFiltered20Q;
   Int_t           hltL3fL1sMu18f0TkFiltered20QL3pfecalIsoRhoFilteredEB0p11EE0p08;
   Int_t           hltL3fL1sMu18f0TkFiltered20QL3pfhcalIsoRhoFilteredHB0p21HE0p22;
   Int_t           hltL3fL1sMu1lqL1f0L2f10L3Filtered17TkIsoFiltered0p4;
   Int_t           hltL3fL1sMu20L1f0L2f10QL3Filtered22Q;
   Int_t           hltL3fL1sMu20L1f0L2f10QL3Filtered22QL3pfecalIsoRhoFilteredEB0p11EE0p08;
   Int_t           hltL3fL1sMu20L1f0L2f10QL3Filtered22QL3pfhcalIsoRhoFilteredHB0p21HE0p22;
   Int_t           hltL3fL1sMu20L1f0Tkf22QL3trkIsoFiltered0p09;
   Int_t           hltL3fL1sMu20f0TkFiltered22Q;
   Int_t           hltL3fL1sMu20f0TkFiltered22QL3pfecalIsoRhoFilteredEB0p11EE0p08;
   Int_t           hltL3fL1sMu20f0TkFiltered22QL3pfhcalIsoRhoFilteredHB0p21HE0p22;
   Int_t           hltL3fL1sMu22L1f0L2f10QL3Filtered24Q;
   Int_t           hltL3fL1sMu22L1f0L2f10QL3Filtered24QL3pfecalIsoRhoFilteredEB0p11EE0p08;
   Int_t           hltL3fL1sMu22L1f0L2f10QL3Filtered24QL3pfhcalIsoRhoFilteredHB0p21HE0p22;
   Int_t           hltL3fL1sMu22L1f0Tkf24QL3trkIsoFiltered0p09;
   Int_t           hltL3fL1sMu22f0TkFiltered24Q;
   Int_t           hltL3fL1sMu22f0TkFiltered24QL3pfecalIsoRhoFilteredEB0p11EE0p08;
   Int_t           hltL3fL1sMu22f0TkFiltered24QL3pfhcalIsoRhoFilteredHB0p21HE0p22;
   Int_t           tkHighPt;
   UInt_t          run;
   UInt_t          lumi;
   ULong64_t       event;
   Int_t           truePU;
   Float_t         mass;
   Int_t           mcTrue;
   Float_t         mcMass;
   Float_t         tag_IP;
   Float_t         tag_IPError;
   Float_t         tag_SIP;
   Float_t         tag_abseta;
   Float_t         tag_caloCompatibility;
   Float_t         tag_charge;
   Float_t         tag_chargedHadIso03;
   Float_t         tag_chargedHadIso04;
   Float_t         tag_chargedParticleIso03;
   Float_t         tag_chargedParticleIso04;
   Float_t         tag_chi2LocMom;
   Float_t         tag_chi2LocPos;
   Float_t         tag_combRelIso;
   Float_t         tag_combRelIsoPF03;
   Float_t         tag_combRelIsoPF03dBeta;
   Float_t         tag_combRelIsoPF04;
   Float_t         tag_combRelIsoPF04dBeta;
   Float_t         tag_dB;
   Float_t         tag_ecalIso;
   Float_t         tag_edB;
   Float_t         tag_emEnergy;
   Float_t         tag_emS9Energy;
   Float_t         tag_eta;
   Float_t         tag_glbChi2;
   Float_t         tag_glbPtError;
   Float_t         tag_glbSigmaPtOverPt;
   Float_t         tag_glbTrackProb;
   Float_t         tag_glbValidMuHits;
   Float_t         tag_hadEnergy;
   Float_t         tag_hadS9Energy;
   Float_t         tag_hcalIso;
   Float_t         tag_l1dphi;
   Float_t         tag_l1dr;
   Float_t         tag_l1drByQ;
   Float_t         tag_l1eta;
   Float_t         tag_l1phi;
   Float_t         tag_l1pt;
   Float_t         tag_l1ptByQ;
   Float_t         tag_l1q;
   Float_t         tag_l1qByQ;
   Float_t         tag_l2dr;
   Float_t         tag_l2eta;
   Float_t         tag_l2pt;
   Float_t         tag_l3dr;
   Float_t         tag_l3pt;
   Float_t         tag_neutralHadIso03;
   Float_t         tag_neutralHadIso04;
   Float_t         tag_numberOfMatchedStations;
   Float_t         tag_numberOfMatches;
   Float_t         tag_p;
   Float_t         tag_phi;
   Float_t         tag_photonIso03;
   Float_t         tag_photonIso04;
   Float_t         tag_pt;
   Float_t         tag_puIso03;
   Float_t         tag_puIso04;
   Float_t         tag_relEcalIso;
   Float_t         tag_relHcalIso;
   Float_t         tag_relTkIso;
   Float_t         tag_segmentCompatibility;
   Float_t         tag_staQoverP;
   Float_t         tag_staQoverPerror;
   Float_t         tag_staValidStations;
   Float_t         tag_tkChi2;
   Float_t         tag_tkExpHitIn;
   Float_t         tag_tkExpHitOut;
   Float_t         tag_tkHitFract;
   Float_t         tag_tkIso;
   Float_t         tag_tkKink;
   Float_t         tag_tkPixelLay;
   Float_t         tag_tkPtError;
   Float_t         tag_tkSigmaPtOverPt;
   Float_t         tag_tkTrackerLay;
   Float_t         tag_tkValidHits;
   Float_t         tag_tkValidPixelHits;
   Float_t         tag_dxyBS;
   Float_t         tag_dxyPVdzmin;
   Float_t         tag_dzPV;
   Float_t         tag_fixedGridRhoFastjetAll;
   Float_t         tag_fixedGridRhoFastjetAllCalo;
   Float_t         tag_fixedGridRhoFastjetCentralCalo;
   Float_t         tag_fixedGridRhoFastjetCentralChargedPileUp;
   Float_t         tag_fixedGridRhoFastjetCentralNeutral;
   Float_t         tag_isoTrk03Abs;
   Float_t         tag_isoTrk03Rel;
   Float_t         tag_met;
   Float_t         tag_mt;
   Float_t         tag_nSplitTk;
   Float_t         tag_nVertices;
   Int_t           tag_DiMuonGlb17Glb8RelTrkIsoFiltered0p4;
   Int_t           tag_DoubleIsoMu17Mu8_IsoMu17leg;
   Int_t           tag_DoubleIsoMu17Mu8_IsoMu8leg;
   Int_t           tag_DoubleIsoMu17Mu8_Mu17leg;
   Int_t           tag_DoubleIsoMu17Mu8_Mu8leg;
   Int_t           tag_DoubleIsoMu17Mu8dZ_Mu17leg;
   Int_t           tag_DoubleIsoMu17TkMu8_IsoMu17leg;
   Int_t           tag_DoubleIsoMu17TkMu8_IsoMu8leg;
   Int_t           tag_DoubleIsoMu17TkMu8_Mu17leg;
   Int_t           tag_DoubleIsoMu17TkMu8_TkMu8leg;
   Int_t           tag_DoubleIsoMu17TkMu8dZ_Mu17;
   Int_t           tag_DoubleMu30TkMu11;
   Int_t           tag_DoubleMu30TkMu11_Mu30leg;
   Int_t           tag_DoubleMu30TkMu11_TkMu11leg;
   Int_t           tag_DoubleIsoMu17Mu8dZ_Mass3p8;
   Int_t           tag_DoubleIsoMu17Mu8dZ_Mass8;
   Int_t           tag_HLT_TkMu50;
   Int_t           tag_IsoMu18;
   Int_t           tag_IsoMu20;
   Int_t           tag_IsoMu22;
   Int_t           tag_IsoMu22_eta2p1;
   Int_t           tag_IsoMu24;
   Int_t           tag_IsoMu24_eta2p1;
   Int_t           tag_IsoMu27;
   Int_t           tag_IsoTkMu18;
   Int_t           tag_IsoTkMu20;
   Int_t           tag_IsoTkMu22;
   Int_t           tag_IsoTkMu22_eta2p1;
   Int_t           tag_IsoTkMu24;
   Int_t           tag_IsoTkMu24_eta2p1;
   Int_t           tag_IsoTkMu27;
   Int_t           tag_L1sMu16;
   Int_t           tag_L1sMu18;
   Int_t           tag_L1sMu20;
   Int_t           tag_L1sMu25Eta2p1;
   Int_t           tag_L2fL1sDoubleMu114L1f0L2Filtered10OneMu;
   Int_t           tag_L2fL1sDoubleMu114L1f0OneMuL2Filtered10;
   Int_t           tag_L2fL1sMu18L1f0L2Filtered10Q;
   Int_t           tag_L2pfL1sDoubleMu114L1f0L2PreFiltered0;
   Int_t           tag_L3crIsoL1sMu16L1f0L2f10QL3f18QL3trkIsoFiltered0p09;
   Int_t           tag_L3fL1sDoubleMu114L1f0L2f10L3Filtered17;
   Int_t           tag_L3fL1sDoubleMu114L1f0L2f10OneMuL3Filtered17;
   Int_t           tag_L3fL1sMu16L1f0L2f10QL3Filtered18Q;
   Int_t           tag_L3fL1sMu16L1f0L2f10QL3Filtered18QL3pfecalIsoRhoFilteredEB0p11EE0p08;
   Int_t           tag_L3fL1sMu16L1f0L2f10QL3Filtered18QL3pfhcalIsoRhoFilteredHB0p21HE0p22;
   Int_t           tag_L3fL1sMu16f0TkFiltered18Q;
   Int_t           tag_L3pfL1sDoubleMu114L1f0L2pf0L3PreFiltered8;
   Int_t           tag_Mu17;
   Int_t           tag_Mu17_IsoTrkVVL;
   Int_t           tag_Mu20;
   Int_t           tag_Mu23_TrkIsoVVL;
   Int_t           tag_Mu45_eta2p1;
   Int_t           tag_Mu50;
   Int_t           tag_hltL2fL1sMu10lqL1f0L2Filtered10;
   Int_t           tag_hltL2fL1sMu20L1f0L2Filtered10Q;
   Int_t           tag_hltL2fL1sMu22L1f0L2Filtered10Q;
   Int_t           tag_hltL3crIsoL1sMu18L1f0L2f10QL3f20QL3trkIsoFiltered0p09;
   Int_t           tag_hltL3crIsoL1sMu20L1f0L2f10QL3f22QL3trkIsoFiltered0p09;
   Int_t           tag_hltL3crIsoL1sMu22L1f0L2f10QL3f24QL3trkIsoFiltered0p09;
   Int_t           tag_hltL3fL1sMu10lqL1f0L2f10L3Filtered17;
   Int_t           tag_hltL3fL1sMu18L1f0L2f10QL3Filtered20Q;
   Int_t           tag_hltL3fL1sMu18L1f0L2f10QL3Filtered20QL3pfecalIsoRhoFilteredEB0p11EE0p08;
   Int_t           tag_hltL3fL1sMu18L1f0L2f10QL3Filtered20QL3pfhcalIsoRhoFilteredHB0p21HE0p22;
   Int_t           tag_hltL3fL1sMu18L1f0Tkf20QL3trkIsoFiltered0p09;
   Int_t           tag_hltL3fL1sMu18f0TkFiltered20Q;
   Int_t           tag_hltL3fL1sMu18f0TkFiltered20QL3pfecalIsoRhoFilteredEB0p11EE0p08;
   Int_t           tag_hltL3fL1sMu18f0TkFiltered20QL3pfhcalIsoRhoFilteredHB0p21HE0p22;
   Int_t           tag_hltL3fL1sMu1lqL1f0L2f10L3Filtered17TkIsoFiltered0p4;
   Int_t           tag_hltL3fL1sMu20L1f0L2f10QL3Filtered22Q;
   Int_t           tag_hltL3fL1sMu20L1f0L2f10QL3Filtered22QL3pfecalIsoRhoFilteredEB0p11EE0p08;
   Int_t           tag_hltL3fL1sMu20L1f0L2f10QL3Filtered22QL3pfhcalIsoRhoFilteredHB0p21HE0p22;
   Int_t           tag_hltL3fL1sMu20L1f0Tkf22QL3trkIsoFiltered0p09;
   Int_t           tag_hltL3fL1sMu20f0TkFiltered22Q;
   Int_t           tag_hltL3fL1sMu20f0TkFiltered22QL3pfecalIsoRhoFilteredEB0p11EE0p08;
   Int_t           tag_hltL3fL1sMu20f0TkFiltered22QL3pfhcalIsoRhoFilteredHB0p21HE0p22;
   Int_t           tag_hltL3fL1sMu22L1f0L2f10QL3Filtered24Q;
   Int_t           tag_hltL3fL1sMu22L1f0L2f10QL3Filtered24QL3pfecalIsoRhoFilteredEB0p11EE0p08;
   Int_t           tag_hltL3fL1sMu22L1f0L2f10QL3Filtered24QL3pfhcalIsoRhoFilteredHB0p21HE0p22;
   Int_t           tag_hltL3fL1sMu22L1f0Tkf24QL3trkIsoFiltered0p09;
   Int_t           tag_hltL3fL1sMu22f0TkFiltered24Q;
   Int_t           tag_hltL3fL1sMu22f0TkFiltered24QL3pfecalIsoRhoFilteredEB0p11EE0p08;
   Int_t           tag_hltL3fL1sMu22f0TkFiltered24QL3pfhcalIsoRhoFilteredHB0p21HE0p22;
   Float_t         pair_deltaR;
   Float_t         pair_dz;
   Float_t         pair_pt;
   Float_t         pair_rapidity;
   Float_t         pair_actualPileUp;
   Float_t         pair_genWeight;
   Float_t         pair_nJets30;
   Float_t         pair_newTuneP_mass;
   Float_t         pair_newTuneP_probe_pt;
   Float_t         pair_newTuneP_probe_sigmaPtOverPt;
   Float_t         pair_newTuneP_probe_trackType;
   Float_t         pair_probeMultiplicity;
   Float_t         pair_probeMultiplicity_Pt10_M60140;
   Float_t         pair_probeMultiplicity_TMGM;
   Float_t         pair_truePileUp;
   Int_t           pair_BestZ;

   // List of branches
   TBranch        *b_IP;   //!
   TBranch        *b_IPError;   //!
   TBranch        *b_SIP;   //!
   TBranch        *b_abseta;   //!
   TBranch        *b_caloCompatibility;   //!
   TBranch        *b_charge;   //!
   TBranch        *b_chargedHadIso03;   //!
   TBranch        *b_chargedHadIso04;   //!
   TBranch        *b_chargedParticleIso03;   //!
   TBranch        *b_chargedParticleIso04;   //!
   TBranch        *b_chi2LocMom;   //!
   TBranch        *b_chi2LocPos;   //!
   TBranch        *b_combRelIso;   //!
   TBranch        *b_combRelIsoPF03;   //!
   TBranch        *b_combRelIsoPF03dBeta;   //!
   TBranch        *b_combRelIsoPF04;   //!
   TBranch        *b_combRelIsoPF04dBeta;   //!
   TBranch        *b_dB;   //!
   TBranch        *b_ecalIso;   //!
   TBranch        *b_edB;   //!
   TBranch        *b_emEnergy;   //!
   TBranch        *b_emS9Energy;   //!
   TBranch        *b_eta;   //!
   TBranch        *b_glbChi2;   //!
   TBranch        *b_glbPtError;   //!
   TBranch        *b_glbSigmaPtOverPt;   //!
   TBranch        *b_glbTrackProb;   //!
   TBranch        *b_glbValidMuHits;   //!
   TBranch        *b_hadEnergy;   //!
   TBranch        *b_hadS9Energy;   //!
   TBranch        *b_hcalIso;   //!
   TBranch        *b_l1dphi;   //!
   TBranch        *b_l1dr;   //!
   TBranch        *b_l1drByQ;   //!
   TBranch        *b_l1eta;   //!
   TBranch        *b_l1phi;   //!
   TBranch        *b_l1pt;   //!
   TBranch        *b_l1ptByQ;   //!
   TBranch        *b_l1q;   //!
   TBranch        *b_l1qByQ;   //!
   TBranch        *b_l2dr;   //!
   TBranch        *b_l2eta;   //!
   TBranch        *b_l2pt;   //!
   TBranch        *b_l3dr;   //!
   TBranch        *b_l3pt;   //!
   TBranch        *b_neutralHadIso03;   //!
   TBranch        *b_neutralHadIso04;   //!
   TBranch        *b_numberOfMatchedStations;   //!
   TBranch        *b_numberOfMatches;   //!
   TBranch        *b_p;   //!
   TBranch        *b_phi;   //!
   TBranch        *b_photonIso03;   //!
   TBranch        *b_photonIso04;   //!
   TBranch        *b_pt;   //!
   TBranch        *b_puIso03;   //!
   TBranch        *b_puIso04;   //!
   TBranch        *b_relEcalIso;   //!
   TBranch        *b_relHcalIso;   //!
   TBranch        *b_relTkIso;   //!
   TBranch        *b_segmentCompatibility;   //!
   TBranch        *b_staQoverP;   //!
   TBranch        *b_staQoverPerror;   //!
   TBranch        *b_staValidStations;   //!
   TBranch        *b_tkChi2;   //!
   TBranch        *b_tkExpHitIn;   //!
   TBranch        *b_tkExpHitOut;   //!
   TBranch        *b_tkHitFract;   //!
   TBranch        *b_tkIso;   //!
   TBranch        *b_tkKink;   //!
   TBranch        *b_tkPixelLay;   //!
   TBranch        *b_tkPtError;   //!
   TBranch        *b_tkSigmaPtOverPt;   //!
   TBranch        *b_tkTrackerLay;   //!
   TBranch        *b_tkValidHits;   //!
   TBranch        *b_tkValidPixelHits;   //!
   TBranch        *b_JetBTagCSV;   //!
   TBranch        *b_JetNDauCharged;   //!
   TBranch        *b_JetPtRatio;   //!
   TBranch        *b_JetPtRel;   //!
   TBranch        *b_activity_miniIsoCharged;   //!
   TBranch        *b_activity_miniIsoNeutrals;   //!
   TBranch        *b_activity_miniIsoPUCharged;   //!
   TBranch        *b_activity_miniIsoPhotons;   //!
   TBranch        *b_dxyBS;   //!
   TBranch        *b_dxyPVdzmin;   //!
   TBranch        *b_dzPV;   //!
   TBranch        *b_fixedGridRhoFastjetAll;   //!
   TBranch        *b_fixedGridRhoFastjetAllCalo;   //!
   TBranch        *b_fixedGridRhoFastjetCentralCalo;   //!
   TBranch        *b_fixedGridRhoFastjetCentralChargedPileUp;   //!
   TBranch        *b_fixedGridRhoFastjetCentralNeutral;   //!
   TBranch        *b_isoTrk03Abs;   //!
   TBranch        *b_isoTrk03Rel;   //!
   TBranch        *b_miniIsoCharged;   //!
   TBranch        *b_miniIsoNeutrals;   //!
   TBranch        *b_miniIsoPUCharged;   //!
   TBranch        *b_miniIsoPhotons;   //!
   TBranch        *b_mt;   //!
   TBranch        *b_muPFIsoValueCHR04PUPPI;   //!
   TBranch        *b_muPFIsoValueCHR04PUPPINoLep;   //!
   TBranch        *b_muPFIsoValueNHR04PUPPI;   //!
   TBranch        *b_muPFIsoValueNHR04PUPPINoLep;   //!
   TBranch        *b_muPFIsoValuePhR04PUPPI;   //!
   TBranch        *b_muPFIsoValuePhR04PUPPINoLep;   //!
   TBranch        *b_nSplitTk;   //!
   TBranch        *b_Calo;   //!
   TBranch        *b_DiMuonGlb17Glb8RelTrkIsoFiltered0p4;   //!
   TBranch        *b_DoubleIsoMu17Mu8_IsoMu17leg;   //!
   TBranch        *b_DoubleIsoMu17Mu8_IsoMu8leg;   //!
   TBranch        *b_DoubleIsoMu17Mu8_Mu17leg;   //!
   TBranch        *b_DoubleIsoMu17Mu8_Mu8leg;   //!
   TBranch        *b_DoubleIsoMu17Mu8dZ_Mu17leg;   //!
   TBranch        *b_DoubleIsoMu17TkMu8_IsoMu17leg;   //!
   TBranch        *b_DoubleIsoMu17TkMu8_IsoMu8leg;   //!
   TBranch        *b_DoubleIsoMu17TkMu8_Mu17leg;   //!
   TBranch        *b_DoubleIsoMu17TkMu8_TkMu8leg;   //!
   TBranch        *b_DoubleIsoMu17TkMu8dZ_Mu17;   //!
   TBranch        *b_DoubleMu30TkMu11;   //!
   TBranch        *b_DoubleMu30TkMu11_Mu30leg;   //!
   TBranch        *b_DoubleMu30TkMu11_TkMu11leg;   //!
   TBranch        *b_Glb;   //!
   TBranch        *b_GlbPT;   //!
   TBranch        *b_HLT_TkMu50;   //!
   TBranch        *b_HWWID;   //!
   TBranch        *b_HighPt;   //!
   TBranch        *b_IsoMu18;   //!
   TBranch        *b_IsoMu20;   //!
   TBranch        *b_IsoMu22;   //!
   TBranch        *b_IsoMu22_eta2p1;   //!
   TBranch        *b_IsoMu24;   //!
   TBranch        *b_IsoMu24_eta2p1;   //!
   TBranch        *b_IsoMu27;   //!
   TBranch        *b_IsoTkMu18;   //!
   TBranch        *b_IsoTkMu20;   //!
   TBranch        *b_IsoTkMu22;   //!
   TBranch        *b_IsoTkMu22_eta2p1;   //!
   TBranch        *b_IsoTkMu24;   //!
   TBranch        *b_IsoTkMu24_eta2p1;   //!
   TBranch        *b_IsoTkMu27;   //!
   TBranch        *b_L1sMu16;   //!
   TBranch        *b_L1sMu18;   //!
   TBranch        *b_L1sMu20;   //!
   TBranch        *b_L1sMu25Eta2p1;   //!
   TBranch        *b_L2fL1sDoubleMu114L1f0L2Filtered10OneMu;   //!
   TBranch        *b_L2fL1sDoubleMu114L1f0OneMuL2Filtered10;   //!
   TBranch        *b_L2fL1sMu18L1f0L2Filtered10Q;   //!
   TBranch        *b_L2pfL1sDoubleMu114L1f0L2PreFiltered0;   //!
   TBranch        *b_L3crIsoL1sMu16L1f0L2f10QL3f18QL3trkIsoFiltered0p09;   //!
   TBranch        *b_L3fL1sDoubleMu114L1f0L2f10L3Filtered17;   //!
   TBranch        *b_L3fL1sDoubleMu114L1f0L2f10OneMuL3Filtered17;   //!
   TBranch        *b_L3fL1sMu16L1f0L2f10QL3Filtered18Q;   //!
   TBranch        *b_L3fL1sMu16L1f0L2f10QL3Filtered18QL3pfecalIsoRhoFilteredEB0p11EE0p08;   //!
   TBranch        *b_L3fL1sMu16L1f0L2f10QL3Filtered18QL3pfhcalIsoRhoFilteredHB0p21HE0p22;   //!
   TBranch        *b_L3fL1sMu16f0TkFiltered18Q;   //!
   TBranch        *b_L3pfL1sDoubleMu114L1f0L2pf0L3PreFiltered8;   //!
   TBranch        *b_Loose;   //!
   TBranch        *b_Medium;   //!
   TBranch        *b_Medium2016;   //!
   TBranch        *b_Mu17;   //!
   TBranch        *b_Mu17_IsoTrkVVL;   //!
   TBranch        *b_Mu20;   //!
   TBranch        *b_Mu23_TrkIsoVVL;   //!
   TBranch        *b_Mu45_eta2p1;   //!
   TBranch        *b_Mu50;   //!
   TBranch        *b_MuIDForOutsideInTk;   //!
   TBranch        *b_PF;   //!
   TBranch        *b_TM;   //!
   TBranch        *b_TMA;   //!
   TBranch        *b_TMLSAT;   //!
   TBranch        *b_TMLST;   //!
   TBranch        *b_TMOSL;   //!
   TBranch        *b_TMOST;   //!
   TBranch        *b_TMOSTQual;   //!
   TBranch        *b_Tight2012;   //!
   TBranch        *b_Track_HP;   //!
   TBranch        *b_VBTF;   //!
   TBranch        *b_VBTF_nL8;   //!
   TBranch        *b_VBTF_nL9;   //!
   TBranch        *b_hltL2fL1sMu10lqL1f0L2Filtered10;   //!
   TBranch        *b_hltL2fL1sMu20L1f0L2Filtered10Q;   //!
   TBranch        *b_hltL2fL1sMu22L1f0L2Filtered10Q;   //!
   TBranch        *b_hltL3crIsoL1sMu18L1f0L2f10QL3f20QL3trkIsoFiltered0p09;   //!
   TBranch        *b_hltL3crIsoL1sMu20L1f0L2f10QL3f22QL3trkIsoFiltered0p09;   //!
   TBranch        *b_hltL3crIsoL1sMu22L1f0L2f10QL3f24QL3trkIsoFiltered0p09;   //!
   TBranch        *b_hltL3fL1sMu10lqL1f0L2f10L3Filtered17;   //!
   TBranch        *b_hltL3fL1sMu18L1f0L2f10QL3Filtered20Q;   //!
   TBranch        *b_hltL3fL1sMu18L1f0L2f10QL3Filtered20QL3pfecalIsoRhoFilteredEB0p11EE0p08;   //!
   TBranch        *b_hltL3fL1sMu18L1f0L2f10QL3Filtered20QL3pfhcalIsoRhoFilteredHB0p21HE0p22;   //!
   TBranch        *b_hltL3fL1sMu18L1f0Tkf20QL3trkIsoFiltered0p09;   //!
   TBranch        *b_hltL3fL1sMu18f0TkFiltered20Q;   //!
   TBranch        *b_hltL3fL1sMu18f0TkFiltered20QL3pfecalIsoRhoFilteredEB0p11EE0p08;   //!
   TBranch        *b_hltL3fL1sMu18f0TkFiltered20QL3pfhcalIsoRhoFilteredHB0p21HE0p22;   //!
   TBranch        *b_hltL3fL1sMu1lqL1f0L2f10L3Filtered17TkIsoFiltered0p4;   //!
   TBranch        *b_hltL3fL1sMu20L1f0L2f10QL3Filtered22Q;   //!
   TBranch        *b_hltL3fL1sMu20L1f0L2f10QL3Filtered22QL3pfecalIsoRhoFilteredEB0p11EE0p08;   //!
   TBranch        *b_hltL3fL1sMu20L1f0L2f10QL3Filtered22QL3pfhcalIsoRhoFilteredHB0p21HE0p22;   //!
   TBranch        *b_hltL3fL1sMu20L1f0Tkf22QL3trkIsoFiltered0p09;   //!
   TBranch        *b_hltL3fL1sMu20f0TkFiltered22Q;   //!
   TBranch        *b_hltL3fL1sMu20f0TkFiltered22QL3pfecalIsoRhoFilteredEB0p11EE0p08;   //!
   TBranch        *b_hltL3fL1sMu20f0TkFiltered22QL3pfhcalIsoRhoFilteredHB0p21HE0p22;   //!
   TBranch        *b_hltL3fL1sMu22L1f0L2f10QL3Filtered24Q;   //!
   TBranch        *b_hltL3fL1sMu22L1f0L2f10QL3Filtered24QL3pfecalIsoRhoFilteredEB0p11EE0p08;   //!
   TBranch        *b_hltL3fL1sMu22L1f0L2f10QL3Filtered24QL3pfhcalIsoRhoFilteredHB0p21HE0p22;   //!
   TBranch        *b_hltL3fL1sMu22L1f0Tkf24QL3trkIsoFiltered0p09;   //!
   TBranch        *b_hltL3fL1sMu22f0TkFiltered24Q;   //!
   TBranch        *b_hltL3fL1sMu22f0TkFiltered24QL3pfecalIsoRhoFilteredEB0p11EE0p08;   //!
   TBranch        *b_hltL3fL1sMu22f0TkFiltered24QL3pfhcalIsoRhoFilteredHB0p21HE0p22;   //!
   TBranch        *b_tkHighPt;   //!
   TBranch        *b_run;   //!
   TBranch        *b_lumi;   //!
   TBranch        *b_event;   //!
   TBranch        *b_truePU;   //!
   TBranch        *b_mass;   //!
   TBranch        *b_mcTrue;   //!
   TBranch        *b_mcMass;   //!
   TBranch        *b_tag_IP;   //!
   TBranch        *b_tag_IPError;   //!
   TBranch        *b_tag_SIP;   //!
   TBranch        *b_tag_abseta;   //!
   TBranch        *b_tag_caloCompatibility;   //!
   TBranch        *b_tag_charge;   //!
   TBranch        *b_tag_chargedHadIso03;   //!
   TBranch        *b_tag_chargedHadIso04;   //!
   TBranch        *b_tag_chargedParticleIso03;   //!
   TBranch        *b_tag_chargedParticleIso04;   //!
   TBranch        *b_tag_chi2LocMom;   //!
   TBranch        *b_tag_chi2LocPos;   //!
   TBranch        *b_tag_combRelIso;   //!
   TBranch        *b_tag_combRelIsoPF03;   //!
   TBranch        *b_tag_combRelIsoPF03dBeta;   //!
   TBranch        *b_tag_combRelIsoPF04;   //!
   TBranch        *b_tag_combRelIsoPF04dBeta;   //!
   TBranch        *b_tag_dB;   //!
   TBranch        *b_tag_ecalIso;   //!
   TBranch        *b_tag_edB;   //!
   TBranch        *b_tag_emEnergy;   //!
   TBranch        *b_tag_emS9Energy;   //!
   TBranch        *b_tag_eta;   //!
   TBranch        *b_tag_glbChi2;   //!
   TBranch        *b_tag_glbPtError;   //!
   TBranch        *b_tag_glbSigmaPtOverPt;   //!
   TBranch        *b_tag_glbTrackProb;   //!
   TBranch        *b_tag_glbValidMuHits;   //!
   TBranch        *b_tag_hadEnergy;   //!
   TBranch        *b_tag_hadS9Energy;   //!
   TBranch        *b_tag_hcalIso;   //!
   TBranch        *b_tag_l1dphi;   //!
   TBranch        *b_tag_l1dr;   //!
   TBranch        *b_tag_l1drByQ;   //!
   TBranch        *b_tag_l1eta;   //!
   TBranch        *b_tag_l1phi;   //!
   TBranch        *b_tag_l1pt;   //!
   TBranch        *b_tag_l1ptByQ;   //!
   TBranch        *b_tag_l1q;   //!
   TBranch        *b_tag_l1qByQ;   //!
   TBranch        *b_tag_l2dr;   //!
   TBranch        *b_tag_l2eta;   //!
   TBranch        *b_tag_l2pt;   //!
   TBranch        *b_tag_l3dr;   //!
   TBranch        *b_tag_l3pt;   //!
   TBranch        *b_tag_neutralHadIso03;   //!
   TBranch        *b_tag_neutralHadIso04;   //!
   TBranch        *b_tag_numberOfMatchedStations;   //!
   TBranch        *b_tag_numberOfMatches;   //!
   TBranch        *b_tag_p;   //!
   TBranch        *b_tag_phi;   //!
   TBranch        *b_tag_photonIso03;   //!
   TBranch        *b_tag_photonIso04;   //!
   TBranch        *b_tag_pt;   //!
   TBranch        *b_tag_puIso03;   //!
   TBranch        *b_tag_puIso04;   //!
   TBranch        *b_tag_relEcalIso;   //!
   TBranch        *b_tag_relHcalIso;   //!
   TBranch        *b_tag_relTkIso;   //!
   TBranch        *b_tag_segmentCompatibility;   //!
   TBranch        *b_tag_staQoverP;   //!
   TBranch        *b_tag_staQoverPerror;   //!
   TBranch        *b_tag_staValidStations;   //!
   TBranch        *b_tag_tkChi2;   //!
   TBranch        *b_tag_tkExpHitIn;   //!
   TBranch        *b_tag_tkExpHitOut;   //!
   TBranch        *b_tag_tkHitFract;   //!
   TBranch        *b_tag_tkIso;   //!
   TBranch        *b_tag_tkKink;   //!
   TBranch        *b_tag_tkPixelLay;   //!
   TBranch        *b_tag_tkPtError;   //!
   TBranch        *b_tag_tkSigmaPtOverPt;   //!
   TBranch        *b_tag_tkTrackerLay;   //!
   TBranch        *b_tag_tkValidHits;   //!
   TBranch        *b_tag_tkValidPixelHits;   //!
   TBranch        *b_tag_dxyBS;   //!
   TBranch        *b_tag_dxyPVdzmin;   //!
   TBranch        *b_tag_dzPV;   //!
   TBranch        *b_tag_fixedGridRhoFastjetAll;   //!
   TBranch        *b_tag_fixedGridRhoFastjetAllCalo;   //!
   TBranch        *b_tag_fixedGridRhoFastjetCentralCalo;   //!
   TBranch        *b_tag_fixedGridRhoFastjetCentralChargedPileUp;   //!
   TBranch        *b_tag_fixedGridRhoFastjetCentralNeutral;   //!
   TBranch        *b_tag_isoTrk03Abs;   //!
   TBranch        *b_tag_isoTrk03Rel;   //!
   TBranch        *b_tag_met;   //!
   TBranch        *b_tag_mt;   //!
   TBranch        *b_tag_nSplitTk;   //!
   TBranch        *b_tag_nVertices;   //!
   TBranch        *b_tag_DiMuonGlb17Glb8RelTrkIsoFiltered0p4;   //!
   TBranch        *b_tag_DoubleIsoMu17Mu8_IsoMu17leg;   //!
   TBranch        *b_tag_DoubleIsoMu17Mu8_IsoMu8leg;   //!
   TBranch        *b_tag_DoubleIsoMu17Mu8_Mu17leg;   //!
   TBranch        *b_tag_DoubleIsoMu17Mu8_Mu8leg;   //!
   TBranch        *b_tag_DoubleIsoMu17Mu8dZ_Mu17leg;   //!
   TBranch        *b_tag_DoubleIsoMu17TkMu8_IsoMu17leg;   //!
   TBranch        *b_tag_DoubleIsoMu17TkMu8_IsoMu8leg;   //!
   TBranch        *b_tag_DoubleIsoMu17TkMu8_Mu17leg;   //!
   TBranch        *b_tag_DoubleIsoMu17TkMu8_TkMu8leg;   //!
   TBranch        *b_tag_DoubleIsoMu17TkMu8dZ_Mu17;   //!
   TBranch        *b_tag_DoubleMu30TkMu11;   //!
   TBranch        *b_tag_DoubleMu30TkMu11_Mu30leg;   //!
   TBranch        *b_tag_DoubleMu30TkMu11_TkMu11leg;   //!
   TBranch        *b_tag_HLT_TkMu50;   //!
   TBranch        *b_tag_IsoMu18;   //!
   TBranch        *b_tag_IsoMu20;   //!
   TBranch        *b_tag_IsoMu22;   //!
   TBranch        *b_tag_IsoMu22_eta2p1;   //!
   TBranch        *b_tag_IsoMu24;   //!
   TBranch        *b_tag_IsoMu24_eta2p1;   //!
   TBranch        *b_tag_IsoMu27;   //!
   TBranch        *b_tag_IsoTkMu18;   //!
   TBranch        *b_tag_IsoTkMu20;   //!
   TBranch        *b_tag_IsoTkMu22;   //!
   TBranch        *b_tag_IsoTkMu22_eta2p1;   //!
   TBranch        *b_tag_IsoTkMu24;   //!
   TBranch        *b_tag_IsoTkMu24_eta2p1;   //!
   TBranch        *b_tag_IsoTkMu27;   //!
   TBranch        *b_tag_L1sMu16;   //!
   TBranch        *b_tag_L1sMu18;   //!
   TBranch        *b_tag_L1sMu20;   //!
   TBranch        *b_tag_L1sMu25Eta2p1;   //!
   TBranch        *b_tag_L2fL1sDoubleMu114L1f0L2Filtered10OneMu;   //!
   TBranch        *b_tag_L2fL1sDoubleMu114L1f0OneMuL2Filtered10;   //!
   TBranch        *b_tag_L2fL1sMu18L1f0L2Filtered10Q;   //!
   TBranch        *b_tag_L2pfL1sDoubleMu114L1f0L2PreFiltered0;   //!
   TBranch        *b_tag_L3crIsoL1sMu16L1f0L2f10QL3f18QL3trkIsoFiltered0p09;   //!
   TBranch        *b_tag_L3fL1sDoubleMu114L1f0L2f10L3Filtered17;   //!
   TBranch        *b_tag_L3fL1sDoubleMu114L1f0L2f10OneMuL3Filtered17;   //!
   TBranch        *b_tag_L3fL1sMu16L1f0L2f10QL3Filtered18Q;   //!
   TBranch        *b_tag_L3fL1sMu16L1f0L2f10QL3Filtered18QL3pfecalIsoRhoFilteredEB0p11EE0p08;   //!
   TBranch        *b_tag_L3fL1sMu16L1f0L2f10QL3Filtered18QL3pfhcalIsoRhoFilteredHB0p21HE0p22;   //!
   TBranch        *b_tag_L3fL1sMu16f0TkFiltered18Q;   //!
   TBranch        *b_tag_L3pfL1sDoubleMu114L1f0L2pf0L3PreFiltered8;   //!
   TBranch        *b_tag_Mu17;   //!
   TBranch        *b_tag_Mu17_IsoTrkVVL;   //!
   TBranch        *b_tag_Mu20;   //!
   TBranch        *b_tag_Mu23_TrkIsoVVL;   //!
   TBranch        *b_tag_Mu45_eta2p1;   //!
   TBranch        *b_tag_Mu50;   //!
   TBranch        *b_tag_hltL2fL1sMu10lqL1f0L2Filtered10;   //!
   TBranch        *b_tag_hltL2fL1sMu20L1f0L2Filtered10Q;   //!
   TBranch        *b_tag_hltL2fL1sMu22L1f0L2Filtered10Q;   //!
   TBranch        *b_tag_hltL3crIsoL1sMu18L1f0L2f10QL3f20QL3trkIsoFiltered0p09;   //!
   TBranch        *b_tag_hltL3crIsoL1sMu20L1f0L2f10QL3f22QL3trkIsoFiltered0p09;   //!
   TBranch        *b_tag_hltL3crIsoL1sMu22L1f0L2f10QL3f24QL3trkIsoFiltered0p09;   //!
   TBranch        *b_tag_hltL3fL1sMu10lqL1f0L2f10L3Filtered17;   //!
   TBranch        *b_tag_hltL3fL1sMu18L1f0L2f10QL3Filtered20Q;   //!
   TBranch        *b_tag_hltL3fL1sMu18L1f0L2f10QL3Filtered20QL3pfecalIsoRhoFilteredEB0p11EE0p08;   //!
   TBranch        *b_tag_hltL3fL1sMu18L1f0L2f10QL3Filtered20QL3pfhcalIsoRhoFilteredHB0p21HE0p22;   //!
   TBranch        *b_tag_hltL3fL1sMu18L1f0Tkf20QL3trkIsoFiltered0p09;   //!
   TBranch        *b_tag_hltL3fL1sMu18f0TkFiltered20Q;   //!
   TBranch        *b_tag_hltL3fL1sMu18f0TkFiltered20QL3pfecalIsoRhoFilteredEB0p11EE0p08;   //!
   TBranch        *b_tag_hltL3fL1sMu18f0TkFiltered20QL3pfhcalIsoRhoFilteredHB0p21HE0p22;   //!
   TBranch        *b_tag_hltL3fL1sMu1lqL1f0L2f10L3Filtered17TkIsoFiltered0p4;   //!
   TBranch        *b_tag_hltL3fL1sMu20L1f0L2f10QL3Filtered22Q;   //!
   TBranch        *b_tag_hltL3fL1sMu20L1f0L2f10QL3Filtered22QL3pfecalIsoRhoFilteredEB0p11EE0p08;   //!
   TBranch        *b_tag_hltL3fL1sMu20L1f0L2f10QL3Filtered22QL3pfhcalIsoRhoFilteredHB0p21HE0p22;   //!
   TBranch        *b_tag_hltL3fL1sMu20L1f0Tkf22QL3trkIsoFiltered0p09;   //!
   TBranch        *b_tag_hltL3fL1sMu20f0TkFiltered22Q;   //!
   TBranch        *b_tag_hltL3fL1sMu20f0TkFiltered22QL3pfecalIsoRhoFilteredEB0p11EE0p08;   //!
   TBranch        *b_tag_hltL3fL1sMu20f0TkFiltered22QL3pfhcalIsoRhoFilteredHB0p21HE0p22;   //!
   TBranch        *b_tag_hltL3fL1sMu22L1f0L2f10QL3Filtered24Q;   //!
   TBranch        *b_tag_hltL3fL1sMu22L1f0L2f10QL3Filtered24QL3pfecalIsoRhoFilteredEB0p11EE0p08;   //!
   TBranch        *b_tag_hltL3fL1sMu22L1f0L2f10QL3Filtered24QL3pfhcalIsoRhoFilteredHB0p21HE0p22;   //!
   TBranch        *b_tag_hltL3fL1sMu22L1f0Tkf24QL3trkIsoFiltered0p09;   //!
   TBranch        *b_tag_hltL3fL1sMu22f0TkFiltered24Q;   //!
   TBranch        *b_tag_hltL3fL1sMu22f0TkFiltered24QL3pfecalIsoRhoFilteredEB0p11EE0p08;   //!
   TBranch        *b_tag_hltL3fL1sMu22f0TkFiltered24QL3pfhcalIsoRhoFilteredHB0p21HE0p22;   //!
   TBranch        *b_pair_deltaR;   //!
   TBranch        *b_pair_dz;   //!
   TBranch        *b_pair_pt;   //!
   TBranch        *b_pair_rapidity;   //!
   TBranch        *b_pair_actualPileUp;   //!
   TBranch        *b_pair_genWeight;   //!
   TBranch        *b_pair_nJets30;   //!
   TBranch        *b_pair_newTuneP_mass;   //!
   TBranch        *b_pair_newTuneP_probe_pt;   //!
   TBranch        *b_pair_newTuneP_probe_sigmaPtOverPt;   //!
   TBranch        *b_pair_newTuneP_probe_trackType;   //!
   TBranch        *b_pair_probeMultiplicity;   //!
   TBranch        *b_pair_probeMultiplicity_Pt10_M60140;   //!
   TBranch        *b_pair_probeMultiplicity_TMGM;   //!
   TBranch        *b_pair_truePileUp;   //!
   TBranch        *b_pair_BestZ;   //!

   MuonPogFitterTree(TTree *tree, bool isMC=false);
   virtual ~MuonPogFitterTree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree, bool isMC=false);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef MuonPogFitterTree_cxx
MuonPogFitterTree::MuonPogFitterTree(TTree *tree, bool isMC) : fChain(0) 
{
   assert(tree);
   Init(tree, isMC);
}

MuonPogFitterTree::~MuonPogFitterTree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t MuonPogFitterTree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t MuonPogFitterTree::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void MuonPogFitterTree::Init(TTree *tree, bool isMC)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   //fChain->SetBranchAddress("IP", &IP, &b_IP);
   //fChain->SetBranchAddress("IPError", &IPError, &b_IPError);
   //fChain->SetBranchAddress("SIP", &SIP, &b_SIP);
   fChain->SetBranchAddress("abseta", &abseta, &b_abseta);
   //fChain->SetBranchAddress("caloCompatibility", &caloCompatibility, &b_caloCompatibility);
   fChain->SetBranchAddress("charge", &charge, &b_charge);
   //fChain->SetBranchAddress("chargedHadIso03", &chargedHadIso03, &b_chargedHadIso03);
   //fChain->SetBranchAddress("chargedHadIso04", &chargedHadIso04, &b_chargedHadIso04);
   //fChain->SetBranchAddress("chargedParticleIso03", &chargedParticleIso03, &b_chargedParticleIso03);
   //fChain->SetBranchAddress("chargedParticleIso04", &chargedParticleIso04, &b_chargedParticleIso04);
   //fChain->SetBranchAddress("chi2LocMom", &chi2LocMom, &b_chi2LocMom);
   //fChain->SetBranchAddress("chi2LocPos", &chi2LocPos, &b_chi2LocPos);
   //fChain->SetBranchAddress("combRelIso", &combRelIso, &b_combRelIso);
   //fChain->SetBranchAddress("combRelIsoPF03", &combRelIsoPF03, &b_combRelIsoPF03);
   //fChain->SetBranchAddress("combRelIsoPF03dBeta", &combRelIsoPF03dBeta, &b_combRelIsoPF03dBeta);
   //fChain->SetBranchAddress("combRelIsoPF04", &combRelIsoPF04, &b_combRelIsoPF04);
   fChain->SetBranchAddress("combRelIsoPF04dBeta", &combRelIsoPF04dBeta, &b_combRelIsoPF04dBeta);
   //fChain->SetBranchAddress("dB", &dB, &b_dB);
   //fChain->SetBranchAddress("ecalIso", &ecalIso, &b_ecalIso);
   //fChain->SetBranchAddress("edB", &edB, &b_edB);
   //fChain->SetBranchAddress("emEnergy", &emEnergy, &b_emEnergy);
   //fChain->SetBranchAddress("emS9Energy", &emS9Energy, &b_emS9Energy);
   fChain->SetBranchAddress("eta", &eta, &b_eta);
   //fChain->SetBranchAddress("glbChi2", &glbChi2, &b_glbChi2);
   //fChain->SetBranchAddress("glbPtError", &glbPtError, &b_glbPtError);
   //fChain->SetBranchAddress("glbSigmaPtOverPt", &glbSigmaPtOverPt, &b_glbSigmaPtOverPt);
   //fChain->SetBranchAddress("glbTrackProb", &glbTrackProb, &b_glbTrackProb);
   //fChain->SetBranchAddress("glbValidMuHits", &glbValidMuHits, &b_glbValidMuHits);
   //fChain->SetBranchAddress("hadEnergy", &hadEnergy, &b_hadEnergy);
   //fChain->SetBranchAddress("hadS9Energy", &hadS9Energy, &b_hadS9Energy);
   //fChain->SetBranchAddress("hcalIso", &hcalIso, &b_hcalIso);
   //fChain->SetBranchAddress("l1dphi", &l1dphi, &b_l1dphi);
   //fChain->SetBranchAddress("l1dr", &l1dr, &b_l1dr);
   //fChain->SetBranchAddress("l1drByQ", &l1drByQ, &b_l1drByQ);
   //fChain->SetBranchAddress("l1eta", &l1eta, &b_l1eta);
   //fChain->SetBranchAddress("l1phi", &l1phi, &b_l1phi);
   //fChain->SetBranchAddress("l1pt", &l1pt, &b_l1pt);
   //fChain->SetBranchAddress("l1ptByQ", &l1ptByQ, &b_l1ptByQ);
   //fChain->SetBranchAddress("l1q", &l1q, &b_l1q);
   //fChain->SetBranchAddress("l1qByQ", &l1qByQ, &b_l1qByQ);
   //fChain->SetBranchAddress("l2dr", &l2dr, &b_l2dr);
   //fChain->SetBranchAddress("l2eta", &l2eta, &b_l2eta);
   //fChain->SetBranchAddress("l2pt", &l2pt, &b_l2pt);
   //fChain->SetBranchAddress("l3dr", &l3dr, &b_l3dr);
   //fChain->SetBranchAddress("l3pt", &l3pt, &b_l3pt);
   //fChain->SetBranchAddress("neutralHadIso03", &neutralHadIso03, &b_neutralHadIso03);
   //fChain->SetBranchAddress("neutralHadIso04", &neutralHadIso04, &b_neutralHadIso04);
   //fChain->SetBranchAddress("numberOfMatchedStations", &numberOfMatchedStations, &b_numberOfMatchedStations);
   //fChain->SetBranchAddress("numberOfMatches", &numberOfMatches, &b_numberOfMatches);
   fChain->SetBranchAddress("p", &p, &b_p);
   fChain->SetBranchAddress("phi", &phi, &b_phi);
   //fChain->SetBranchAddress("photonIso03", &photonIso03, &b_photonIso03);
   //fChain->SetBranchAddress("photonIso04", &photonIso04, &b_photonIso04);
   fChain->SetBranchAddress("pt", &pt, &b_pt);
   //fChain->SetBranchAddress("puIso03", &puIso03, &b_puIso03);
   //fChain->SetBranchAddress("puIso04", &puIso04, &b_puIso04);
   //fChain->SetBranchAddress("relEcalIso", &relEcalIso, &b_relEcalIso);
   //fChain->SetBranchAddress("relHcalIso", &relHcalIso, &b_relHcalIso);
   //fChain->SetBranchAddress("relTkIso", &relTkIso, &b_relTkIso);
   //fChain->SetBranchAddress("segmentCompatibility", &segmentCompatibility, &b_segmentCompatibility);
   //fChain->SetBranchAddress("staQoverP", &staQoverP, &b_staQoverP);
   //fChain->SetBranchAddress("staQoverPerror", &staQoverPerror, &b_staQoverPerror);
   //fChain->SetBranchAddress("staValidStations", &staValidStations, &b_staValidStations);
   //fChain->SetBranchAddress("tkChi2", &tkChi2, &b_tkChi2);
   //fChain->SetBranchAddress("tkExpHitIn", &tkExpHitIn, &b_tkExpHitIn);
   //fChain->SetBranchAddress("tkExpHitOut", &tkExpHitOut, &b_tkExpHitOut);
   //fChain->SetBranchAddress("tkHitFract", &tkHitFract, &b_tkHitFract);
   //fChain->SetBranchAddress("tkIso", &tkIso, &b_tkIso);
   //fChain->SetBranchAddress("tkKink", &tkKink, &b_tkKink);
   //fChain->SetBranchAddress("tkPixelLay", &tkPixelLay, &b_tkPixelLay);
   //fChain->SetBranchAddress("tkPtError", &tkPtError, &b_tkPtError);
   //fChain->SetBranchAddress("tkSigmaPtOverPt", &tkSigmaPtOverPt, &b_tkSigmaPtOverPt);
   //fChain->SetBranchAddress("tkTrackerLay", &tkTrackerLay, &b_tkTrackerLay);
   //fChain->SetBranchAddress("tkValidHits", &tkValidHits, &b_tkValidHits);
   //fChain->SetBranchAddress("tkValidPixelHits", &tkValidPixelHits, &b_tkValidPixelHits);
   //fChain->SetBranchAddress("JetBTagCSV", &JetBTagCSV, &b_JetBTagCSV);
   //fChain->SetBranchAddress("JetNDauCharged", &JetNDauCharged, &b_JetNDauCharged);
   //fChain->SetBranchAddress("JetPtRatio", &JetPtRatio, &b_JetPtRatio);
   //fChain->SetBranchAddress("JetPtRel", &JetPtRel, &b_JetPtRel);
   //fChain->SetBranchAddress("activity_miniIsoCharged", &activity_miniIsoCharged, &b_activity_miniIsoCharged);
   //fChain->SetBranchAddress("activity_miniIsoNeutrals", &activity_miniIsoNeutrals, &b_activity_miniIsoNeutrals);
   //fChain->SetBranchAddress("activity_miniIsoPUCharged", &activity_miniIsoPUCharged, &b_activity_miniIsoPUCharged);
   //fChain->SetBranchAddress("activity_miniIsoPhotons", &activity_miniIsoPhotons, &b_activity_miniIsoPhotons);
   //fChain->SetBranchAddress("dxyBS", &dxyBS, &b_dxyBS);
   //fChain->SetBranchAddress("dxyPVdzmin", &dxyPVdzmin, &b_dxyPVdzmin);
   //fChain->SetBranchAddress("dzPV", &dzPV, &b_dzPV);
   //fChain->SetBranchAddress("fixedGridRhoFastjetAll", &fixedGridRhoFastjetAll, &b_fixedGridRhoFastjetAll);
   //fChain->SetBranchAddress("fixedGridRhoFastjetAllCalo", &fixedGridRhoFastjetAllCalo, &b_fixedGridRhoFastjetAllCalo);
   //fChain->SetBranchAddress("fixedGridRhoFastjetCentralCalo", &fixedGridRhoFastjetCentralCalo, &b_fixedGridRhoFastjetCentralCalo);
   //fChain->SetBranchAddress("fixedGridRhoFastjetCentralChargedPileUp", &fixedGridRhoFastjetCentralChargedPileUp, &b_fixedGridRhoFastjetCentralChargedPileUp);
   //fChain->SetBranchAddress("fixedGridRhoFastjetCentralNeutral", &fixedGridRhoFastjetCentralNeutral, &b_fixedGridRhoFastjetCentralNeutral);
   //fChain->SetBranchAddress("isoTrk03Abs", &isoTrk03Abs, &b_isoTrk03Abs);
   //fChain->SetBranchAddress("isoTrk03Rel", &isoTrk03Rel, &b_isoTrk03Rel);
   //fChain->SetBranchAddress("miniIsoCharged", &miniIsoCharged, &b_miniIsoCharged);
   //fChain->SetBranchAddress("miniIsoNeutrals", &miniIsoNeutrals, &b_miniIsoNeutrals);
   //fChain->SetBranchAddress("miniIsoPUCharged", &miniIsoPUCharged, &b_miniIsoPUCharged);
   //fChain->SetBranchAddress("miniIsoPhotons", &miniIsoPhotons, &b_miniIsoPhotons);
   //fChain->SetBranchAddress("mt", &mt, &b_mt);
   //fChain->SetBranchAddress("muPFIsoValueCHR04PUPPI", &muPFIsoValueCHR04PUPPI, &b_muPFIsoValueCHR04PUPPI);
   //fChain->SetBranchAddress("muPFIsoValueCHR04PUPPINoLep", &muPFIsoValueCHR04PUPPINoLep, &b_muPFIsoValueCHR04PUPPINoLep);
   //fChain->SetBranchAddress("muPFIsoValueNHR04PUPPI", &muPFIsoValueNHR04PUPPI, &b_muPFIsoValueNHR04PUPPI);
   //fChain->SetBranchAddress("muPFIsoValueNHR04PUPPINoLep", &muPFIsoValueNHR04PUPPINoLep, &b_muPFIsoValueNHR04PUPPINoLep);
   //fChain->SetBranchAddress("muPFIsoValuePhR04PUPPI", &muPFIsoValuePhR04PUPPI, &b_muPFIsoValuePhR04PUPPI);
   //fChain->SetBranchAddress("muPFIsoValuePhR04PUPPINoLep", &muPFIsoValuePhR04PUPPINoLep, &b_muPFIsoValuePhR04PUPPINoLep);
   //fChain->SetBranchAddress("nSplitTk", &nSplitTk, &b_nSplitTk);
   //fChain->SetBranchAddress("Calo", &Calo, &b_Calo);
   //fChain->SetBranchAddress("DiMuonGlb17Glb8RelTrkIsoFiltered0p4", &DiMuonGlb17Glb8RelTrkIsoFiltered0p4, &b_DiMuonGlb17Glb8RelTrkIsoFiltered0p4);
   //fChain->SetBranchAddress("DoubleIsoMu17Mu8_IsoMu17leg", &DoubleIsoMu17Mu8_IsoMu17leg, &b_DoubleIsoMu17Mu8_IsoMu17leg);
   //fChain->SetBranchAddress("DoubleIsoMu17Mu8_IsoMu8leg", &DoubleIsoMu17Mu8_IsoMu8leg, &b_DoubleIsoMu17Mu8_IsoMu8leg);
   //fChain->SetBranchAddress("DoubleIsoMu17Mu8_Mu17leg", &DoubleIsoMu17Mu8_Mu17leg, &b_DoubleIsoMu17Mu8_Mu17leg);
   //fChain->SetBranchAddress("DoubleIsoMu17Mu8_Mu8leg", &DoubleIsoMu17Mu8_Mu8leg, &b_DoubleIsoMu17Mu8_Mu8leg);
   //fChain->SetBranchAddress("DoubleIsoMu17Mu8dZ_Mu17leg", &DoubleIsoMu17Mu8dZ_Mu17leg, &b_DoubleIsoMu17Mu8dZ_Mu17leg);
   //fChain->SetBranchAddress("DoubleIsoMu17TkMu8_IsoMu17leg", &DoubleIsoMu17TkMu8_IsoMu17leg, &b_DoubleIsoMu17TkMu8_IsoMu17leg);
   //fChain->SetBranchAddress("DoubleIsoMu17TkMu8_IsoMu8leg", &DoubleIsoMu17TkMu8_IsoMu8leg, &b_DoubleIsoMu17TkMu8_IsoMu8leg);
   //fChain->SetBranchAddress("DoubleIsoMu17TkMu8_Mu17leg", &DoubleIsoMu17TkMu8_Mu17leg, &b_DoubleIsoMu17TkMu8_Mu17leg);
   //fChain->SetBranchAddress("DoubleIsoMu17TkMu8_TkMu8leg", &DoubleIsoMu17TkMu8_TkMu8leg, &b_DoubleIsoMu17TkMu8_TkMu8leg);
   //fChain->SetBranchAddress("DoubleIsoMu17TkMu8dZ_Mu17", &DoubleIsoMu17TkMu8dZ_Mu17, &b_DoubleIsoMu17TkMu8dZ_Mu17);
   //fChain->SetBranchAddress("DoubleMu30TkMu11", &DoubleMu30TkMu11, &b_DoubleMu30TkMu11);
   //fChain->SetBranchAddress("DoubleMu30TkMu11_Mu30leg", &DoubleMu30TkMu11_Mu30leg, &b_DoubleMu30TkMu11_Mu30leg);
   //fChain->SetBranchAddress("DoubleMu30TkMu11_TkMu11leg", &DoubleMu30TkMu11_TkMu11leg, &b_DoubleMu30TkMu11_TkMu11leg);
   fChain->SetBranchAddress("Glb", &Glb, &b_Glb);
   fChain->SetBranchAddress("GlbPT", &GlbPT, &b_GlbPT);
   //fChain->SetBranchAddress("HLT_TkMu50", &HLT_TkMu50, &b_HLT_TkMu50);
   //fChain->SetBranchAddress("HWWID", &HWWID, &b_HWWID);
   //fChain->SetBranchAddress("HighPt", &HighPt, &b_HighPt);
   //fChain->SetBranchAddress("IsoMu18", &IsoMu18, &b_IsoMu18);
   //fChain->SetBranchAddress("IsoMu20", &IsoMu20, &b_IsoMu20);
   //fChain->SetBranchAddress("IsoMu22", &IsoMu22, &b_IsoMu22);
   //fChain->SetBranchAddress("IsoMu22_eta2p1", &IsoMu22_eta2p1, &b_IsoMu22_eta2p1);
   fChain->SetBranchAddress("IsoMu24", &IsoMu24, &b_IsoMu24);
   //fChain->SetBranchAddress("IsoMu24_eta2p1", &IsoMu24_eta2p1, &b_IsoMu24_eta2p1);
   //fChain->SetBranchAddress("IsoMu27", &IsoMu27, &b_IsoMu27);
   //fChain->SetBranchAddress("IsoTkMu18", &IsoTkMu18, &b_IsoTkMu18);
   //fChain->SetBranchAddress("IsoTkMu20", &IsoTkMu20, &b_IsoTkMu20);
   //fChain->SetBranchAddress("IsoTkMu22", &IsoTkMu22, &b_IsoTkMu22);
   //fChain->SetBranchAddress("IsoTkMu22_eta2p1", &IsoTkMu22_eta2p1, &b_IsoTkMu22_eta2p1);
   fChain->SetBranchAddress("IsoTkMu24", &IsoTkMu24, &b_IsoTkMu24);
   //fChain->SetBranchAddress("IsoTkMu24_eta2p1", &IsoTkMu24_eta2p1, &b_IsoTkMu24_eta2p1);
   //fChain->SetBranchAddress("IsoTkMu27", &IsoTkMu27, &b_IsoTkMu27);
   //fChain->SetBranchAddress("L1sMu16", &L1sMu16, &b_L1sMu16);
   //fChain->SetBranchAddress("L1sMu18", &L1sMu18, &b_L1sMu18);
   //fChain->SetBranchAddress("L1sMu20", &L1sMu20, &b_L1sMu20);
   //fChain->SetBranchAddress("L1sMu25Eta2p1", &L1sMu25Eta2p1, &b_L1sMu25Eta2p1);
   //fChain->SetBranchAddress("L2fL1sDoubleMu114L1f0L2Filtered10OneMu", &L2fL1sDoubleMu114L1f0L2Filtered10OneMu, &b_L2fL1sDoubleMu114L1f0L2Filtered10OneMu);
   //fChain->SetBranchAddress("L2fL1sDoubleMu114L1f0OneMuL2Filtered10", &L2fL1sDoubleMu114L1f0OneMuL2Filtered10, &b_L2fL1sDoubleMu114L1f0OneMuL2Filtered10);
   //fChain->SetBranchAddress("L2fL1sMu18L1f0L2Filtered10Q", &L2fL1sMu18L1f0L2Filtered10Q, &b_L2fL1sMu18L1f0L2Filtered10Q);
   //fChain->SetBranchAddress("L2pfL1sDoubleMu114L1f0L2PreFiltered0", &L2pfL1sDoubleMu114L1f0L2PreFiltered0, &b_L2pfL1sDoubleMu114L1f0L2PreFiltered0);
   //fChain->SetBranchAddress("L3crIsoL1sMu16L1f0L2f10QL3f18QL3trkIsoFiltered0p09", &L3crIsoL1sMu16L1f0L2f10QL3f18QL3trkIsoFiltered0p09, &b_L3crIsoL1sMu16L1f0L2f10QL3f18QL3trkIsoFiltered0p09);
   //fChain->SetBranchAddress("L3fL1sDoubleMu114L1f0L2f10L3Filtered17", &L3fL1sDoubleMu114L1f0L2f10L3Filtered17, &b_L3fL1sDoubleMu114L1f0L2f10L3Filtered17);
   //fChain->SetBranchAddress("L3fL1sDoubleMu114L1f0L2f10OneMuL3Filtered17", &L3fL1sDoubleMu114L1f0L2f10OneMuL3Filtered17, &b_L3fL1sDoubleMu114L1f0L2f10OneMuL3Filtered17);
   //fChain->SetBranchAddress("L3fL1sMu16L1f0L2f10QL3Filtered18Q", &L3fL1sMu16L1f0L2f10QL3Filtered18Q, &b_L3fL1sMu16L1f0L2f10QL3Filtered18Q);
   //fChain->SetBranchAddress("L3fL1sMu16L1f0L2f10QL3Filtered18QL3pfecalIsoRhoFilteredEB0p11EE0p08", &L3fL1sMu16L1f0L2f10QL3Filtered18QL3pfecalIsoRhoFilteredEB0p11EE0p08, &b_L3fL1sMu16L1f0L2f10QL3Filtered18QL3pfecalIsoRhoFilteredEB0p11EE0p08);
   //fChain->SetBranchAddress("L3fL1sMu16L1f0L2f10QL3Filtered18QL3pfhcalIsoRhoFilteredHB0p21HE0p22", &L3fL1sMu16L1f0L2f10QL3Filtered18QL3pfhcalIsoRhoFilteredHB0p21HE0p22, &b_L3fL1sMu16L1f0L2f10QL3Filtered18QL3pfhcalIsoRhoFilteredHB0p21HE0p22);
   //fChain->SetBranchAddress("L3fL1sMu16f0TkFiltered18Q", &L3fL1sMu16f0TkFiltered18Q, &b_L3fL1sMu16f0TkFiltered18Q);
   //fChain->SetBranchAddress("L3pfL1sDoubleMu114L1f0L2pf0L3PreFiltered8", &L3pfL1sDoubleMu114L1f0L2pf0L3PreFiltered8, &b_L3pfL1sDoubleMu114L1f0L2pf0L3PreFiltered8);
   fChain->SetBranchAddress("Loose", &Loose, &b_Loose);
   fChain->SetBranchAddress("Medium", &Medium, &b_Medium);
   fChain->SetBranchAddress("Medium2016", &Medium2016, &b_Medium2016);
   //fChain->SetBranchAddress("Mu17", &Mu17, &b_Mu17);
   //fChain->SetBranchAddress("Mu17_IsoTrkVVL", &Mu17_IsoTrkVVL, &b_Mu17_IsoTrkVVL);
   //fChain->SetBranchAddress("Mu20", &Mu20, &b_Mu20);
   //fChain->SetBranchAddress("Mu23_TrkIsoVVL", &Mu23_TrkIsoVVL, &b_Mu23_TrkIsoVVL);
   //fChain->SetBranchAddress("Mu45_eta2p1", &Mu45_eta2p1, &b_Mu45_eta2p1);
   //fChain->SetBranchAddress("Mu50", &Mu50, &b_Mu50);
   //fChain->SetBranchAddress("MuIDForOutsideInTk", &MuIDForOutsideInTk, &b_MuIDForOutsideInTk);
   fChain->SetBranchAddress("PF", &PF, &b_PF);
   fChain->SetBranchAddress("TM", &TM, &b_TM);
   //fChain->SetBranchAddress("TMA", &TMA, &b_TMA);
   //fChain->SetBranchAddress("TMLSAT", &TMLSAT, &b_TMLSAT);
   //fChain->SetBranchAddress("TMLST", &TMLST, &b_TMLST);
   //fChain->SetBranchAddress("TMOSL", &TMOSL, &b_TMOSL);
   //fChain->SetBranchAddress("TMOST", &TMOST, &b_TMOST);
   //fChain->SetBranchAddress("TMOSTQual", &TMOSTQual, &b_TMOSTQual);
   fChain->SetBranchAddress("Tight2012", &Tight2012, &b_Tight2012);
   //fChain->SetBranchAddress("Track_HP", &Track_HP, &b_Track_HP);
   //fChain->SetBranchAddress("VBTF", &VBTF, &b_VBTF);
   //fChain->SetBranchAddress("VBTF_nL8", &VBTF_nL8, &b_VBTF_nL8);
   //fChain->SetBranchAddress("VBTF_nL9", &VBTF_nL9, &b_VBTF_nL9);
   //fChain->SetBranchAddress("hltL2fL1sMu10lqL1f0L2Filtered10", &hltL2fL1sMu10lqL1f0L2Filtered10, &b_hltL2fL1sMu10lqL1f0L2Filtered10);
   //fChain->SetBranchAddress("hltL2fL1sMu20L1f0L2Filtered10Q", &hltL2fL1sMu20L1f0L2Filtered10Q, &b_hltL2fL1sMu20L1f0L2Filtered10Q);
   //fChain->SetBranchAddress("hltL2fL1sMu22L1f0L2Filtered10Q", &hltL2fL1sMu22L1f0L2Filtered10Q, &b_hltL2fL1sMu22L1f0L2Filtered10Q);
   //fChain->SetBranchAddress("hltL3crIsoL1sMu18L1f0L2f10QL3f20QL3trkIsoFiltered0p09", &hltL3crIsoL1sMu18L1f0L2f10QL3f20QL3trkIsoFiltered0p09, &b_hltL3crIsoL1sMu18L1f0L2f10QL3f20QL3trkIsoFiltered0p09);
   //fChain->SetBranchAddress("hltL3crIsoL1sMu20L1f0L2f10QL3f22QL3trkIsoFiltered0p09", &hltL3crIsoL1sMu20L1f0L2f10QL3f22QL3trkIsoFiltered0p09, &b_hltL3crIsoL1sMu20L1f0L2f10QL3f22QL3trkIsoFiltered0p09);
   //fChain->SetBranchAddress("hltL3crIsoL1sMu22L1f0L2f10QL3f24QL3trkIsoFiltered0p09", &hltL3crIsoL1sMu22L1f0L2f10QL3f24QL3trkIsoFiltered0p09, &b_hltL3crIsoL1sMu22L1f0L2f10QL3f24QL3trkIsoFiltered0p09);
   //fChain->SetBranchAddress("hltL3fL1sMu10lqL1f0L2f10L3Filtered17", &hltL3fL1sMu10lqL1f0L2f10L3Filtered17, &b_hltL3fL1sMu10lqL1f0L2f10L3Filtered17);
   //fChain->SetBranchAddress("hltL3fL1sMu18L1f0L2f10QL3Filtered20Q", &hltL3fL1sMu18L1f0L2f10QL3Filtered20Q, &b_hltL3fL1sMu18L1f0L2f10QL3Filtered20Q);
   //fChain->SetBranchAddress("hltL3fL1sMu18L1f0L2f10QL3Filtered20QL3pfecalIsoRhoFilteredEB0p11EE0p08", &hltL3fL1sMu18L1f0L2f10QL3Filtered20QL3pfecalIsoRhoFilteredEB0p11EE0p08, &b_hltL3fL1sMu18L1f0L2f10QL3Filtered20QL3pfecalIsoRhoFilteredEB0p11EE0p08);
   //fChain->SetBranchAddress("hltL3fL1sMu18L1f0L2f10QL3Filtered20QL3pfhcalIsoRhoFilteredHB0p21HE0p22", &hltL3fL1sMu18L1f0L2f10QL3Filtered20QL3pfhcalIsoRhoFilteredHB0p21HE0p22, &b_hltL3fL1sMu18L1f0L2f10QL3Filtered20QL3pfhcalIsoRhoFilteredHB0p21HE0p22);
   //fChain->SetBranchAddress("hltL3fL1sMu18L1f0Tkf20QL3trkIsoFiltered0p09", &hltL3fL1sMu18L1f0Tkf20QL3trkIsoFiltered0p09, &b_hltL3fL1sMu18L1f0Tkf20QL3trkIsoFiltered0p09);
   //fChain->SetBranchAddress("hltL3fL1sMu18f0TkFiltered20Q", &hltL3fL1sMu18f0TkFiltered20Q, &b_hltL3fL1sMu18f0TkFiltered20Q);
   //fChain->SetBranchAddress("hltL3fL1sMu18f0TkFiltered20QL3pfecalIsoRhoFilteredEB0p11EE0p08", &hltL3fL1sMu18f0TkFiltered20QL3pfecalIsoRhoFilteredEB0p11EE0p08, &b_hltL3fL1sMu18f0TkFiltered20QL3pfecalIsoRhoFilteredEB0p11EE0p08);
   //fChain->SetBranchAddress("hltL3fL1sMu18f0TkFiltered20QL3pfhcalIsoRhoFilteredHB0p21HE0p22", &hltL3fL1sMu18f0TkFiltered20QL3pfhcalIsoRhoFilteredHB0p21HE0p22, &b_hltL3fL1sMu18f0TkFiltered20QL3pfhcalIsoRhoFilteredHB0p21HE0p22);
   //fChain->SetBranchAddress("hltL3fL1sMu1lqL1f0L2f10L3Filtered17TkIsoFiltered0p4", &hltL3fL1sMu1lqL1f0L2f10L3Filtered17TkIsoFiltered0p4, &b_hltL3fL1sMu1lqL1f0L2f10L3Filtered17TkIsoFiltered0p4);
   //fChain->SetBranchAddress("hltL3fL1sMu20L1f0L2f10QL3Filtered22Q", &hltL3fL1sMu20L1f0L2f10QL3Filtered22Q, &b_hltL3fL1sMu20L1f0L2f10QL3Filtered22Q);
   //fChain->SetBranchAddress("hltL3fL1sMu20L1f0L2f10QL3Filtered22QL3pfecalIsoRhoFilteredEB0p11EE0p08", &hltL3fL1sMu20L1f0L2f10QL3Filtered22QL3pfecalIsoRhoFilteredEB0p11EE0p08, &b_hltL3fL1sMu20L1f0L2f10QL3Filtered22QL3pfecalIsoRhoFilteredEB0p11EE0p08);
   //fChain->SetBranchAddress("hltL3fL1sMu20L1f0L2f10QL3Filtered22QL3pfhcalIsoRhoFilteredHB0p21HE0p22", &hltL3fL1sMu20L1f0L2f10QL3Filtered22QL3pfhcalIsoRhoFilteredHB0p21HE0p22, &b_hltL3fL1sMu20L1f0L2f10QL3Filtered22QL3pfhcalIsoRhoFilteredHB0p21HE0p22);
   //fChain->SetBranchAddress("hltL3fL1sMu20L1f0Tkf22QL3trkIsoFiltered0p09", &hltL3fL1sMu20L1f0Tkf22QL3trkIsoFiltered0p09, &b_hltL3fL1sMu20L1f0Tkf22QL3trkIsoFiltered0p09);
   //fChain->SetBranchAddress("hltL3fL1sMu20f0TkFiltered22Q", &hltL3fL1sMu20f0TkFiltered22Q, &b_hltL3fL1sMu20f0TkFiltered22Q);
   //fChain->SetBranchAddress("hltL3fL1sMu20f0TkFiltered22QL3pfecalIsoRhoFilteredEB0p11EE0p08", &hltL3fL1sMu20f0TkFiltered22QL3pfecalIsoRhoFilteredEB0p11EE0p08, &b_hltL3fL1sMu20f0TkFiltered22QL3pfecalIsoRhoFilteredEB0p11EE0p08);
   //fChain->SetBranchAddress("hltL3fL1sMu20f0TkFiltered22QL3pfhcalIsoRhoFilteredHB0p21HE0p22", &hltL3fL1sMu20f0TkFiltered22QL3pfhcalIsoRhoFilteredHB0p21HE0p22, &b_hltL3fL1sMu20f0TkFiltered22QL3pfhcalIsoRhoFilteredHB0p21HE0p22);
   //fChain->SetBranchAddress("hltL3fL1sMu22L1f0L2f10QL3Filtered24Q", &hltL3fL1sMu22L1f0L2f10QL3Filtered24Q, &b_hltL3fL1sMu22L1f0L2f10QL3Filtered24Q);
   //fChain->SetBranchAddress("hltL3fL1sMu22L1f0L2f10QL3Filtered24QL3pfecalIsoRhoFilteredEB0p11EE0p08", &hltL3fL1sMu22L1f0L2f10QL3Filtered24QL3pfecalIsoRhoFilteredEB0p11EE0p08, &b_hltL3fL1sMu22L1f0L2f10QL3Filtered24QL3pfecalIsoRhoFilteredEB0p11EE0p08);
   //fChain->SetBranchAddress("hltL3fL1sMu22L1f0L2f10QL3Filtered24QL3pfhcalIsoRhoFilteredHB0p21HE0p22", &hltL3fL1sMu22L1f0L2f10QL3Filtered24QL3pfhcalIsoRhoFilteredHB0p21HE0p22, &b_hltL3fL1sMu22L1f0L2f10QL3Filtered24QL3pfhcalIsoRhoFilteredHB0p21HE0p22);
   //fChain->SetBranchAddress("hltL3fL1sMu22L1f0Tkf24QL3trkIsoFiltered0p09", &hltL3fL1sMu22L1f0Tkf24QL3trkIsoFiltered0p09, &b_hltL3fL1sMu22L1f0Tkf24QL3trkIsoFiltered0p09);
   //fChain->SetBranchAddress("hltL3fL1sMu22f0TkFiltered24Q", &hltL3fL1sMu22f0TkFiltered24Q, &b_hltL3fL1sMu22f0TkFiltered24Q);
   //fChain->SetBranchAddress("hltL3fL1sMu22f0TkFiltered24QL3pfecalIsoRhoFilteredEB0p11EE0p08", &hltL3fL1sMu22f0TkFiltered24QL3pfecalIsoRhoFilteredEB0p11EE0p08, &b_hltL3fL1sMu22f0TkFiltered24QL3pfecalIsoRhoFilteredEB0p11EE0p08);
   //fChain->SetBranchAddress("hltL3fL1sMu22f0TkFiltered24QL3pfhcalIsoRhoFilteredHB0p21HE0p22", &hltL3fL1sMu22f0TkFiltered24QL3pfhcalIsoRhoFilteredHB0p21HE0p22, &b_hltL3fL1sMu22f0TkFiltered24QL3pfhcalIsoRhoFilteredHB0p21HE0p22);
   //fChain->SetBranchAddress("tkHighPt", &tkHighPt, &b_tkHighPt);
   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("lumi", &lumi, &b_lumi);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("truePU", &truePU, &b_truePU);
   fChain->SetBranchAddress("mass", &mass, &b_mass);
   if(isMC) fChain->SetBranchAddress("mcTrue", &mcTrue, &b_mcTrue);
   if(isMC) fChain->SetBranchAddress("mcMass", &mcMass, &b_mcMass);
   //fChain->SetBranchAddress("tag_IP", &tag_IP, &b_tag_IP);
   //fChain->SetBranchAddress("tag_IPError", &tag_IPError, &b_tag_IPError);
   //fChain->SetBranchAddress("tag_SIP", &tag_SIP, &b_tag_SIP);
   fChain->SetBranchAddress("tag_abseta", &tag_abseta, &b_tag_abseta);
   //fChain->SetBranchAddress("tag_caloCompatibility", &tag_caloCompatibility, &b_tag_caloCompatibility);
   fChain->SetBranchAddress("tag_charge", &tag_charge, &b_tag_charge);
   //fChain->SetBranchAddress("tag_chargedHadIso03", &tag_chargedHadIso03, &b_tag_chargedHadIso03);
   //fChain->SetBranchAddress("tag_chargedHadIso04", &tag_chargedHadIso04, &b_tag_chargedHadIso04);
   //fChain->SetBranchAddress("tag_chargedParticleIso03", &tag_chargedParticleIso03, &b_tag_chargedParticleIso03);
   //fChain->SetBranchAddress("tag_chargedParticleIso04", &tag_chargedParticleIso04, &b_tag_chargedParticleIso04);
   //fChain->SetBranchAddress("tag_chi2LocMom", &tag_chi2LocMom, &b_tag_chi2LocMom);
   //fChain->SetBranchAddress("tag_chi2LocPos", &tag_chi2LocPos, &b_tag_chi2LocPos);
   //fChain->SetBranchAddress("tag_combRelIso", &tag_combRelIso, &b_tag_combRelIso);
   //fChain->SetBranchAddress("tag_combRelIsoPF03", &tag_combRelIsoPF03, &b_tag_combRelIsoPF03);
   //fChain->SetBranchAddress("tag_combRelIsoPF03dBeta", &tag_combRelIsoPF03dBeta, &b_tag_combRelIsoPF03dBeta);
   //fChain->SetBranchAddress("tag_combRelIsoPF04", &tag_combRelIsoPF04, &b_tag_combRelIsoPF04);
   fChain->SetBranchAddress("tag_combRelIsoPF04dBeta", &tag_combRelIsoPF04dBeta, &b_tag_combRelIsoPF04dBeta);
   //fChain->SetBranchAddress("tag_dB", &tag_dB, &b_tag_dB);
   //fChain->SetBranchAddress("tag_ecalIso", &tag_ecalIso, &b_tag_ecalIso);
   //fChain->SetBranchAddress("tag_edB", &tag_edB, &b_tag_edB);
   //fChain->SetBranchAddress("tag_emEnergy", &tag_emEnergy, &b_tag_emEnergy);
   //fChain->SetBranchAddress("tag_emS9Energy", &tag_emS9Energy, &b_tag_emS9Energy);
   fChain->SetBranchAddress("tag_eta", &tag_eta, &b_tag_eta);
   //fChain->SetBranchAddress("tag_glbChi2", &tag_glbChi2, &b_tag_glbChi2);
   //fChain->SetBranchAddress("tag_glbPtError", &tag_glbPtError, &b_tag_glbPtError);
   //fChain->SetBranchAddress("tag_glbSigmaPtOverPt", &tag_glbSigmaPtOverPt, &b_tag_glbSigmaPtOverPt);
   //fChain->SetBranchAddress("tag_glbTrackProb", &tag_glbTrackProb, &b_tag_glbTrackProb);
   //fChain->SetBranchAddress("tag_glbValidMuHits", &tag_glbValidMuHits, &b_tag_glbValidMuHits);
   //fChain->SetBranchAddress("tag_hadEnergy", &tag_hadEnergy, &b_tag_hadEnergy);
   //fChain->SetBranchAddress("tag_hadS9Energy", &tag_hadS9Energy, &b_tag_hadS9Energy);
   //fChain->SetBranchAddress("tag_hcalIso", &tag_hcalIso, &b_tag_hcalIso);
   //fChain->SetBranchAddress("tag_l1dphi", &tag_l1dphi, &b_tag_l1dphi);
   //fChain->SetBranchAddress("tag_l1dr", &tag_l1dr, &b_tag_l1dr);
   //fChain->SetBranchAddress("tag_l1drByQ", &tag_l1drByQ, &b_tag_l1drByQ);
   //fChain->SetBranchAddress("tag_l1eta", &tag_l1eta, &b_tag_l1eta);
   //fChain->SetBranchAddress("tag_l1phi", &tag_l1phi, &b_tag_l1phi);
   //fChain->SetBranchAddress("tag_l1pt", &tag_l1pt, &b_tag_l1pt);
   //fChain->SetBranchAddress("tag_l1ptByQ", &tag_l1ptByQ, &b_tag_l1ptByQ);
   //fChain->SetBranchAddress("tag_l1q", &tag_l1q, &b_tag_l1q);
   //fChain->SetBranchAddress("tag_l1qByQ", &tag_l1qByQ, &b_tag_l1qByQ);
   //fChain->SetBranchAddress("tag_l2dr", &tag_l2dr, &b_tag_l2dr);
   //fChain->SetBranchAddress("tag_l2eta", &tag_l2eta, &b_tag_l2eta);
   //fChain->SetBranchAddress("tag_l2pt", &tag_l2pt, &b_tag_l2pt);
   //fChain->SetBranchAddress("tag_l3dr", &tag_l3dr, &b_tag_l3dr);
   //fChain->SetBranchAddress("tag_l3pt", &tag_l3pt, &b_tag_l3pt);
   //fChain->SetBranchAddress("tag_neutralHadIso03", &tag_neutralHadIso03, &b_tag_neutralHadIso03);
   //fChain->SetBranchAddress("tag_neutralHadIso04", &tag_neutralHadIso04, &b_tag_neutralHadIso04);
   //fChain->SetBranchAddress("tag_numberOfMatchedStations", &tag_numberOfMatchedStations, &b_tag_numberOfMatchedStations);
   //fChain->SetBranchAddress("tag_numberOfMatches", &tag_numberOfMatches, &b_tag_numberOfMatches);
   fChain->SetBranchAddress("tag_p", &tag_p, &b_tag_p);
   fChain->SetBranchAddress("tag_phi", &tag_phi, &b_tag_phi);
   //fChain->SetBranchAddress("tag_photonIso03", &tag_photonIso03, &b_tag_photonIso03);
   //fChain->SetBranchAddress("tag_photonIso04", &tag_photonIso04, &b_tag_photonIso04);
   fChain->SetBranchAddress("tag_pt", &tag_pt, &b_tag_pt);
   //fChain->SetBranchAddress("tag_puIso03", &tag_puIso03, &b_tag_puIso03);
   //fChain->SetBranchAddress("tag_puIso04", &tag_puIso04, &b_tag_puIso04);
   //fChain->SetBranchAddress("tag_relEcalIso", &tag_relEcalIso, &b_tag_relEcalIso);
   //fChain->SetBranchAddress("tag_relHcalIso", &tag_relHcalIso, &b_tag_relHcalIso);
   //fChain->SetBranchAddress("tag_relTkIso", &tag_relTkIso, &b_tag_relTkIso);
   //fChain->SetBranchAddress("tag_segmentCompatibility", &tag_segmentCompatibility, &b_tag_segmentCompatibility);
   //fChain->SetBranchAddress("tag_staQoverP", &tag_staQoverP, &b_tag_staQoverP);
   //fChain->SetBranchAddress("tag_staQoverPerror", &tag_staQoverPerror, &b_tag_staQoverPerror);
   //fChain->SetBranchAddress("tag_staValidStations", &tag_staValidStations, &b_tag_staValidStations);
   //fChain->SetBranchAddress("tag_tkChi2", &tag_tkChi2, &b_tag_tkChi2);
   //fChain->SetBranchAddress("tag_tkExpHitIn", &tag_tkExpHitIn, &b_tag_tkExpHitIn);
   //fChain->SetBranchAddress("tag_tkExpHitOut", &tag_tkExpHitOut, &b_tag_tkExpHitOut);
   //fChain->SetBranchAddress("tag_tkHitFract", &tag_tkHitFract, &b_tag_tkHitFract);
   //fChain->SetBranchAddress("tag_tkIso", &tag_tkIso, &b_tag_tkIso);
   //fChain->SetBranchAddress("tag_tkKink", &tag_tkKink, &b_tag_tkKink);
   //fChain->SetBranchAddress("tag_tkPixelLay", &tag_tkPixelLay, &b_tag_tkPixelLay);
   //fChain->SetBranchAddress("tag_tkPtError", &tag_tkPtError, &b_tag_tkPtError);
   //fChain->SetBranchAddress("tag_tkSigmaPtOverPt", &tag_tkSigmaPtOverPt, &b_tag_tkSigmaPtOverPt);
   //fChain->SetBranchAddress("tag_tkTrackerLay", &tag_tkTrackerLay, &b_tag_tkTrackerLay);
   //fChain->SetBranchAddress("tag_tkValidHits", &tag_tkValidHits, &b_tag_tkValidHits);
   //fChain->SetBranchAddress("tag_tkValidPixelHits", &tag_tkValidPixelHits, &b_tag_tkValidPixelHits);
   //fChain->SetBranchAddress("tag_dxyBS", &tag_dxyBS, &b_tag_dxyBS);
   //fChain->SetBranchAddress("tag_dxyPVdzmin", &tag_dxyPVdzmin, &b_tag_dxyPVdzmin);
   //fChain->SetBranchAddress("tag_dzPV", &tag_dzPV, &b_tag_dzPV);
   //fChain->SetBranchAddress("tag_fixedGridRhoFastjetAll", &tag_fixedGridRhoFastjetAll, &b_tag_fixedGridRhoFastjetAll);
   //fChain->SetBranchAddress("tag_fixedGridRhoFastjetAllCalo", &tag_fixedGridRhoFastjetAllCalo, &b_tag_fixedGridRhoFastjetAllCalo);
   //fChain->SetBranchAddress("tag_fixedGridRhoFastjetCentralCalo", &tag_fixedGridRhoFastjetCentralCalo, &b_tag_fixedGridRhoFastjetCentralCalo);
   //fChain->SetBranchAddress("tag_fixedGridRhoFastjetCentralChargedPileUp", &tag_fixedGridRhoFastjetCentralChargedPileUp, &b_tag_fixedGridRhoFastjetCentralChargedPileUp);
   //fChain->SetBranchAddress("tag_fixedGridRhoFastjetCentralNeutral", &tag_fixedGridRhoFastjetCentralNeutral, &b_tag_fixedGridRhoFastjetCentralNeutral);
   //fChain->SetBranchAddress("tag_isoTrk03Abs", &tag_isoTrk03Abs, &b_tag_isoTrk03Abs);
   //fChain->SetBranchAddress("tag_isoTrk03Rel", &tag_isoTrk03Rel, &b_tag_isoTrk03Rel);
   //fChain->SetBranchAddress("tag_met", &tag_met, &b_tag_met);
   //fChain->SetBranchAddress("tag_mt", &tag_mt, &b_tag_mt);
   //fChain->SetBranchAddress("tag_nSplitTk", &tag_nSplitTk, &b_tag_nSplitTk);
   fChain->SetBranchAddress("tag_nVertices", &tag_nVertices, &b_tag_nVertices);
   //fChain->SetBranchAddress("tag_DiMuonGlb17Glb8RelTrkIsoFiltered0p4", &tag_DiMuonGlb17Glb8RelTrkIsoFiltered0p4, &b_tag_DiMuonGlb17Glb8RelTrkIsoFiltered0p4);
   //fChain->SetBranchAddress("tag_DoubleIsoMu17Mu8_IsoMu17leg", &tag_DoubleIsoMu17Mu8_IsoMu17leg, &b_tag_DoubleIsoMu17Mu8_IsoMu17leg);
   //fChain->SetBranchAddress("tag_DoubleIsoMu17Mu8_IsoMu8leg", &tag_DoubleIsoMu17Mu8_IsoMu8leg, &b_tag_DoubleIsoMu17Mu8_IsoMu8leg);
   //fChain->SetBranchAddress("tag_DoubleIsoMu17Mu8_Mu17leg", &tag_DoubleIsoMu17Mu8_Mu17leg, &b_tag_DoubleIsoMu17Mu8_Mu17leg);
   //fChain->SetBranchAddress("tag_DoubleIsoMu17Mu8_Mu8leg", &tag_DoubleIsoMu17Mu8_Mu8leg, &b_tag_DoubleIsoMu17Mu8_Mu8leg);
   //fChain->SetBranchAddress("tag_DoubleIsoMu17Mu8dZ_Mu17leg", &tag_DoubleIsoMu17Mu8dZ_Mu17leg, &b_tag_DoubleIsoMu17Mu8dZ_Mu17leg);
   //fChain->SetBranchAddress("tag_DoubleIsoMu17TkMu8_IsoMu17leg", &tag_DoubleIsoMu17TkMu8_IsoMu17leg, &b_tag_DoubleIsoMu17TkMu8_IsoMu17leg);
   //fChain->SetBranchAddress("tag_DoubleIsoMu17TkMu8_IsoMu8leg", &tag_DoubleIsoMu17TkMu8_IsoMu8leg, &b_tag_DoubleIsoMu17TkMu8_IsoMu8leg);
   //fChain->SetBranchAddress("tag_DoubleIsoMu17TkMu8_Mu17leg", &tag_DoubleIsoMu17TkMu8_Mu17leg, &b_tag_DoubleIsoMu17TkMu8_Mu17leg);
   //fChain->SetBranchAddress("tag_DoubleIsoMu17TkMu8_TkMu8leg", &tag_DoubleIsoMu17TkMu8_TkMu8leg, &b_tag_DoubleIsoMu17TkMu8_TkMu8leg);
   //fChain->SetBranchAddress("tag_DoubleIsoMu17TkMu8dZ_Mu17", &tag_DoubleIsoMu17TkMu8dZ_Mu17, &b_tag_DoubleIsoMu17TkMu8dZ_Mu17);
   //fChain->SetBranchAddress("tag_DoubleMu30TkMu11", &tag_DoubleMu30TkMu11, &b_tag_DoubleMu30TkMu11);
   //fChain->SetBranchAddress("tag_DoubleMu30TkMu11_Mu30leg", &tag_DoubleMu30TkMu11_Mu30leg, &b_tag_DoubleMu30TkMu11_Mu30leg);
   //fChain->SetBranchAddress("tag_DoubleMu30TkMu11_TkMu11leg", &tag_DoubleMu30TkMu11_TkMu11leg, &b_tag_DoubleMu30TkMu11_TkMu11leg);
   //fChain->SetBranchAddress("tag_HLT_TkMu50", &tag_HLT_TkMu50, &b_tag_HLT_TkMu50);
   //fChain->SetBranchAddress("tag_IsoMu18", &tag_IsoMu18, &b_tag_IsoMu18);
   //fChain->SetBranchAddress("tag_IsoMu20", &tag_IsoMu20, &b_tag_IsoMu20);
   //fChain->SetBranchAddress("tag_IsoMu22", &tag_IsoMu22, &b_tag_IsoMu22);
   //fChain->SetBranchAddress("tag_IsoMu22_eta2p1", &tag_IsoMu22_eta2p1, &b_tag_IsoMu22_eta2p1);
   fChain->SetBranchAddress("tag_IsoMu24", &tag_IsoMu24, &b_tag_IsoMu24);
   //fChain->SetBranchAddress("tag_IsoMu24_eta2p1", &tag_IsoMu24_eta2p1, &b_tag_IsoMu24_eta2p1);
   //fChain->SetBranchAddress("tag_IsoMu27", &tag_IsoMu27, &b_tag_IsoMu27);
   //fChain->SetBranchAddress("tag_IsoTkMu18", &tag_IsoTkMu18, &b_tag_IsoTkMu18);
   //fChain->SetBranchAddress("tag_IsoTkMu20", &tag_IsoTkMu20, &b_tag_IsoTkMu20);
   //fChain->SetBranchAddress("tag_IsoTkMu22", &tag_IsoTkMu22, &b_tag_IsoTkMu22);
   //fChain->SetBranchAddress("tag_IsoTkMu22_eta2p1", &tag_IsoTkMu22_eta2p1, &b_tag_IsoTkMu22_eta2p1);
   fChain->SetBranchAddress("tag_IsoTkMu24", &tag_IsoTkMu24, &b_tag_IsoTkMu24);
   //fChain->SetBranchAddress("tag_IsoTkMu24_eta2p1", &tag_IsoTkMu24_eta2p1, &b_tag_IsoTkMu24_eta2p1);
   //fChain->SetBranchAddress("tag_IsoTkMu27", &tag_IsoTkMu27, &b_tag_IsoTkMu27);
   //fChain->SetBranchAddress("tag_L1sMu16", &tag_L1sMu16, &b_tag_L1sMu16);
   //fChain->SetBranchAddress("tag_L1sMu18", &tag_L1sMu18, &b_tag_L1sMu18);
   //fChain->SetBranchAddress("tag_L1sMu20", &tag_L1sMu20, &b_tag_L1sMu20);
   //fChain->SetBranchAddress("tag_L1sMu25Eta2p1", &tag_L1sMu25Eta2p1, &b_tag_L1sMu25Eta2p1);
   //fChain->SetBranchAddress("tag_L2fL1sDoubleMu114L1f0L2Filtered10OneMu", &tag_L2fL1sDoubleMu114L1f0L2Filtered10OneMu, &b_tag_L2fL1sDoubleMu114L1f0L2Filtered10OneMu);
   //fChain->SetBranchAddress("tag_L2fL1sDoubleMu114L1f0OneMuL2Filtered10", &tag_L2fL1sDoubleMu114L1f0OneMuL2Filtered10, &b_tag_L2fL1sDoubleMu114L1f0OneMuL2Filtered10);
   //fChain->SetBranchAddress("tag_L2fL1sMu18L1f0L2Filtered10Q", &tag_L2fL1sMu18L1f0L2Filtered10Q, &b_tag_L2fL1sMu18L1f0L2Filtered10Q);
   //fChain->SetBranchAddress("tag_L2pfL1sDoubleMu114L1f0L2PreFiltered0", &tag_L2pfL1sDoubleMu114L1f0L2PreFiltered0, &b_tag_L2pfL1sDoubleMu114L1f0L2PreFiltered0);
   //fChain->SetBranchAddress("tag_L3crIsoL1sMu16L1f0L2f10QL3f18QL3trkIsoFiltered0p09", &tag_L3crIsoL1sMu16L1f0L2f10QL3f18QL3trkIsoFiltered0p09, &b_tag_L3crIsoL1sMu16L1f0L2f10QL3f18QL3trkIsoFiltered0p09);
   //fChain->SetBranchAddress("tag_L3fL1sDoubleMu114L1f0L2f10L3Filtered17", &tag_L3fL1sDoubleMu114L1f0L2f10L3Filtered17, &b_tag_L3fL1sDoubleMu114L1f0L2f10L3Filtered17);
   //fChain->SetBranchAddress("tag_L3fL1sDoubleMu114L1f0L2f10OneMuL3Filtered17", &tag_L3fL1sDoubleMu114L1f0L2f10OneMuL3Filtered17, &b_tag_L3fL1sDoubleMu114L1f0L2f10OneMuL3Filtered17);
   //fChain->SetBranchAddress("tag_L3fL1sMu16L1f0L2f10QL3Filtered18Q", &tag_L3fL1sMu16L1f0L2f10QL3Filtered18Q, &b_tag_L3fL1sMu16L1f0L2f10QL3Filtered18Q);
   //fChain->SetBranchAddress("tag_L3fL1sMu16L1f0L2f10QL3Filtered18QL3pfecalIsoRhoFilteredEB0p11EE0p08", &tag_L3fL1sMu16L1f0L2f10QL3Filtered18QL3pfecalIsoRhoFilteredEB0p11EE0p08, &b_tag_L3fL1sMu16L1f0L2f10QL3Filtered18QL3pfecalIsoRhoFilteredEB0p11EE0p08);
   //fChain->SetBranchAddress("tag_L3fL1sMu16L1f0L2f10QL3Filtered18QL3pfhcalIsoRhoFilteredHB0p21HE0p22", &tag_L3fL1sMu16L1f0L2f10QL3Filtered18QL3pfhcalIsoRhoFilteredHB0p21HE0p22, &b_tag_L3fL1sMu16L1f0L2f10QL3Filtered18QL3pfhcalIsoRhoFilteredHB0p21HE0p22);
   //fChain->SetBranchAddress("tag_L3fL1sMu16f0TkFiltered18Q", &tag_L3fL1sMu16f0TkFiltered18Q, &b_tag_L3fL1sMu16f0TkFiltered18Q);
   //fChain->SetBranchAddress("tag_L3pfL1sDoubleMu114L1f0L2pf0L3PreFiltered8", &tag_L3pfL1sDoubleMu114L1f0L2pf0L3PreFiltered8, &b_tag_L3pfL1sDoubleMu114L1f0L2pf0L3PreFiltered8);
   //fChain->SetBranchAddress("tag_Mu17", &tag_Mu17, &b_tag_Mu17);
   //fChain->SetBranchAddress("tag_Mu17_IsoTrkVVL", &tag_Mu17_IsoTrkVVL, &b_tag_Mu17_IsoTrkVVL);
   //fChain->SetBranchAddress("tag_Mu20", &tag_Mu20, &b_tag_Mu20);
   //fChain->SetBranchAddress("tag_Mu23_TrkIsoVVL", &tag_Mu23_TrkIsoVVL, &b_tag_Mu23_TrkIsoVVL);
   //fChain->SetBranchAddress("tag_Mu45_eta2p1", &tag_Mu45_eta2p1, &b_tag_Mu45_eta2p1);
   //fChain->SetBranchAddress("tag_Mu50", &tag_Mu50, &b_tag_Mu50);
   //fChain->SetBranchAddress("tag_hltL2fL1sMu10lqL1f0L2Filtered10", &tag_hltL2fL1sMu10lqL1f0L2Filtered10, &b_tag_hltL2fL1sMu10lqL1f0L2Filtered10);
   //fChain->SetBranchAddress("tag_hltL2fL1sMu20L1f0L2Filtered10Q", &tag_hltL2fL1sMu20L1f0L2Filtered10Q, &b_tag_hltL2fL1sMu20L1f0L2Filtered10Q);
   //fChain->SetBranchAddress("tag_hltL2fL1sMu22L1f0L2Filtered10Q", &tag_hltL2fL1sMu22L1f0L2Filtered10Q, &b_tag_hltL2fL1sMu22L1f0L2Filtered10Q);
   //fChain->SetBranchAddress("tag_hltL3crIsoL1sMu18L1f0L2f10QL3f20QL3trkIsoFiltered0p09", &tag_hltL3crIsoL1sMu18L1f0L2f10QL3f20QL3trkIsoFiltered0p09, &b_tag_hltL3crIsoL1sMu18L1f0L2f10QL3f20QL3trkIsoFiltered0p09);
   //fChain->SetBranchAddress("tag_hltL3crIsoL1sMu20L1f0L2f10QL3f22QL3trkIsoFiltered0p09", &tag_hltL3crIsoL1sMu20L1f0L2f10QL3f22QL3trkIsoFiltered0p09, &b_tag_hltL3crIsoL1sMu20L1f0L2f10QL3f22QL3trkIsoFiltered0p09);
   //fChain->SetBranchAddress("tag_hltL3crIsoL1sMu22L1f0L2f10QL3f24QL3trkIsoFiltered0p09", &tag_hltL3crIsoL1sMu22L1f0L2f10QL3f24QL3trkIsoFiltered0p09, &b_tag_hltL3crIsoL1sMu22L1f0L2f10QL3f24QL3trkIsoFiltered0p09);
   //fChain->SetBranchAddress("tag_hltL3fL1sMu10lqL1f0L2f10L3Filtered17", &tag_hltL3fL1sMu10lqL1f0L2f10L3Filtered17, &b_tag_hltL3fL1sMu10lqL1f0L2f10L3Filtered17);
   //fChain->SetBranchAddress("tag_hltL3fL1sMu18L1f0L2f10QL3Filtered20Q", &tag_hltL3fL1sMu18L1f0L2f10QL3Filtered20Q, &b_tag_hltL3fL1sMu18L1f0L2f10QL3Filtered20Q);
   //fChain->SetBranchAddress("tag_hltL3fL1sMu18L1f0L2f10QL3Filtered20QL3pfecalIsoRhoFilteredEB0p11EE0p08", &tag_hltL3fL1sMu18L1f0L2f10QL3Filtered20QL3pfecalIsoRhoFilteredEB0p11EE0p08, &b_tag_hltL3fL1sMu18L1f0L2f10QL3Filtered20QL3pfecalIsoRhoFilteredEB0p11EE0p08);
   //fChain->SetBranchAddress("tag_hltL3fL1sMu18L1f0L2f10QL3Filtered20QL3pfhcalIsoRhoFilteredHB0p21HE0p22", &tag_hltL3fL1sMu18L1f0L2f10QL3Filtered20QL3pfhcalIsoRhoFilteredHB0p21HE0p22, &b_tag_hltL3fL1sMu18L1f0L2f10QL3Filtered20QL3pfhcalIsoRhoFilteredHB0p21HE0p22);
   //fChain->SetBranchAddress("tag_hltL3fL1sMu18L1f0Tkf20QL3trkIsoFiltered0p09", &tag_hltL3fL1sMu18L1f0Tkf20QL3trkIsoFiltered0p09, &b_tag_hltL3fL1sMu18L1f0Tkf20QL3trkIsoFiltered0p09);
   //fChain->SetBranchAddress("tag_hltL3fL1sMu18f0TkFiltered20Q", &tag_hltL3fL1sMu18f0TkFiltered20Q, &b_tag_hltL3fL1sMu18f0TkFiltered20Q);
   //fChain->SetBranchAddress("tag_hltL3fL1sMu18f0TkFiltered20QL3pfecalIsoRhoFilteredEB0p11EE0p08", &tag_hltL3fL1sMu18f0TkFiltered20QL3pfecalIsoRhoFilteredEB0p11EE0p08, &b_tag_hltL3fL1sMu18f0TkFiltered20QL3pfecalIsoRhoFilteredEB0p11EE0p08);
   //fChain->SetBranchAddress("tag_hltL3fL1sMu18f0TkFiltered20QL3pfhcalIsoRhoFilteredHB0p21HE0p22", &tag_hltL3fL1sMu18f0TkFiltered20QL3pfhcalIsoRhoFilteredHB0p21HE0p22, &b_tag_hltL3fL1sMu18f0TkFiltered20QL3pfhcalIsoRhoFilteredHB0p21HE0p22);
   //fChain->SetBranchAddress("tag_hltL3fL1sMu1lqL1f0L2f10L3Filtered17TkIsoFiltered0p4", &tag_hltL3fL1sMu1lqL1f0L2f10L3Filtered17TkIsoFiltered0p4, &b_tag_hltL3fL1sMu1lqL1f0L2f10L3Filtered17TkIsoFiltered0p4);
   //fChain->SetBranchAddress("tag_hltL3fL1sMu20L1f0L2f10QL3Filtered22Q", &tag_hltL3fL1sMu20L1f0L2f10QL3Filtered22Q, &b_tag_hltL3fL1sMu20L1f0L2f10QL3Filtered22Q);
   //fChain->SetBranchAddress("tag_hltL3fL1sMu20L1f0L2f10QL3Filtered22QL3pfecalIsoRhoFilteredEB0p11EE0p08", &tag_hltL3fL1sMu20L1f0L2f10QL3Filtered22QL3pfecalIsoRhoFilteredEB0p11EE0p08, &b_tag_hltL3fL1sMu20L1f0L2f10QL3Filtered22QL3pfecalIsoRhoFilteredEB0p11EE0p08);
   //fChain->SetBranchAddress("tag_hltL3fL1sMu20L1f0L2f10QL3Filtered22QL3pfhcalIsoRhoFilteredHB0p21HE0p22", &tag_hltL3fL1sMu20L1f0L2f10QL3Filtered22QL3pfhcalIsoRhoFilteredHB0p21HE0p22, &b_tag_hltL3fL1sMu20L1f0L2f10QL3Filtered22QL3pfhcalIsoRhoFilteredHB0p21HE0p22);
   //fChain->SetBranchAddress("tag_hltL3fL1sMu20L1f0Tkf22QL3trkIsoFiltered0p09", &tag_hltL3fL1sMu20L1f0Tkf22QL3trkIsoFiltered0p09, &b_tag_hltL3fL1sMu20L1f0Tkf22QL3trkIsoFiltered0p09);
   //fChain->SetBranchAddress("tag_hltL3fL1sMu20f0TkFiltered22Q", &tag_hltL3fL1sMu20f0TkFiltered22Q, &b_tag_hltL3fL1sMu20f0TkFiltered22Q);
   //fChain->SetBranchAddress("tag_hltL3fL1sMu20f0TkFiltered22QL3pfecalIsoRhoFilteredEB0p11EE0p08", &tag_hltL3fL1sMu20f0TkFiltered22QL3pfecalIsoRhoFilteredEB0p11EE0p08, &b_tag_hltL3fL1sMu20f0TkFiltered22QL3pfecalIsoRhoFilteredEB0p11EE0p08);
   //fChain->SetBranchAddress("tag_hltL3fL1sMu20f0TkFiltered22QL3pfhcalIsoRhoFilteredHB0p21HE0p22", &tag_hltL3fL1sMu20f0TkFiltered22QL3pfhcalIsoRhoFilteredHB0p21HE0p22, &b_tag_hltL3fL1sMu20f0TkFiltered22QL3pfhcalIsoRhoFilteredHB0p21HE0p22);
   //fChain->SetBranchAddress("tag_hltL3fL1sMu22L1f0L2f10QL3Filtered24Q", &tag_hltL3fL1sMu22L1f0L2f10QL3Filtered24Q, &b_tag_hltL3fL1sMu22L1f0L2f10QL3Filtered24Q);
   //fChain->SetBranchAddress("tag_hltL3fL1sMu22L1f0L2f10QL3Filtered24QL3pfecalIsoRhoFilteredEB0p11EE0p08", &tag_hltL3fL1sMu22L1f0L2f10QL3Filtered24QL3pfecalIsoRhoFilteredEB0p11EE0p08, &b_tag_hltL3fL1sMu22L1f0L2f10QL3Filtered24QL3pfecalIsoRhoFilteredEB0p11EE0p08);
   //fChain->SetBranchAddress("tag_hltL3fL1sMu22L1f0L2f10QL3Filtered24QL3pfhcalIsoRhoFilteredHB0p21HE0p22", &tag_hltL3fL1sMu22L1f0L2f10QL3Filtered24QL3pfhcalIsoRhoFilteredHB0p21HE0p22, &b_tag_hltL3fL1sMu22L1f0L2f10QL3Filtered24QL3pfhcalIsoRhoFilteredHB0p21HE0p22);
   //fChain->SetBranchAddress("tag_hltL3fL1sMu22L1f0Tkf24QL3trkIsoFiltered0p09", &tag_hltL3fL1sMu22L1f0Tkf24QL3trkIsoFiltered0p09, &b_tag_hltL3fL1sMu22L1f0Tkf24QL3trkIsoFiltered0p09);
   //fChain->SetBranchAddress("tag_hltL3fL1sMu22f0TkFiltered24Q", &tag_hltL3fL1sMu22f0TkFiltered24Q, &b_tag_hltL3fL1sMu22f0TkFiltered24Q);
   //fChain->SetBranchAddress("tag_hltL3fL1sMu22f0TkFiltered24QL3pfecalIsoRhoFilteredEB0p11EE0p08", &tag_hltL3fL1sMu22f0TkFiltered24QL3pfecalIsoRhoFilteredEB0p11EE0p08, &b_tag_hltL3fL1sMu22f0TkFiltered24QL3pfecalIsoRhoFilteredEB0p11EE0p08);
   //fChain->SetBranchAddress("tag_hltL3fL1sMu22f0TkFiltered24QL3pfhcalIsoRhoFilteredHB0p21HE0p22", &tag_hltL3fL1sMu22f0TkFiltered24QL3pfhcalIsoRhoFilteredHB0p21HE0p22, &b_tag_hltL3fL1sMu22f0TkFiltered24QL3pfhcalIsoRhoFilteredHB0p21HE0p22);
   //fChain->SetBranchAddress("pair_deltaR", &pair_deltaR, &b_pair_deltaR);
   //fChain->SetBranchAddress("pair_dz", &pair_dz, &b_pair_dz);
   fChain->SetBranchAddress("pair_pt", &pair_pt, &b_pair_pt);
   fChain->SetBranchAddress("pair_rapidity", &pair_rapidity, &b_pair_rapidity);
   if(isMC) fChain->SetBranchAddress("pair_actualPileUp", &pair_actualPileUp, &b_pair_actualPileUp);
   if(isMC) fChain->SetBranchAddress("pair_genWeight", &pair_genWeight, &b_pair_genWeight);
   //fChain->SetBranchAddress("pair_nJets30", &pair_nJets30, &b_pair_nJets30);
   //fChain->SetBranchAddress("pair_newTuneP_mass", &pair_newTuneP_mass, &b_pair_newTuneP_mass);
   //fChain->SetBranchAddress("pair_newTuneP_probe_pt", &pair_newTuneP_probe_pt, &b_pair_newTuneP_probe_pt);
   //fChain->SetBranchAddress("pair_newTuneP_probe_sigmaPtOverPt", &pair_newTuneP_probe_sigmaPtOverPt, &b_pair_newTuneP_probe_sigmaPtOverPt);
   //fChain->SetBranchAddress("pair_newTuneP_probe_trackType", &pair_newTuneP_probe_trackType, &b_pair_newTuneP_probe_trackType);
   fChain->SetBranchAddress("pair_probeMultiplicity", &pair_probeMultiplicity, &b_pair_probeMultiplicity);
   //fChain->SetBranchAddress("pair_probeMultiplicity_Pt10_M60140", &pair_probeMultiplicity_Pt10_M60140, &b_pair_probeMultiplicity_Pt10_M60140);
   //fChain->SetBranchAddress("pair_probeMultiplicity_TMGM", &pair_probeMultiplicity_TMGM, &b_pair_probeMultiplicity_TMGM);
   if(isMC) fChain->SetBranchAddress("pair_truePileUp", &pair_truePileUp, &b_pair_truePileUp);
   //fChain->SetBranchAddress("pair_BestZ", &pair_BestZ, &b_pair_BestZ);
   Notify();
}

Bool_t MuonPogFitterTree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void MuonPogFitterTree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t MuonPogFitterTree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef MuonPogFitterTree_cxx
