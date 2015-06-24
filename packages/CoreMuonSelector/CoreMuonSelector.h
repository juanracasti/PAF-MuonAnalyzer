///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////                                                                                             /////////////
/////////////                                     CORE MUON SELECTOR                                      /////////////
/////////////                                                                                             /////////////
/////////////                                  Juan R. Castiñeiras (IFCA)                                 /////////////
/////////////                                          Jun 2016                                           /////////////
/////////////                                                                                             /////////////
/////////////                              -> Adjust to a 120 width window <-                             /////////////
/////////////                                                                                             /////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include "PAF/computing/PAFChainItemSelector.h"

#include <TH1F.h>
#include <TMatrix.h>
#include <TH2F.h>
#include <TLorentzVector.h>
#include "Riostream.h"  


class CoreMuonSelector: public PAFChainItemSelector{
  
 public:
  virtual ~CoreMuonSelector() {}
  
  virtual void                Initialise();
  virtual void                InsideLoop();
  virtual void                Summary();
  
 protected:
  // See description in implementation
  void                        CheckMuons();
  void                        SetGenInfo();
  void                        GetMatching();
  bool                        passMediumID(int);
  bool                        passISO(int, string, float);
  float                       getISO(int, string);
  void                        SetEventFlags();
  void                        Counting();
  void                        doEffsRECO(int, int);
  void                        doEffsRECODilep();
  void                        ISORocCurve();


  // My Declarations:
  // Define data members

  // VARIABLES FOR EACH EVENT (to be initialized after every event)
 
  // GEN Info
  std::vector<TLorentzVector> G_GEN_PromptMuon_4vec;  //Lorentz vector for all GEN prompt muons, or coming from 
                                                      //a prompt tau decay, ordered by Pt
  std::vector<TLorentzVector> G_GEN_Muon_4vec;        //Lorentz vector for muons that are not included above
  bool                        G_GEN_isMuMu;      //There are 2 GEN prompt muons that come directly from a boson  
  bool                        G_GEN_isMuTau;     // "                        "  but the 2nd (in Pt) comes from a tau
  bool                        G_GEN_isTauMu;     // "                        "  but the 1st comes from a tau
  bool                        G_GEN_isTauTau;    // "                        "  and both come from a tau
  bool                        G_GEN_Pass;        // 1 of the 4 above is true
  bool                        G_GEN_isNonPrompt; //There is at least a GEN non prompt muon
  
  // RECO muons
  std::vector<TLorentzVector> G_Muon_4vec;            //Lorentz vector for all RECO muons ordered by Pt

  std::vector<bool>           G_MuonID_Tight;         //Which RECO muons pass Tight ID 
  std::vector<bool>           G_MuonID_Medium;        // ""          ""  pass Medium ID
  std::vector<bool>           G_MuonID_HWW;           // ""          ""  pass HWW ID
  std::vector<bool>           G_MuonID_Tight_GoT;     // ""          ""  pass Tight ID with GLB or TRK arbitrated 
  std::vector<bool>           G_MuonID_IPs_HWW;       // ""          ""  pass IP cuts (dxy and dz) from HWW ID
  std::vector<bool>           G_MuonID_GLBorTRKArb;   // ""          ""  are GLB or TRK Arbitrated
  std::vector<bool>           G_MuonID_Fiducial;      // ""          ""  have:
                                                      //   * For the first muon,  Pt > 20 & |Eta| < 2.4 
                                                      //   * For the other muons, Pt > 10 & |Eta| < 2.4    

  std::vector<bool>           G_MuonISO03;            //Which RECO muons pass PF Rel. ISO, dR=0.3 (w.p. = 0.12)
  std::vector<bool>           G_MuonISO03_dBeta;      // ""                            "", dR=0.3 and dBeta corr.
  std::vector<bool>           G_MuonISO03_PFWeighted; // ""                            "", dR=0.3 and PFweighted corr.
  std::vector<bool>           G_MuonISO03_PUPPI;      // ""                            "", dR=0.3 and PUPPI corr.
  std::vector<bool>           G_MuonISO04;            // ""                            "", dR=0.4
  std::vector<bool>           G_MuonISO04_dBeta;      // ""                            "", dR=0.4 and dBeta corr.
  std::vector<bool>           G_MuonISO04_PFWeighted; // ""                            "", dR=0.4 and PFweighted corr.
  std::vector<bool>           G_MuonISO04_PUPPI;      // ""                            "", dR=0.4 and PUPPI corr.

  std::vector<bool>           G_Muon_ChCompatible;    //Which RECO muons have opposite charge than the first one
                                                      //(true for the first one too)
  std::vector<int>            G_Muon_Matching;        //To which GEN prompt muon the RECO muons are matched
                                                      //  * 1: matched to the 1st GEN prompt muon
                                                      //  * 2: matched to the 2nd GEN prompt muon 
                                                      //  * 0: not matched to any GEN prompt muon
  bool                        G_PassMatching;         //The 1st and the 2nd RECO muons are matched to GEN muons
                                                      //This provides the maximum ID and ISO efficiency

  // Sizes
  UInt_t                      G_RecoMuSize;      //How many RECO muons there are
  UInt_t                      G_NPV;             //How many Primary Vtx. there are

  // Event Flags
  bool                        EvtFlag_Fiducial; // 1st and 2nd muon pass Fiducial Selection
  bool                        EvtFlag_Gen;      // Same as G_GEN_Pass
  bool                        EvtFlag_Matching; // Same as G_PassMatching


  // VARIABLES FOR ALL EVENTS (to be initialized only once)

  // Counting events

  UInt_t                      GCount_AllEvents;
  UInt_t                      GCount_GenEvents;
  UInt_t                      GCount_Fiducial_AtLeast2;
  UInt_t                      GCount_Fiducial_2;
  UInt_t                      GCount_Fiducial_1st2nd;
  UInt_t                      GCount_Match_1st2nd;
  UInt_t                      GCount_NoMatch_1st2nd;
  UInt_t                      GCount_MatchTight_1st2nd;
  UInt_t                      GCount_MatchTightIso_1st2nd;
  UInt_t                      GCount_MatchTightIso_Only1st;
  UInt_t                      GCount_MatchTightIso_Only2nd;
  UInt_t                      GCount_MatchTightIso_None;
  UInt_t                      GCount_MatchTight_Only1st;
  UInt_t                      GCount_MatchTight_Only2nd;
  UInt_t                      GCount_MatchTight_None;
  UInt_t                      GCount_Tight_1st2nd;
  UInt_t                      GCount_TightIso_1st2nd;
  UInt_t                      GCount_TightIso_Only1st;
  UInt_t                      GCount_TightIso_Only2nd;
  UInt_t                      GCount_TightIso_None;
  UInt_t                      GCount_Tight_Only1st;
  UInt_t                      GCount_Tight_Only2nd;
  UInt_t                      GCount_Tight_None;
  UInt_t                      GCount_Fiducial_1st3rd;
  UInt_t                      GCount_Match_1st3rd;
  UInt_t                      GCount_NoMatch_1st3rd;
  UInt_t                      GCount_MatchTight_1st3rd;
  UInt_t                      GCount_MatchTightIso_1st3rd;
  UInt_t                      GCount_MatchTightIso_Only3rd;
  UInt_t                      GCount_MatchTight_Only3rd;
  UInt_t                      GCount_Tight_1st3rd;
  UInt_t                      GCount_TightIso_1st3rd;
  UInt_t                      GCount_TightIso_Only3rd;
  UInt_t                      GCount_Tight_Only3rd;
  UInt_t                      GCount_Fiducial_1stOther;
  UInt_t                      GCount_Fiducial_MoreThan2;
  UInt_t                      GCount_Match_MoreThan2_OK;
  UInt_t                      GCount_Match_MoreThan2;
  UInt_t                      GCount_NoMatch_MoreThan2;
  UInt_t                      GCount_Fiducial_3;
  UInt_t                      GCount_Fiducial_MoreThan3;
  UInt_t                      GCount_Fiducial_Only1st;
  UInt_t                      GCount_Fiducial_None;

 
  // Histograms 
  TH1F                        *h_N_PV;

  // ID and ISO Efficiencies vs pt, eta, and npv

  // Single muons
  TH1F       *h_Eff_pt_Matched[2];
  TH1F       *h_Eff_pt_TightID[2];
  TH1F       *h_Eff_pt_MediumID[2];
  TH1F       *h_Eff_pt_HWWID[2];
  TH1F       *h_Eff_pt_TightIDGoT[2];
  TH1F       *h_Eff_pt_TightIDipsHWW[2];
  TH1F       *h_Eff_pt_MediumIDipsHWW[2];
  TH1F       *h_Eff_pt_TightID_ISO03[2];
  TH1F       *h_Eff_pt_TightID_ISO04[2];
  TH1F       *h_Eff_pt_TightID_ISO03dBeta[2];
  TH1F       *h_Eff_pt_TightID_ISO04dBeta[2];
  TH1F       *h_Eff_pt_TightID_ISO03PFWeighted[2];
  TH1F       *h_Eff_pt_TightID_ISO04PFWeighted[2];
  TH1F       *h_Eff_pt_TightID_ISO03PUPPI[2];
  TH1F       *h_Eff_pt_TightID_ISO04PUPPI[2];
  TH1F       *h_Eff_pt_MediumID_ISO03[2];
  TH1F       *h_Eff_pt_MediumID_ISO04[2];
  TH1F       *h_Eff_pt_MediumID_ISO03dBeta[2];
  TH1F       *h_Eff_pt_MediumID_ISO04dBeta[2];
  TH1F       *h_Eff_pt_MediumID_ISO03PFWeighted[2];
  TH1F       *h_Eff_pt_MediumID_ISO04PFWeighted[2];
  TH1F       *h_Eff_pt_MediumID_ISO03PUPPI[2];
  TH1F       *h_Eff_pt_MediumID_ISO04PUPPI[2];
  
  TH1F       *h_Eff_eta_Matched[2];
  TH1F       *h_Eff_eta_TightID[2];
  TH1F       *h_Eff_eta_MediumID[2];
  TH1F       *h_Eff_eta_HWWID[2];
  TH1F       *h_Eff_eta_TightIDGoT[2];
  TH1F       *h_Eff_eta_TightIDipsHWW[2];
  TH1F       *h_Eff_eta_MediumIDipsHWW[2];
  TH1F       *h_Eff_eta_TightID_ISO03[2];
  TH1F       *h_Eff_eta_TightID_ISO04[2];
  TH1F       *h_Eff_eta_TightID_ISO03dBeta[2];
  TH1F       *h_Eff_eta_TightID_ISO04dBeta[2];
  TH1F       *h_Eff_eta_TightID_ISO03PFWeighted[2];
  TH1F       *h_Eff_eta_TightID_ISO04PFWeighted[2];
  TH1F       *h_Eff_eta_TightID_ISO03PUPPI[2];
  TH1F       *h_Eff_eta_TightID_ISO04PUPPI[2];
  TH1F       *h_Eff_eta_MediumID_ISO03[2];
  TH1F       *h_Eff_eta_MediumID_ISO04[2];
  TH1F       *h_Eff_eta_MediumID_ISO03dBeta[2];
  TH1F       *h_Eff_eta_MediumID_ISO04dBeta[2];
  TH1F       *h_Eff_eta_MediumID_ISO03PFWeighted[2];
  TH1F       *h_Eff_eta_MediumID_ISO04PFWeighted[2];
  TH1F       *h_Eff_eta_MediumID_ISO03PUPPI[2];
  TH1F       *h_Eff_eta_MediumID_ISO04PUPPI[2];
  
  TH1F       *h_Eff_npv_Matched[2];
  TH1F       *h_Eff_npv_TightID[2];
  TH1F       *h_Eff_npv_MediumID[2];
  TH1F       *h_Eff_npv_HWWID[2];
  TH1F       *h_Eff_npv_TightIDGoT[2];
  TH1F       *h_Eff_npv_TightIDipsHWW[2];
  TH1F       *h_Eff_npv_MediumIDipsHWW[2];
  TH1F       *h_Eff_npv_TightID_ISO03[2];
  TH1F       *h_Eff_npv_TightID_ISO04[2];
  TH1F       *h_Eff_npv_TightID_ISO03dBeta[2];
  TH1F       *h_Eff_npv_TightID_ISO04dBeta[2];
  TH1F       *h_Eff_npv_TightID_ISO03PFWeighted[2];
  TH1F       *h_Eff_npv_TightID_ISO04PFWeighted[2];
  TH1F       *h_Eff_npv_TightID_ISO03PUPPI[2];
  TH1F       *h_Eff_npv_TightID_ISO04PUPPI[2];
  TH1F       *h_Eff_npv_MediumID_ISO03[2];
  TH1F       *h_Eff_npv_MediumID_ISO04[2];
  TH1F       *h_Eff_npv_MediumID_ISO03dBeta[2];
  TH1F       *h_Eff_npv_MediumID_ISO04dBeta[2];
  TH1F       *h_Eff_npv_MediumID_ISO03PFWeighted[2];
  TH1F       *h_Eff_npv_MediumID_ISO04PFWeighted[2];
  TH1F       *h_Eff_npv_MediumID_ISO03PUPPI[2];
  TH1F       *h_Eff_npv_MediumID_ISO04PUPPI[2];

  // Dileptonic ISO efficiencies
  TH1F       *h_Eff_pt_Matched_Dilep;
  TH1F       *h_Eff_pt_HWWID_Dilep;
  TH1F       *h_Eff_pt_TightIDGoT_Dilep;
  TH1F       *h_Eff_pt_TightIDipsHWW_Dilep;
  TH1F       *h_Eff_pt_MediumIDipsHWW_Dilep;
  TH1F       *h_Eff_pt_TightID_Dilep;
  TH1F       *h_Eff_pt_TightID_ISO03_Dilep;
  TH1F       *h_Eff_pt_TightID_ISO04_Dilep;
  TH1F       *h_Eff_pt_TightID_ISO03dBeta_Dilep;
  TH1F       *h_Eff_pt_TightID_ISO04dBeta_Dilep;
  TH1F       *h_Eff_pt_TightID_ISO03PFWeighted_Dilep;
  TH1F       *h_Eff_pt_TightID_ISO04PFWeighted_Dilep;
  TH1F       *h_Eff_pt_TightID_ISO03PUPPI_Dilep;
  TH1F       *h_Eff_pt_TightID_ISO04PUPPI_Dilep;
  TH1F       *h_Eff_pt_MediumID_Dilep;
  TH1F       *h_Eff_pt_MediumID_ISO03_Dilep;
  TH1F       *h_Eff_pt_MediumID_ISO04_Dilep;
  TH1F       *h_Eff_pt_MediumID_ISO03dBeta_Dilep;
  TH1F       *h_Eff_pt_MediumID_ISO04dBeta_Dilep;
  TH1F       *h_Eff_pt_MediumID_ISO03PFWeighted_Dilep;
  TH1F       *h_Eff_pt_MediumID_ISO04PFWeighted_Dilep;
  TH1F       *h_Eff_pt_MediumID_ISO03PUPPI_Dilep;
  TH1F       *h_Eff_pt_MediumID_ISO04PUPPI_Dilep;

  TH1F       *h_Eff_eta_Matched_Dilep;
  TH1F       *h_Eff_eta_HWWID_Dilep;
  TH1F       *h_Eff_eta_TightIDGoT_Dilep;
  TH1F       *h_Eff_eta_TightIDipsHWW_Dilep;
  TH1F       *h_Eff_eta_MediumIDipsHWW_Dilep;
  TH1F       *h_Eff_eta_TightID_Dilep;
  TH1F       *h_Eff_eta_TightID_ISO03_Dilep;
  TH1F       *h_Eff_eta_TightID_ISO04_Dilep;
  TH1F       *h_Eff_eta_TightID_ISO03dBeta_Dilep;
  TH1F       *h_Eff_eta_TightID_ISO04dBeta_Dilep;
  TH1F       *h_Eff_eta_TightID_ISO03PFWeighted_Dilep;
  TH1F       *h_Eff_eta_TightID_ISO04PFWeighted_Dilep;
  TH1F       *h_Eff_eta_TightID_ISO03PUPPI_Dilep;
  TH1F       *h_Eff_eta_TightID_ISO04PUPPI_Dilep;
  TH1F       *h_Eff_eta_MediumID_Dilep;
  TH1F       *h_Eff_eta_MediumID_ISO03_Dilep;
  TH1F       *h_Eff_eta_MediumID_ISO04_Dilep;
  TH1F       *h_Eff_eta_MediumID_ISO03dBeta_Dilep;
  TH1F       *h_Eff_eta_MediumID_ISO04dBeta_Dilep;
  TH1F       *h_Eff_eta_MediumID_ISO03PFWeighted_Dilep;
  TH1F       *h_Eff_eta_MediumID_ISO04PFWeighted_Dilep;
  TH1F       *h_Eff_eta_MediumID_ISO03PUPPI_Dilep;
  TH1F       *h_Eff_eta_MediumID_ISO04PUPPI_Dilep;

  TH1F       *h_Eff_npv_Matched_Dilep;
  TH1F       *h_Eff_npv_HWWID_Dilep;
  TH1F       *h_Eff_npv_TightIDGoT_Dilep;
  TH1F       *h_Eff_npv_TightIDipsHWW_Dilep;
  TH1F       *h_Eff_npv_MediumIDipsHWW_Dilep;
  TH1F       *h_Eff_npv_TightID_Dilep;
  TH1F       *h_Eff_npv_TightID_ISO03_Dilep;
  TH1F       *h_Eff_npv_TightID_ISO04_Dilep;
  TH1F       *h_Eff_npv_TightID_ISO03dBeta_Dilep;
  TH1F       *h_Eff_npv_TightID_ISO04dBeta_Dilep;
  TH1F       *h_Eff_npv_TightID_ISO03PFWeighted_Dilep;
  TH1F       *h_Eff_npv_TightID_ISO04PFWeighted_Dilep;
  TH1F       *h_Eff_npv_TightID_ISO03PUPPI_Dilep;
  TH1F       *h_Eff_npv_TightID_ISO04PUPPI_Dilep;
  TH1F       *h_Eff_npv_MediumID_Dilep;
  TH1F       *h_Eff_npv_MediumID_ISO03_Dilep;
  TH1F       *h_Eff_npv_MediumID_ISO04_Dilep;
  TH1F       *h_Eff_npv_MediumID_ISO03dBeta_Dilep;
  TH1F       *h_Eff_npv_MediumID_ISO04dBeta_Dilep;
  TH1F       *h_Eff_npv_MediumID_ISO03PFWeighted_Dilep;
  TH1F       *h_Eff_npv_MediumID_ISO04PFWeighted_Dilep;
  TH1F       *h_Eff_npv_MediumID_ISO03PUPPI_Dilep;
  TH1F       *h_Eff_npv_MediumID_ISO04PUPPI_Dilep;

  //ISO ROC Curves
  TH1F       *h_RC_TightID_ISO03_Dilep;
  TH1F       *h_RC_TightID_ISO04_Dilep;
  TH1F       *h_RC_TightID_ISO03dBeta_Dilep;
  TH1F       *h_RC_TightID_ISO04dBeta_Dilep;
  TH1F       *h_RC_TightID_ISO03PFWeighted_Dilep;
  TH1F       *h_RC_TightID_ISO04PFWeighted_Dilep;
  TH1F       *h_RC_TightID_ISO03PUPPI_Dilep;
  TH1F       *h_RC_TightID_ISO04PUPPI_Dilep;
  TH1F       *h_RC_MediumID_ISO03_Dilep;
  TH1F       *h_RC_MediumID_ISO04_Dilep;
  TH1F       *h_RC_MediumID_ISO03dBeta_Dilep;
  TH1F       *h_RC_MediumID_ISO04dBeta_Dilep;
  TH1F       *h_RC_MediumID_ISO03PFWeighted_Dilep;
  TH1F       *h_RC_MediumID_ISO04PFWeighted_Dilep;
  TH1F       *h_RC_MediumID_ISO03PUPPI_Dilep;
  TH1F       *h_RC_MediumID_ISO04PUPPI_Dilep;
  
  // Input parameters
  TString                     _Signal;       // Type of Signal
  int                         _NEvents;      // Total number of events in the sample before skim
  float                       _Luminosity;   // Total luminosity
  float                       _XSection;     // Process cross section
  bool                        _IsDATA;       // True if is Data, False in case MC
  int                         _WhichRun;     // 1 in case of RunI samples. 2 In case of RunII samples.
  bool                        _Debug;        // True for verbose while debugging
  bool                        _Report;       // Count events and print final report

  // Weights
  float                       _factN;        // Normalization factor



 public:  
 CoreMuonSelector() : 
     PAFChainItemSelector(),
     G_GEN_isMuMu(),
     G_GEN_isMuTau(),
     G_GEN_isTauMu(),
     G_GEN_isTauTau(),
     G_GEN_Pass(),
     G_GEN_isNonPrompt(),
     G_PassMatching(),
     G_RecoMuSize(),
     G_NPV(),
     EvtFlag_Fiducial(),
     EvtFlag_Gen(),
     EvtFlag_Matching(),
       
     GCount_AllEvents(),
     GCount_GenEvents(),
     GCount_Fiducial_AtLeast2(),
     GCount_Fiducial_2(),
     GCount_Fiducial_1st2nd(),
     GCount_Match_1st2nd(),
     GCount_NoMatch_1st2nd(),
     GCount_MatchTight_1st2nd(),
     GCount_MatchTightIso_1st2nd(),
     GCount_MatchTightIso_Only1st(),
     GCount_MatchTightIso_Only2nd(),
     GCount_MatchTightIso_None(),
     GCount_MatchTight_Only1st(),
     GCount_MatchTight_Only2nd(),
     GCount_MatchTight_None(),
     GCount_Tight_1st2nd(),
     GCount_TightIso_1st2nd(),
     GCount_TightIso_Only1st(),
     GCount_TightIso_Only2nd(),
     GCount_TightIso_None(),
     GCount_Tight_Only1st(),
     GCount_Tight_Only2nd(),
     GCount_Tight_None(),
     GCount_Fiducial_1st3rd(),
     GCount_Match_1st3rd(),
     GCount_NoMatch_1st3rd(),
     GCount_MatchTight_1st3rd(),
     GCount_MatchTightIso_1st3rd(),
     GCount_MatchTightIso_Only3rd(),
     GCount_MatchTight_Only3rd(),
     GCount_Tight_1st3rd(),
     GCount_TightIso_1st3rd(),
     GCount_TightIso_Only3rd(),
     GCount_Tight_Only3rd(),
     GCount_Fiducial_1stOther(),
     GCount_Fiducial_MoreThan2(),
     GCount_Match_MoreThan2_OK(),
     GCount_Match_MoreThan2(),
     GCount_NoMatch_MoreThan2(),
     GCount_Fiducial_3(),
     GCount_Fiducial_MoreThan3(),
     GCount_Fiducial_Only1st(),
     GCount_Fiducial_None(),
     
     h_N_PV(),

     h_Eff_pt_Matched(),
     h_Eff_pt_TightID(),
     h_Eff_pt_MediumID(),
     h_Eff_pt_HWWID(),
     h_Eff_pt_TightIDGoT(),
     h_Eff_pt_TightIDipsHWW(),
     h_Eff_pt_MediumIDipsHWW(),
     h_Eff_pt_TightID_ISO03(),
     h_Eff_pt_TightID_ISO04(),
     h_Eff_pt_TightID_ISO03dBeta(),
     h_Eff_pt_TightID_ISO04dBeta(),
     h_Eff_pt_TightID_ISO03PFWeighted(),
     h_Eff_pt_TightID_ISO04PFWeighted(),
     h_Eff_pt_TightID_ISO03PUPPI(),
     h_Eff_pt_TightID_ISO04PUPPI(),
     h_Eff_pt_MediumID_ISO03(),
     h_Eff_pt_MediumID_ISO04(),
     h_Eff_pt_MediumID_ISO03dBeta(),
     h_Eff_pt_MediumID_ISO04dBeta(),
     h_Eff_pt_MediumID_ISO03PFWeighted(),
     h_Eff_pt_MediumID_ISO04PFWeighted(),
     h_Eff_pt_MediumID_ISO03PUPPI(),
     h_Eff_pt_MediumID_ISO04PUPPI(),

     h_Eff_eta_Matched(),
     h_Eff_eta_TightID(),
     h_Eff_eta_MediumID(),
     h_Eff_eta_HWWID(),
     h_Eff_eta_TightIDGoT(),
     h_Eff_eta_TightIDipsHWW(),
     h_Eff_eta_MediumIDipsHWW(),
     h_Eff_eta_TightID_ISO03(),
     h_Eff_eta_TightID_ISO04(),
     h_Eff_eta_TightID_ISO03dBeta(),
     h_Eff_eta_TightID_ISO04dBeta(),
     h_Eff_eta_TightID_ISO03PFWeighted(),
     h_Eff_eta_TightID_ISO04PFWeighted(),
     h_Eff_eta_TightID_ISO03PUPPI(),
     h_Eff_eta_TightID_ISO04PUPPI(),
     h_Eff_eta_MediumID_ISO03(),
     h_Eff_eta_MediumID_ISO04(),
     h_Eff_eta_MediumID_ISO03dBeta(),
     h_Eff_eta_MediumID_ISO04dBeta(),
     h_Eff_eta_MediumID_ISO03PFWeighted(),
     h_Eff_eta_MediumID_ISO04PFWeighted(),
     h_Eff_eta_MediumID_ISO03PUPPI(),
     h_Eff_eta_MediumID_ISO04PUPPI(),
  
     h_Eff_npv_Matched(),
     h_Eff_npv_TightID(),
     h_Eff_npv_MediumID(),
     h_Eff_npv_HWWID(),
     h_Eff_npv_TightIDGoT(),
     h_Eff_npv_TightIDipsHWW(),
     h_Eff_npv_MediumIDipsHWW(),
     h_Eff_npv_TightID_ISO03(),
     h_Eff_npv_TightID_ISO04(),
     h_Eff_npv_TightID_ISO03dBeta(),
     h_Eff_npv_TightID_ISO04dBeta(),
     h_Eff_npv_TightID_ISO03PFWeighted(),
     h_Eff_npv_TightID_ISO04PFWeighted(),
     h_Eff_npv_TightID_ISO03PUPPI(),
     h_Eff_npv_TightID_ISO04PUPPI(),
     h_Eff_npv_MediumID_ISO03(),
     h_Eff_npv_MediumID_ISO04(),
     h_Eff_npv_MediumID_ISO03dBeta(),
     h_Eff_npv_MediumID_ISO04dBeta(),
     h_Eff_npv_MediumID_ISO03PFWeighted(),
     h_Eff_npv_MediumID_ISO04PFWeighted(),
     h_Eff_npv_MediumID_ISO03PUPPI(),
     h_Eff_npv_MediumID_ISO04PUPPI(),

     h_Eff_pt_Matched_Dilep(),
     h_Eff_pt_HWWID_Dilep(),
     h_Eff_pt_TightIDGoT_Dilep(),
     h_Eff_pt_TightIDipsHWW_Dilep(),
     h_Eff_pt_MediumIDipsHWW_Dilep(),
     h_Eff_pt_TightID_Dilep(),
     h_Eff_pt_TightID_ISO03_Dilep(),
     h_Eff_pt_TightID_ISO04_Dilep(),
     h_Eff_pt_TightID_ISO03dBeta_Dilep(),
     h_Eff_pt_TightID_ISO04dBeta_Dilep(),
     h_Eff_pt_TightID_ISO03PFWeighted_Dilep(),
     h_Eff_pt_TightID_ISO04PFWeighted_Dilep(),
     h_Eff_pt_TightID_ISO03PUPPI_Dilep(),
     h_Eff_pt_TightID_ISO04PUPPI_Dilep(),
     h_Eff_pt_MediumID_Dilep(),
     h_Eff_pt_MediumID_ISO03_Dilep(),
     h_Eff_pt_MediumID_ISO04_Dilep(),
     h_Eff_pt_MediumID_ISO03dBeta_Dilep(),
     h_Eff_pt_MediumID_ISO04dBeta_Dilep(),
     h_Eff_pt_MediumID_ISO03PFWeighted_Dilep(),
     h_Eff_pt_MediumID_ISO04PFWeighted_Dilep(),
     h_Eff_pt_MediumID_ISO03PUPPI_Dilep(),
     h_Eff_pt_MediumID_ISO04PUPPI_Dilep(),

     h_Eff_eta_Matched_Dilep(),
     h_Eff_eta_HWWID_Dilep(),
     h_Eff_eta_TightIDGoT_Dilep(),
     h_Eff_eta_TightIDipsHWW_Dilep(),
     h_Eff_eta_MediumIDipsHWW_Dilep(),
     h_Eff_eta_TightID_Dilep(),
     h_Eff_eta_TightID_ISO03_Dilep(),
     h_Eff_eta_TightID_ISO04_Dilep(),
     h_Eff_eta_TightID_ISO03dBeta_Dilep(),
     h_Eff_eta_TightID_ISO04dBeta_Dilep(),
     h_Eff_eta_TightID_ISO03PFWeighted_Dilep(),
     h_Eff_eta_TightID_ISO04PFWeighted_Dilep(),
     h_Eff_eta_TightID_ISO03PUPPI_Dilep(),
     h_Eff_eta_TightID_ISO04PUPPI_Dilep(),
     h_Eff_eta_MediumID_Dilep(),
     h_Eff_eta_MediumID_ISO03_Dilep(),
     h_Eff_eta_MediumID_ISO04_Dilep(),
     h_Eff_eta_MediumID_ISO03dBeta_Dilep(),
     h_Eff_eta_MediumID_ISO04dBeta_Dilep(),
     h_Eff_eta_MediumID_ISO03PFWeighted_Dilep(),
     h_Eff_eta_MediumID_ISO04PFWeighted_Dilep(),
     h_Eff_eta_MediumID_ISO03PUPPI_Dilep(),
     h_Eff_eta_MediumID_ISO04PUPPI_Dilep(),

     h_Eff_npv_Matched_Dilep(),
     h_Eff_npv_HWWID_Dilep(),
     h_Eff_npv_TightIDGoT_Dilep(),
     h_Eff_npv_TightIDipsHWW_Dilep(),
     h_Eff_npv_MediumIDipsHWW_Dilep(),
     h_Eff_npv_TightID_Dilep(),
     h_Eff_npv_TightID_ISO03_Dilep(),
     h_Eff_npv_TightID_ISO04_Dilep(),
     h_Eff_npv_TightID_ISO03dBeta_Dilep(),
     h_Eff_npv_TightID_ISO04dBeta_Dilep(),
     h_Eff_npv_TightID_ISO03PFWeighted_Dilep(),
     h_Eff_npv_TightID_ISO04PFWeighted_Dilep(),
     h_Eff_npv_TightID_ISO03PUPPI_Dilep(),
     h_Eff_npv_TightID_ISO04PUPPI_Dilep(),
     h_Eff_npv_MediumID_Dilep(),
     h_Eff_npv_MediumID_ISO03_Dilep(),
     h_Eff_npv_MediumID_ISO04_Dilep(),
     h_Eff_npv_MediumID_ISO03dBeta_Dilep(),
     h_Eff_npv_MediumID_ISO04dBeta_Dilep(),
     h_Eff_npv_MediumID_ISO03PFWeighted_Dilep(),
     h_Eff_npv_MediumID_ISO04PFWeighted_Dilep(),
     h_Eff_npv_MediumID_ISO03PUPPI_Dilep(),
     h_Eff_npv_MediumID_ISO04PUPPI_Dilep(),

     h_RC_TightID_ISO03_Dilep(),
     h_RC_TightID_ISO04_Dilep(),
     h_RC_TightID_ISO03dBeta_Dilep(),
     h_RC_TightID_ISO04dBeta_Dilep(),
     h_RC_TightID_ISO03PFWeighted_Dilep(),
     h_RC_TightID_ISO04PFWeighted_Dilep(),
     h_RC_TightID_ISO03PUPPI_Dilep(),
     h_RC_TightID_ISO04PUPPI_Dilep(),
     h_RC_MediumID_ISO03_Dilep(),
     h_RC_MediumID_ISO04_Dilep(),
     h_RC_MediumID_ISO03dBeta_Dilep(),
     h_RC_MediumID_ISO04dBeta_Dilep(),
     h_RC_MediumID_ISO03PFWeighted_Dilep(),
     h_RC_MediumID_ISO04PFWeighted_Dilep(),
     h_RC_MediumID_ISO03PUPPI_Dilep(),
     h_RC_MediumID_ISO04PUPPI_Dilep(),
     
     _Signal(),
     _NEvents(),
     _Luminosity(),
     _XSection(),
     _IsDATA(),
     _WhichRun(),
     _Debug(),
     _Report(),
     _factN()
       { }

   ClassDef(CoreMuonSelector,0);
};
