///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////                                                                                             /////////////
/////////////                                  MUON EFFICIENCIES SELECTOR                                 /////////////
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


class MuonEfficienciesSelector: public PAFChainItemSelector{
  
 public:
  virtual ~MuonEfficienciesSelector() {}
  
  virtual void                Initialise();
  virtual void                InsideLoop();
  virtual void                Summary();
  
 protected:
  // See description in implementation
  void                        EffsSingleMu(UInt_t);
  void                        EffsAllMu();


  // My Declarations:
  // Define data members

  // VARIABLES FOR EACH EVENT (to be initialized after every event)
 
  // GEN Info
  std::vector<TLorentzVector> G_GEN_PromptMuon_4vec;  //Lorentz vector for all GEN prompt muons, or coming from 
                                                      //a prompt tau decay, ordered by Pt
  bool                        G_GEN_isMuMu;      //There are 2 GEN prompt muons that come directly from a boson  
  bool                        G_GEN_isMuTau;     // "                        "  but the 2nd (in Pt) comes from a tau
  bool                        G_GEN_isTauMu;     // "                        "  but the 1st comes from a tau
  bool                        G_GEN_isTauTau;    // "                        "  and both come from a tau
  
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

  // Sizes
  UInt_t                      G_RecoMuSize;      //How many RECO muons there are
  UInt_t                      G_NPV;             //How many Primary Vtx. there are

  // Event Flags
  bool                        EvtFlag_Fiducial; // 1st and 2nd muon pass Fiducial Selection
  bool                        EvtFlag_Gen;      // Same as G_GEN_Pass
  bool                        EvtFlag_Matching; // Same as G_PassMatching


  // VARIABLES FOR ALL EVENTS (to be initialized only once)

  // Histograms 

  // Counting histos for ID and ISO Efficiencies vs pt, eta, and npv

  // Single muons
  TH1F       *h_Eff_pt_NoID[2];
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
  
  TH1F       *h_Eff_eta_NoID[2];
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
  
  TH1F       *h_Eff_npv_NoID[2];
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

  // ISO efficiencies for All (matched) muons (in case of MC like DY, WW, TTbar...)
  TH1F       *h_Eff_pt_NoID_AllMu;
  TH1F       *h_Eff_pt_HWWID_AllMu;
  TH1F       *h_Eff_pt_TightIDGoT_AllMu;
  TH1F       *h_Eff_pt_TightIDipsHWW_AllMu;
  TH1F       *h_Eff_pt_MediumIDipsHWW_AllMu;
  TH1F       *h_Eff_pt_TightID_AllMu;
  TH1F       *h_Eff_pt_TightID_ISO03_AllMu;
  TH1F       *h_Eff_pt_TightID_ISO04_AllMu;
  TH1F       *h_Eff_pt_TightID_ISO03dBeta_AllMu;
  TH1F       *h_Eff_pt_TightID_ISO04dBeta_AllMu;
  TH1F       *h_Eff_pt_TightID_ISO03PFWeighted_AllMu;
  TH1F       *h_Eff_pt_TightID_ISO04PFWeighted_AllMu;
  TH1F       *h_Eff_pt_TightID_ISO03PUPPI_AllMu;
  TH1F       *h_Eff_pt_TightID_ISO04PUPPI_AllMu;
  TH1F       *h_Eff_pt_MediumID_AllMu;
  TH1F       *h_Eff_pt_MediumID_ISO03_AllMu;
  TH1F       *h_Eff_pt_MediumID_ISO04_AllMu;
  TH1F       *h_Eff_pt_MediumID_ISO03dBeta_AllMu;
  TH1F       *h_Eff_pt_MediumID_ISO04dBeta_AllMu;
  TH1F       *h_Eff_pt_MediumID_ISO03PFWeighted_AllMu;
  TH1F       *h_Eff_pt_MediumID_ISO04PFWeighted_AllMu;
  TH1F       *h_Eff_pt_MediumID_ISO03PUPPI_AllMu;
  TH1F       *h_Eff_pt_MediumID_ISO04PUPPI_AllMu;

  TH1F       *h_Eff_eta_NoID_AllMu;
  TH1F       *h_Eff_eta_HWWID_AllMu;
  TH1F       *h_Eff_eta_TightIDGoT_AllMu;
  TH1F       *h_Eff_eta_TightIDipsHWW_AllMu;
  TH1F       *h_Eff_eta_MediumIDipsHWW_AllMu;
  TH1F       *h_Eff_eta_TightID_AllMu;
  TH1F       *h_Eff_eta_TightID_ISO03_AllMu;
  TH1F       *h_Eff_eta_TightID_ISO04_AllMu;
  TH1F       *h_Eff_eta_TightID_ISO03dBeta_AllMu;
  TH1F       *h_Eff_eta_TightID_ISO04dBeta_AllMu;
  TH1F       *h_Eff_eta_TightID_ISO03PFWeighted_AllMu;
  TH1F       *h_Eff_eta_TightID_ISO04PFWeighted_AllMu;
  TH1F       *h_Eff_eta_TightID_ISO03PUPPI_AllMu;
  TH1F       *h_Eff_eta_TightID_ISO04PUPPI_AllMu;
  TH1F       *h_Eff_eta_MediumID_AllMu;
  TH1F       *h_Eff_eta_MediumID_ISO03_AllMu;
  TH1F       *h_Eff_eta_MediumID_ISO04_AllMu;
  TH1F       *h_Eff_eta_MediumID_ISO03dBeta_AllMu;
  TH1F       *h_Eff_eta_MediumID_ISO04dBeta_AllMu;
  TH1F       *h_Eff_eta_MediumID_ISO03PFWeighted_AllMu;
  TH1F       *h_Eff_eta_MediumID_ISO04PFWeighted_AllMu;
  TH1F       *h_Eff_eta_MediumID_ISO03PUPPI_AllMu;
  TH1F       *h_Eff_eta_MediumID_ISO04PUPPI_AllMu;

  TH1F       *h_Eff_npv_NoID_AllMu;
  TH1F       *h_Eff_npv_HWWID_AllMu;
  TH1F       *h_Eff_npv_TightIDGoT_AllMu;
  TH1F       *h_Eff_npv_TightIDipsHWW_AllMu;
  TH1F       *h_Eff_npv_MediumIDipsHWW_AllMu;
  TH1F       *h_Eff_npv_TightID_AllMu;
  TH1F       *h_Eff_npv_TightID_ISO03_AllMu;
  TH1F       *h_Eff_npv_TightID_ISO04_AllMu;
  TH1F       *h_Eff_npv_TightID_ISO03dBeta_AllMu;
  TH1F       *h_Eff_npv_TightID_ISO04dBeta_AllMu;
  TH1F       *h_Eff_npv_TightID_ISO03PFWeighted_AllMu;
  TH1F       *h_Eff_npv_TightID_ISO04PFWeighted_AllMu;
  TH1F       *h_Eff_npv_TightID_ISO03PUPPI_AllMu;
  TH1F       *h_Eff_npv_TightID_ISO04PUPPI_AllMu;
  TH1F       *h_Eff_npv_MediumID_AllMu;
  TH1F       *h_Eff_npv_MediumID_ISO03_AllMu;
  TH1F       *h_Eff_npv_MediumID_ISO04_AllMu;
  TH1F       *h_Eff_npv_MediumID_ISO03dBeta_AllMu;
  TH1F       *h_Eff_npv_MediumID_ISO04dBeta_AllMu;
  TH1F       *h_Eff_npv_MediumID_ISO03PFWeighted_AllMu;
  TH1F       *h_Eff_npv_MediumID_ISO04PFWeighted_AllMu;
  TH1F       *h_Eff_npv_MediumID_ISO03PUPPI_AllMu;
  TH1F       *h_Eff_npv_MediumID_ISO04PUPPI_AllMu;

  
  // Input parameters
  TString                     _Signal;       // Type of Signal
  int                         _NEvents;      // Total number of events in the sample before skim
  bool                        _IsDATA;       // True if is Data, False in case MC
  int                         _WhichRun;     // 1 in case of RunI samples. 2 In case of RunII samples.
  bool                        _Debug;        // True for verbose while debugging


 public:  
 MuonEfficienciesSelector() : 
     PAFChainItemSelector(),
     G_GEN_isMuMu(),
     G_GEN_isMuTau(),
     G_GEN_isTauMu(),
     G_GEN_isTauTau(),
     G_RecoMuSize(),
     G_NPV(),
     EvtFlag_Fiducial(),
     EvtFlag_Gen(),
     EvtFlag_Matching(),

     h_Eff_pt_NoID(),
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

     h_Eff_eta_NoID(),
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
  
     h_Eff_npv_NoID(),
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

     h_Eff_pt_NoID_AllMu(),
     h_Eff_pt_HWWID_AllMu(),
     h_Eff_pt_TightIDGoT_AllMu(),
     h_Eff_pt_TightIDipsHWW_AllMu(),
     h_Eff_pt_MediumIDipsHWW_AllMu(),
     h_Eff_pt_TightID_AllMu(),
     h_Eff_pt_TightID_ISO03_AllMu(),
     h_Eff_pt_TightID_ISO04_AllMu(),
     h_Eff_pt_TightID_ISO03dBeta_AllMu(),
     h_Eff_pt_TightID_ISO04dBeta_AllMu(),
     h_Eff_pt_TightID_ISO03PFWeighted_AllMu(),
     h_Eff_pt_TightID_ISO04PFWeighted_AllMu(),
     h_Eff_pt_TightID_ISO03PUPPI_AllMu(),
     h_Eff_pt_TightID_ISO04PUPPI_AllMu(),
     h_Eff_pt_MediumID_AllMu(),
     h_Eff_pt_MediumID_ISO03_AllMu(),
     h_Eff_pt_MediumID_ISO04_AllMu(),
     h_Eff_pt_MediumID_ISO03dBeta_AllMu(),
     h_Eff_pt_MediumID_ISO04dBeta_AllMu(),
     h_Eff_pt_MediumID_ISO03PFWeighted_AllMu(),
     h_Eff_pt_MediumID_ISO04PFWeighted_AllMu(),
     h_Eff_pt_MediumID_ISO03PUPPI_AllMu(),
     h_Eff_pt_MediumID_ISO04PUPPI_AllMu(),

     h_Eff_eta_NoID_AllMu(),
     h_Eff_eta_HWWID_AllMu(),
     h_Eff_eta_TightIDGoT_AllMu(),
     h_Eff_eta_TightIDipsHWW_AllMu(),
     h_Eff_eta_MediumIDipsHWW_AllMu(),
     h_Eff_eta_TightID_AllMu(),
     h_Eff_eta_TightID_ISO03_AllMu(),
     h_Eff_eta_TightID_ISO04_AllMu(),
     h_Eff_eta_TightID_ISO03dBeta_AllMu(),
     h_Eff_eta_TightID_ISO04dBeta_AllMu(),
     h_Eff_eta_TightID_ISO03PFWeighted_AllMu(),
     h_Eff_eta_TightID_ISO04PFWeighted_AllMu(),
     h_Eff_eta_TightID_ISO03PUPPI_AllMu(),
     h_Eff_eta_TightID_ISO04PUPPI_AllMu(),
     h_Eff_eta_MediumID_AllMu(),
     h_Eff_eta_MediumID_ISO03_AllMu(),
     h_Eff_eta_MediumID_ISO04_AllMu(),
     h_Eff_eta_MediumID_ISO03dBeta_AllMu(),
     h_Eff_eta_MediumID_ISO04dBeta_AllMu(),
     h_Eff_eta_MediumID_ISO03PFWeighted_AllMu(),
     h_Eff_eta_MediumID_ISO04PFWeighted_AllMu(),
     h_Eff_eta_MediumID_ISO03PUPPI_AllMu(),
     h_Eff_eta_MediumID_ISO04PUPPI_AllMu(),

     h_Eff_npv_NoID_AllMu(),
     h_Eff_npv_HWWID_AllMu(),
     h_Eff_npv_TightIDGoT_AllMu(),
     h_Eff_npv_TightIDipsHWW_AllMu(),
     h_Eff_npv_MediumIDipsHWW_AllMu(),
     h_Eff_npv_TightID_AllMu(),
     h_Eff_npv_TightID_ISO03_AllMu(),
     h_Eff_npv_TightID_ISO04_AllMu(),
     h_Eff_npv_TightID_ISO03dBeta_AllMu(),
     h_Eff_npv_TightID_ISO04dBeta_AllMu(),
     h_Eff_npv_TightID_ISO03PFWeighted_AllMu(),
     h_Eff_npv_TightID_ISO04PFWeighted_AllMu(),
     h_Eff_npv_TightID_ISO03PUPPI_AllMu(),
     h_Eff_npv_TightID_ISO04PUPPI_AllMu(),
     h_Eff_npv_MediumID_AllMu(),
     h_Eff_npv_MediumID_ISO03_AllMu(),
     h_Eff_npv_MediumID_ISO04_AllMu(),
     h_Eff_npv_MediumID_ISO03dBeta_AllMu(),
     h_Eff_npv_MediumID_ISO04dBeta_AllMu(),
     h_Eff_npv_MediumID_ISO03PFWeighted_AllMu(),
     h_Eff_npv_MediumID_ISO04PFWeighted_AllMu(),
     h_Eff_npv_MediumID_ISO03PUPPI_AllMu(),
     h_Eff_npv_MediumID_ISO04PUPPI_AllMu(),
     
     _Signal(),
     _NEvents(),
     _IsDATA(),
     _WhichRun(),
     _Debug()
       { }

   ClassDef(MuonEfficienciesSelector,0);
};
