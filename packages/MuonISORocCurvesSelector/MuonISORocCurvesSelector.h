///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////                                                                                             /////////////
/////////////                                 MUON ISO ROC CURVES SELECTOR                                /////////////
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


class MuonISORocCurvesSelector: public PAFChainItemSelector{
  
 public:
  virtual ~MuonISORocCurvesSelector() {}
  
  virtual void                Initialise();
  virtual void                InsideLoop();
  virtual void                Summary();
  
 protected:
  // See description in implementation
  bool                        passISO(UInt_t, TString, float);
  float                       getISO(UInt_t, TString);
  void                        ISORocCurve();


  // My Declarations:
  // Define data members

  // VARIABLES FOR EACH EVENT (to be initialized after every event)
  
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

  //ISO ROC Curves
  TH1F       *h_RC_TightID_ISO03_AllMu;
  TH1F       *h_RC_TightID_ISO04_AllMu;
  TH1F       *h_RC_TightID_ISO03dBeta_AllMu;
  TH1F       *h_RC_TightID_ISO04dBeta_AllMu;
  TH1F       *h_RC_TightID_ISO03PFWeighted_AllMu;
  TH1F       *h_RC_TightID_ISO04PFWeighted_AllMu;
  TH1F       *h_RC_TightID_ISO03PUPPI_AllMu;
  TH1F       *h_RC_TightID_ISO04PUPPI_AllMu;
  TH1F       *h_RC_MediumID_ISO03_AllMu;
  TH1F       *h_RC_MediumID_ISO04_AllMu;
  TH1F       *h_RC_MediumID_ISO03dBeta_AllMu;
  TH1F       *h_RC_MediumID_ISO04dBeta_AllMu;
  TH1F       *h_RC_MediumID_ISO03PFWeighted_AllMu;
  TH1F       *h_RC_MediumID_ISO04PFWeighted_AllMu;
  TH1F       *h_RC_MediumID_ISO03PUPPI_AllMu;
  TH1F       *h_RC_MediumID_ISO04PUPPI_AllMu;
  
  // Input parameters
  TString                     _Signal;       // Type of Signal
  int                         _NEvents;      // Total number of events in the sample before skim
  bool                        _IsDATA;       // True if is Data, False in case MC
  int                         _WhichRun;     // 1 in case of RunI samples. 2 In case of RunII samples.
  bool                        _Debug;        // True for verbose while debugging


 public:  
 MuonISORocCurvesSelector() : 
     PAFChainItemSelector(),
     G_RecoMuSize(),
     G_NPV(),
     EvtFlag_Fiducial(),
     EvtFlag_Gen(),
     EvtFlag_Matching(),

     h_RC_TightID_ISO03_AllMu(),
     h_RC_TightID_ISO04_AllMu(),
     h_RC_TightID_ISO03dBeta_AllMu(),
     h_RC_TightID_ISO04dBeta_AllMu(),
     h_RC_TightID_ISO03PFWeighted_AllMu(),
     h_RC_TightID_ISO04PFWeighted_AllMu(),
     h_RC_TightID_ISO03PUPPI_AllMu(),
     h_RC_TightID_ISO04PUPPI_AllMu(),
     h_RC_MediumID_ISO03_AllMu(),
     h_RC_MediumID_ISO04_AllMu(),
     h_RC_MediumID_ISO03dBeta_AllMu(),
     h_RC_MediumID_ISO04dBeta_AllMu(),
     h_RC_MediumID_ISO03PFWeighted_AllMu(),
     h_RC_MediumID_ISO04PFWeighted_AllMu(),
     h_RC_MediumID_ISO03PUPPI_AllMu(),
     h_RC_MediumID_ISO04PUPPI_AllMu(),
     
     _Signal(),
     _NEvents(),
     _IsDATA(),
     _WhichRun(),
     _Debug()
       { }

   ClassDef(MuonISORocCurvesSelector,0);
};
