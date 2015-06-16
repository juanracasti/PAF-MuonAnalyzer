////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

/////            Muon Studies for CSA14                     ////

/////        A. Calderón (IFCA)   18 / 08 / 2014            ////

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////



#pragma once

#include "PAF/computing/PAFChainItemSelector.h"

#include <TH1F.h>
#include <TMatrix.h>
#include <TH2F.h>
#include <TLorentzVector.h>
#include "Riostream.h"  


class BasicmuonAnalyzer: public PAFChainItemSelector{
  
 public:
  virtual ~BasicmuonAnalyzer() {}
  
  virtual void                Initialise();
  virtual void                InsideLoop();
  virtual void                Summary();
  
 protected:
  void                        GetAllMuons();
  bool                        passMediumID(int);
  void                        SetGenInfo();
  bool                        MatchGenToReco(int&, int&);
  bool                        passPFIso (int, string, float);
  float                       getPFRelIso (int, string);

  // My Declarations:
  // Define global variables
  
 
  ///** GEN INFORMATION
  std::vector<TLorentzVector> G_GEN_PromptMuon_4vec;
  std::vector<TLorentzVector> G_GEN_Muon_4vec;
  
  //** MUONS
  std::vector<TLorentzVector> G_Muon_4vec;

  std::vector<bool>           G_MuonID_Tight;
  std::vector<bool>           G_MuonID_Medium;
  std::vector<bool>           G_MuonID_HWW;
  std::vector<bool>           G_MuonID_MediumHWW;
  std::vector<bool>           G_MuonID_GLBorTRKArb;
  std::vector<bool>           G_MuonID_Fiducial;
  
  std::vector<bool>           G_MuonISO03;
  std::vector<bool>           G_MuonISO03_dBeta;
  std::vector<bool>           G_MuonISO03_PFWeights;
  std::vector<bool>           G_MuonISO04;
  std::vector<bool>           G_MuonISO04_dBeta;
  std::vector<bool>           G_MuonISO04_PFWeights;
   
  bool                        G_isMuMu;
  bool                        G_isMuTau;
  bool                        G_isTauMu;
  bool                        G_isTauTau;
  bool                        G_isNonPrompt;
   
  bool                        G_Debug_DefineAnalysisVariables;

 
   //* Histograms 
  TH1F                        *h_N_PV;
  TH1F                        *h_N_PV2;
  TH1F                        *h_N_PV3;
  
  
  // * Input parameters
  TString                     _Signal;       // Type of Signal
  int                         _NEvents;      // Total number of events in the sample before skim
  float                       _Luminosity;   // Total luminosity
  float                       _XSection;     // Process cross section
  bool                        _IsDATA; // True if is Data, False in case MC
  int                         _WhichRun; // 1 in case of RunI samples. 2 In case of RunII samples.;

 public:  
 BasicmuonAnalyzer() : 
     PAFChainItemSelector(),
     G_isMuMu(),
     G_isMuTau(),
     G_isTauMu(),
     G_isTauTau(),
     G_isNonPrompt(),
     G_Debug_DefineAnalysisVariables(),
     
     h_N_PV(),
     h_N_PV2(),
     h_N_PV3(),
     
     _Signal(),
     _NEvents(),
     _Luminosity(),
     _XSection(),
     _IsDATA(),
     _WhichRun()
       { }

   ClassDef(BasicmuonAnalyzer,0);
};
