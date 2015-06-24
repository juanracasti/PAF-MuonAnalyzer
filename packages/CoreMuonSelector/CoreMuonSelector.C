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

#include "CoreMuonSelector.h"

#include "TH1.h"
#include "TH2.h"
#include "TKey.h"
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TRandom.h"
#include "TLorentzVector.h"
#include <vector>
#include "TROOT.h"
#include <iostream>

#include "TDatabasePDG.h"

ClassImp(CoreMuonSelector)

// Initialise input parameters and data members for all events
void CoreMuonSelector::Initialise() {

  _Signal     = GetParam<TString>("Signal");
  _IsDATA     = GetParam<bool>("IsDATA");
  _NEvents    = GetParam<int>("NEvents");
  _Luminosity = GetParam<float>("Luminosity");
  _XSection   = GetParam<float>("XSection");
  _WhichRun   = GetParam<int>("WhichRun"); 
  _Debug      = GetParam<bool>("Debug");
  _Report     = GetParam<bool>("Report");

  //Define weights
  _factN = 1.;
  if (!_IsDATA && _XSection > 0) _factN = _XSection * _Luminosity / _NEvents;

  //For counting
  if (_Report) {
    GCount_AllEvents = 0;
    GCount_GenEvents = 0;
    GCount_Fiducial_AtLeast2 = 0;
    GCount_Fiducial_2 = 0;
    GCount_Fiducial_1st2nd = 0;
    GCount_Match_1st2nd = 0;
    GCount_NoMatch_1st2nd = 0;
    GCount_MatchTight_1st2nd = 0;
    GCount_MatchTightIso_1st2nd = 0;
    GCount_MatchTightIso_Only1st = 0;
    GCount_MatchTightIso_Only2nd = 0;
    GCount_MatchTightIso_None = 0;
    GCount_MatchTight_Only1st = 0;
    GCount_MatchTight_Only2nd = 0;
    GCount_MatchTight_None = 0;
    GCount_Tight_1st2nd = 0;
    GCount_TightIso_1st2nd = 0;
    GCount_TightIso_Only1st = 0;
    GCount_TightIso_Only2nd = 0;
    GCount_TightIso_None = 0;
    GCount_Tight_Only1st = 0;
    GCount_Tight_Only2nd = 0;
    GCount_Tight_None = 0;
    GCount_Fiducial_1st3rd = 0;
    GCount_Match_1st3rd = 0;
    GCount_NoMatch_1st3rd = 0;
    GCount_MatchTight_1st3rd = 0;
    GCount_MatchTightIso_1st3rd = 0;
    GCount_MatchTightIso_Only3rd = 0;
    GCount_MatchTight_Only3rd = 0;
    GCount_Tight_1st3rd = 0;
    GCount_TightIso_1st3rd = 0;
    GCount_TightIso_Only3rd = 0;
    GCount_Tight_Only3rd = 0;
    GCount_Fiducial_1stOther = 0;
    GCount_Fiducial_MoreThan2 = 0;
    GCount_Match_MoreThan2_OK = 0;
    GCount_Match_MoreThan2 = 0;
    GCount_NoMatch_MoreThan2 = 0;
    GCount_Fiducial_3 = 0;
    GCount_Fiducial_MoreThan3 = 0;
    GCount_Fiducial_Only1st = 0;
    GCount_Fiducial_None = 0;
  }

 
//------------------------------------------------------------------------------
// Create histos
//------------------------------------------------------------------------------
 

  h_N_PV  = CreateH1F ("h_N_PV","h_N_PV",50,0,50); 
  h_N_PV->TH1::SetDefaultSumw2();

  // Single Muon ID and ISO efficiencies vs pt, eta and npv

  h_Eff_pt_Matched[0]                   = CreateH1F("h_Eff_pt_Matched_Mu1",
						    "h_Eff_pt_Matched_Mu1",                    200, 0, 200);
  h_Eff_pt_Matched[1]                   = CreateH1F("h_Eff_pt_Matched_Mu2",
						    "h_Eff_pt_Matched_Mu2",                    200, 0, 200);
  h_Eff_pt_TightID[0]                   = CreateH1F("h_Eff_pt_TightID_Mu1",
						    "h_Eff_pt_TightID_Mu1",                    200, 0, 200);
  h_Eff_pt_TightID[1]                   = CreateH1F("h_Eff_pt_TightID_Mu2",
						    "h_Eff_pt_TightID_Mu2",                    200, 0, 200);
  h_Eff_pt_MediumID[0]                  = CreateH1F("h_Eff_pt_MediumID_Mu1",
						    "h_Eff_pt_MediumID_Mu1",                   200, 0, 200);
  h_Eff_pt_MediumID[1]                  = CreateH1F("h_Eff_pt_MediumID_Mu2",
						    "h_Eff_pt_MediumID_Mu2",                   200, 0, 200);
  h_Eff_pt_HWWID[0]                     = CreateH1F("h_Eff_pt_HWWID_Mu1",
						    "h_Eff_pt_HWWID_Mu1",                      200, 0, 200);
  h_Eff_pt_HWWID[1]                     = CreateH1F("h_Eff_pt_HWWID_Mu2",
						    "h_Eff_pt_HWWID_Mu2",                      200, 0, 200);
  h_Eff_pt_TightIDGoT[0]                = CreateH1F("h_Eff_pt_TightIDGoT_Mu1",
						    "h_Eff_pt_TightIDGoT_Mu1",                 200, 0, 200);
  h_Eff_pt_TightIDGoT[1]                = CreateH1F("h_Eff_pt_TightIDGoT_Mu2",
						    "h_Eff_pt_TightIDGoT_Mu2",                 200, 0, 200);
  h_Eff_pt_TightIDipsHWW[0]             = CreateH1F("h_Eff_pt_TightIDipsHWW_Mu1",
						    "h_Eff_pt_TightIDipsHWW_Mu1",              200, 0, 200);
  h_Eff_pt_TightIDipsHWW[1]             = CreateH1F("h_Eff_pt_TightIDipsHWW_Mu2",
						    "h_Eff_pt_TightIDipsHWW_Mu2",              200, 0, 200);
  h_Eff_pt_MediumIDipsHWW[0]            = CreateH1F("h_Eff_pt_MediumIDipsHWW_Mu1",
						    "h_Eff_pt_MediumIDipsHWW_Mu1",             200, 0, 200);
  h_Eff_pt_MediumIDipsHWW[1]            = CreateH1F("h_Eff_pt_MediumIDipsHWW_Mu2",
						    "h_Eff_pt_MediumIDipsHWW_Mu2",             200, 0, 200);
  h_Eff_pt_TightID_ISO03[0]             = CreateH1F("h_Eff_pt_TightID_ISO03_Mu1", 
						    "h_Eff_pt_TightID_ISO03_Mu1",              200, 0, 200);
  h_Eff_pt_TightID_ISO03[1]             = CreateH1F("h_Eff_pt_TightID_ISO03_Mu2", 
						    "h_Eff_pt_TightID_ISO03_Mu2",              200, 0, 200);
  h_Eff_pt_TightID_ISO04[0]             = CreateH1F("h_Eff_pt_TightID_ISO04_Mu1", 
						    "h_Eff_pt_TightID_ISO04_Mu1",              200, 0, 200);
  h_Eff_pt_TightID_ISO04[1]             = CreateH1F("h_Eff_pt_TightID_ISO04_Mu2", 
						    "h_Eff_pt_TightID_ISO04_Mu2",              200, 0, 200);
  h_Eff_pt_TightID_ISO03dBeta[0]        = CreateH1F("h_Eff_pt_TightID_ISO03dBeta_Mu1", 
						    "h_Eff_pt_TightID_ISO03dBeta_Mu1",         200, 0, 200);
  h_Eff_pt_TightID_ISO03dBeta[1]        = CreateH1F("h_Eff_pt_TightID_ISO03dBeta_Mu2", 
						    "h_Eff_pt_TightID_ISO03dBeta_Mu2",         200, 0, 200);
  h_Eff_pt_TightID_ISO04dBeta[0]        = CreateH1F("h_Eff_pt_TightID_ISO04dBeta_Mu1", 
						    "h_Eff_pt_TightID_ISO04dBeta_Mu1",         200, 0, 200);
  h_Eff_pt_TightID_ISO04dBeta[1]        = CreateH1F("h_Eff_pt_TightID_ISO04dBeta_Mu2", 
						    "h_Eff_pt_TightID_ISO04dBeta_Mu2",         200, 0, 200);
  h_Eff_pt_TightID_ISO03PFWeighted[0]   = CreateH1F("h_Eff_pt_TightID_ISO03PFWeighted_Mu1", 
						    "h_Eff_pt_TightID_ISO03PFWeighted_Mu1",    200, 0, 200);
  h_Eff_pt_TightID_ISO03PFWeighted[1]   = CreateH1F("h_Eff_pt_TightID_ISO03PFWeighted_Mu2", 
						    "h_Eff_pt_TightID_ISO03PFWeighted_Mu2",    200, 0, 200);
  h_Eff_pt_TightID_ISO04PFWeighted[0]   = CreateH1F("h_Eff_pt_TightID_ISO04PFWeighted_Mu1", 
						    "h_Eff_pt_TightID_ISO04PFWeighted_Mu1",    200, 0, 200);
  h_Eff_pt_TightID_ISO04PFWeighted[1]   = CreateH1F("h_Eff_pt_TightID_ISO04PFWeighted_Mu2", 
						    "h_Eff_pt_TightID_ISO04PFWeighted_Mu2",    200, 0, 200);
  h_Eff_pt_TightID_ISO03PUPPI[0]        = CreateH1F("h_Eff_pt_TightID_ISO03PUPPI_Mu1", 
						    "h_Eff_pt_TightID_ISO03PUPPI_Mu1",         200, 0, 200);
  h_Eff_pt_TightID_ISO03PUPPI[1]        = CreateH1F("h_Eff_pt_TightID_ISO03PUPPI_Mu2", 
						    "h_Eff_pt_TightID_ISO03PUPPI_Mu2",         200, 0, 200);
  h_Eff_pt_TightID_ISO04PUPPI[0]        = CreateH1F("h_Eff_pt_TightID_ISO04PUPPI_Mu1", 
						    "h_Eff_pt_TightID_ISO04PUPPI_Mu1",         200, 0, 200);
  h_Eff_pt_TightID_ISO04PUPPI[1]        = CreateH1F("h_Eff_pt_TightID_ISO04PUPPI_Mu2", 
						    "h_Eff_pt_TightID_ISO04PUPPI_Mu2",         200, 0, 200);
  h_Eff_pt_MediumID_ISO03[0]            = CreateH1F("h_Eff_pt_MediumID_ISO03_Mu1", 
						    "h_Eff_pt_MediumID_ISO03_Mu1",             200, 0, 200);
  h_Eff_pt_MediumID_ISO03[1]            = CreateH1F("h_Eff_pt_MediumID_ISO03_Mu2", 
						    "h_Eff_pt_MediumID_ISO03_Mu2",             200, 0, 200);
  h_Eff_pt_MediumID_ISO04[0]            = CreateH1F("h_Eff_pt_MediumID_ISO04_Mu1", 
						    "h_Eff_pt_MediumID_ISO04_Mu1",             200, 0, 200);
  h_Eff_pt_MediumID_ISO04[1]            = CreateH1F("h_Eff_pt_MediumID_ISO04_Mu2", 
						    "h_Eff_pt_MediumID_ISO04_Mu2",             200, 0, 200);
  h_Eff_pt_MediumID_ISO03dBeta[0]       = CreateH1F("h_Eff_pt_MediumID_ISO03dBeta_Mu1", 
						    "h_Eff_pt_MediumID_ISO03dBeta_Mu1",        200, 0, 200);
  h_Eff_pt_MediumID_ISO03dBeta[1]       = CreateH1F("h_Eff_pt_MediumID_ISO03dBeta_Mu2", 
						    "h_Eff_pt_MediumID_ISO03dBeta_Mu2",        200, 0, 200);
  h_Eff_pt_MediumID_ISO04dBeta[0]       = CreateH1F("h_Eff_pt_MediumID_ISO04dBeta_Mu1", 
						    "h_Eff_pt_MediumID_ISO04dBeta_Mu1",        200, 0, 200);
  h_Eff_pt_MediumID_ISO04dBeta[1]       = CreateH1F("h_Eff_pt_MediumID_ISO04dBeta_Mu2", 
						    "h_Eff_pt_MediumID_ISO04dBeta_Mu2",        200, 0, 200);
  h_Eff_pt_MediumID_ISO03PFWeighted[0]  = CreateH1F("h_Eff_pt_MediumID_ISO03PFWeighted_Mu1", 
						    "h_Eff_pt_MediumID_ISO03PFWeighted_Mu1",   200, 0, 200);
  h_Eff_pt_MediumID_ISO03PFWeighted[1]  = CreateH1F("h_Eff_pt_MediumID_ISO03PFWeighted_Mu2", 
						    "h_Eff_pt_MediumID_ISO03PFWeighted_Mu2",   200, 0, 200);
  h_Eff_pt_MediumID_ISO04PFWeighted[0]  = CreateH1F("h_Eff_pt_MediumID_ISO04PFWeighted_Mu1", 
						    "h_Eff_pt_MediumID_ISO04PFWeighted_Mu1",   200, 0, 200);
  h_Eff_pt_MediumID_ISO04PFWeighted[1]  = CreateH1F("h_Eff_pt_MediumID_ISO04PFWeighted_Mu2", 
						    "h_Eff_pt_MediumID_ISO04PFWeighted_Mu2",   200, 0, 200);
  h_Eff_pt_MediumID_ISO03PUPPI[0]       = CreateH1F("h_Eff_pt_MediumID_ISO03PUPPI_Mu1", 
						    "h_Eff_pt_MediumID_ISO03PUPPI_Mu1",        200, 0, 200);
  h_Eff_pt_MediumID_ISO03PUPPI[1]       = CreateH1F("h_Eff_pt_MediumID_ISO03PUPPI_Mu2", 
						    "h_Eff_pt_MediumID_ISO03PUPPI_Mu2",        200, 0, 200);
  h_Eff_pt_MediumID_ISO04PUPPI[0]       = CreateH1F("h_Eff_pt_MediumID_ISO04PUPPI_Mu1", 
						    "h_Eff_pt_MediumID_ISO04PUPPI_Mu1",        200, 0, 200);
  h_Eff_pt_MediumID_ISO04PUPPI[1]       = CreateH1F("h_Eff_pt_MediumID_ISO04PUPPI_Mu2", 
						    "h_Eff_pt_MediumID_ISO04PUPPI_Mu2",        200, 0, 200);
  
  h_Eff_eta_Matched[0]                  = CreateH1F("h_Eff_eta_Matched_Mu1",
						    "h_Eff_eta_Matched_Mu1",                   50, -2.5, 2.5);
  h_Eff_eta_Matched[1]                  = CreateH1F("h_Eff_eta_Matched_Mu2",
						    "h_Eff_eta_Matched_Mu2",                   50, -2.5, 2.5);
  h_Eff_eta_TightID[0]                  = CreateH1F("h_Eff_eta_TightID_Mu1",
						    "h_Eff_eta_TightID_Mu1",                   50, -2.5, 2.5);
  h_Eff_eta_TightID[1]                  = CreateH1F("h_Eff_eta_TightID_Mu2",
						    "h_Eff_eta_TightID_Mu2",                   50, -2.5, 2.5);
  h_Eff_eta_MediumID[0]                 = CreateH1F("h_Eff_eta_MediumID_Mu1",
						    "h_Eff_eta_MediumID_Mu1",                  50, -2.5, 2.5);
  h_Eff_eta_MediumID[1]                 = CreateH1F("h_Eff_eta_MediumID_Mu2",
						    "h_Eff_eta_MediumID_Mu2",                  50, -2.5, 2.5);
  h_Eff_eta_HWWID[0]                    = CreateH1F("h_Eff_eta_HWWID_Mu1",
						    "h_Eff_eta_HWWID_Mu1",                     50, -2.5, 2.5);
  h_Eff_eta_HWWID[1]                    = CreateH1F("h_Eff_eta_HWWID_Mu2",
						    "h_Eff_eta_HWWID_Mu2",                     50, -2.5, 2.5);
  h_Eff_eta_TightIDGoT[0]               = CreateH1F("h_Eff_eta_TightIDGoT_Mu1",
						    "h_Eff_eta_TightIDGoT_Mu1",                50, -2.5, 2.5);
  h_Eff_eta_TightIDGoT[1]               = CreateH1F("h_Eff_eta_TightIDGoT_Mu2",
						    "h_Eff_eta_TightIDGoT_Mu2",                50, -2.5, 2.5);
  h_Eff_eta_TightIDipsHWW[0]            = CreateH1F("h_Eff_eta_TightIDipsHWW_Mu1",
						    "h_Eff_eta_TightIDipsHWW_Mu1",             50, -2.5, 2.5);
  h_Eff_eta_TightIDipsHWW[1]            = CreateH1F("h_Eff_eta_TightIDipsHWW_Mu2",
						    "h_Eff_eta_TightIDipsHWW_Mu2",             50, -2.5, 2.5);
  h_Eff_eta_MediumIDipsHWW[0]           = CreateH1F("h_Eff_eta_MediumIDipsHWW_Mu1",
						    "h_Eff_eta_MediumIDipsHWW_Mu1",            50, -2.5, 2.5);
  h_Eff_eta_MediumIDipsHWW[1]           = CreateH1F("h_Eff_eta_MediumIDipsHWW_Mu2",
						    "h_Eff_eta_MediumIDipsHWW_Mu2",            50, -2.5, 2.5);
  h_Eff_eta_TightID_ISO03[0]            = CreateH1F("h_Eff_eta_TightID_ISO03_Mu1", 
						    "h_Eff_eta_TightID_ISO03_Mu1",             50, -2.5, 2.5);
  h_Eff_eta_TightID_ISO03[1]            = CreateH1F("h_Eff_eta_TightID_ISO03_Mu2", 
						    "h_Eff_eta_TightID_ISO03_Mu2",             50, -2.5, 2.5);
  h_Eff_eta_TightID_ISO04[0]            = CreateH1F("h_Eff_eta_TightID_ISO04_Mu1", 
						    "h_Eff_eta_TightID_ISO04_Mu1",             50, -2.5, 2.5);
  h_Eff_eta_TightID_ISO04[1]            = CreateH1F("h_Eff_eta_TightID_ISO04_Mu2", 
						    "h_Eff_eta_TightID_ISO04_Mu2",             50, -2.5, 2.5);
  h_Eff_eta_TightID_ISO03dBeta[0]       = CreateH1F("h_Eff_eta_TightID_ISO03dBeta_Mu1", 
						    "h_Eff_eta_TightID_ISO03dBeta_Mu1",        50, -2.5, 2.5);
  h_Eff_eta_TightID_ISO03dBeta[1]       = CreateH1F("h_Eff_eta_TightID_ISO03dBeta_Mu2", 
						    "h_Eff_eta_TightID_ISO03dBeta_Mu2",        50, -2.5, 2.5);
  h_Eff_eta_TightID_ISO04dBeta[0]       = CreateH1F("h_Eff_eta_TightID_ISO04dBeta_Mu1", 
						    "h_Eff_eta_TightID_ISO04dBeta_Mu1",        50, -2.5, 2.5);
  h_Eff_eta_TightID_ISO04dBeta[1]       = CreateH1F("h_Eff_eta_TightID_ISO04dBeta_Mu2", 
						    "h_Eff_eta_TightID_ISO04dBeta_Mu2",        50, -2.5, 2.5);
  h_Eff_eta_TightID_ISO03PFWeighted[0]  = CreateH1F("h_Eff_eta_TightID_ISO03PFWeighted_Mu1", 
						    "h_Eff_eta_TightID_ISO03PFWeighted_Mu1",   50, -2.5, 2.5);
  h_Eff_eta_TightID_ISO03PFWeighted[1]  = CreateH1F("h_Eff_eta_TightID_ISO03PFWeighted_Mu2", 
						    "h_Eff_eta_TightID_ISO03PFWeighted_Mu2",   50, -2.5, 2.5);
  h_Eff_eta_TightID_ISO04PFWeighted[0]  = CreateH1F("h_Eff_eta_TightID_ISO04PFWeighted_Mu1", 
						    "h_Eff_eta_TightID_ISO04PFWeighted_Mu1",   50, -2.5, 2.5);
  h_Eff_eta_TightID_ISO04PFWeighted[1]  = CreateH1F("h_Eff_eta_TightID_ISO04PFWeighted_Mu2", 
						    "h_Eff_eta_TightID_ISO04PFWeighted_Mu2",   50, -2.5, 2.5);
  h_Eff_eta_TightID_ISO03PUPPI[0]       = CreateH1F("h_Eff_eta_TightID_ISO03PUPPI_Mu1", 
						    "h_Eff_eta_TightID_ISO03PUPPI_Mu1",        50, -2.5, 2.5);
  h_Eff_eta_TightID_ISO03PUPPI[1]       = CreateH1F("h_Eff_eta_TightID_ISO03PUPPI_Mu2", 
						    "h_Eff_eta_TightID_ISO03PUPPI_Mu2",        50, -2.5, 2.5);
  h_Eff_eta_TightID_ISO04PUPPI[0]       = CreateH1F("h_Eff_eta_TightID_ISO04PUPPI_Mu1", 
						    "h_Eff_eta_TightID_ISO04PUPPI_Mu1",        50, -2.5, 2.5);
  h_Eff_eta_TightID_ISO04PUPPI[1]       = CreateH1F("h_Eff_eta_TightID_ISO04PUPPI_Mu2", 
						    "h_Eff_eta_TightID_ISO04PUPPI_Mu2",        50, -2.5, 2.5);
  h_Eff_eta_MediumID_ISO03[0]           = CreateH1F("h_Eff_eta_MediumID_ISO03_Mu1", 
						    "h_Eff_eta_MediumID_ISO03_Mu1",            50, -2.5, 2.5);
  h_Eff_eta_MediumID_ISO03[1]           = CreateH1F("h_Eff_eta_MediumID_ISO03_Mu2", 
						    "h_Eff_eta_MediumID_ISO03_Mu2",            50, -2.5, 2.5);
  h_Eff_eta_MediumID_ISO04[0]           = CreateH1F("h_Eff_eta_MediumID_ISO04_Mu1", 
						    "h_Eff_eta_MediumID_ISO04_Mu1",            50, -2.5, 2.5);
  h_Eff_eta_MediumID_ISO04[1]           = CreateH1F("h_Eff_eta_MediumID_ISO04_Mu2", 
						    "h_Eff_eta_MediumID_ISO04_Mu2",            50, -2.5, 2.5);
  h_Eff_eta_MediumID_ISO03dBeta[0]      = CreateH1F("h_Eff_eta_MediumID_ISO03dBeta_Mu1", 
						    "h_Eff_eta_MediumID_ISO03dBeta_Mu1",       50, -2.5, 2.5);
  h_Eff_eta_MediumID_ISO03dBeta[1]      = CreateH1F("h_Eff_eta_MediumID_ISO03dBeta_Mu2", 
						    "h_Eff_eta_MediumID_ISO03dBeta_Mu2",       50, -2.5, 2.5);
  h_Eff_eta_MediumID_ISO04dBeta[0]      = CreateH1F("h_Eff_eta_MediumID_ISO04dBeta_Mu1", 
						    "h_Eff_eta_MediumID_ISO04dBeta_Mu1",       50, -2.5, 2.5);
  h_Eff_eta_MediumID_ISO04dBeta[1]      = CreateH1F("h_Eff_eta_MediumID_ISO04dBeta_Mu2", 
						    "h_Eff_eta_MediumID_ISO04dBeta_Mu2",       50, -2.5, 2.5);
  h_Eff_eta_MediumID_ISO03PFWeighted[0] = CreateH1F("h_Eff_eta_MediumID_ISO03PFWeighted_Mu1", 
						    "h_Eff_eta_MediumID_ISO03PFWeighted_Mu1",  50, -2.5, 2.5);
  h_Eff_eta_MediumID_ISO03PFWeighted[1] = CreateH1F("h_Eff_eta_MediumID_ISO03PFWeighted_Mu2", 
						    "h_Eff_eta_MediumID_ISO03PFWeighted_Mu2",  50, -2.5, 2.5);
  h_Eff_eta_MediumID_ISO04PFWeighted[0] = CreateH1F("h_Eff_eta_MediumID_ISO04PFWeighted_Mu1", 
						    "h_Eff_eta_MediumID_ISO04PFWeighted_Mu1",  50, -2.5, 2.5);
  h_Eff_eta_MediumID_ISO04PFWeighted[1] = CreateH1F("h_Eff_eta_MediumID_ISO04PFWeighted_Mu2", 
						    "h_Eff_eta_MediumID_ISO04PFWeighted_Mu2",  50, -2.5, 2.5);
  h_Eff_eta_MediumID_ISO03PUPPI[0]      = CreateH1F("h_Eff_eta_MediumID_ISO03PUPPI_Mu1", 
						    "h_Eff_eta_MediumID_ISO03PUPPI_Mu1",       50, -2.5, 2.5);
  h_Eff_eta_MediumID_ISO03PUPPI[1]      = CreateH1F("h_Eff_eta_MediumID_ISO03PUPPI_Mu2", 
						    "h_Eff_eta_MediumID_ISO03PUPPI_Mu2",       50, -2.5, 2.5);
  h_Eff_eta_MediumID_ISO04PUPPI[0]      = CreateH1F("h_Eff_eta_MediumID_ISO04PUPPI_Mu1", 
						    "h_Eff_eta_MediumID_ISO04PUPPI_Mu1",       50, -2.5, 2.5);
  h_Eff_eta_MediumID_ISO04PUPPI[1]      = CreateH1F("h_Eff_eta_MediumID_ISO04PUPPI_Mu2", 
						    "h_Eff_eta_MediumID_ISO04PUPPI_Mu2",       50, -2.5, 2.5);
  
  h_Eff_npv_Matched[0]                  = CreateH1F("h_Eff_npv_Matched_Mu1",
						    "h_Eff_npv_Matched_Mu1",                   45, 0, 45);
  h_Eff_npv_Matched[1]                  = CreateH1F("h_Eff_npv_Matched_Mu2",
						    "h_Eff_npv_Matched_Mu2",                   45, 0, 45);
  h_Eff_npv_TightID[0]                  = CreateH1F("h_Eff_npv_TightID_Mu1",
						    "h_Eff_npv_TightID_Mu1",                   45, 0, 45);
  h_Eff_npv_TightID[1]                  = CreateH1F("h_Eff_npv_TightID_Mu2",
						    "h_Eff_npv_TightID_Mu2",                   45, 0, 45);
  h_Eff_npv_MediumID[0]                 = CreateH1F("h_Eff_npv_MediumID_Mu1",
						    "h_Eff_npv_MediumID_Mu1",                  45, 0, 45);
  h_Eff_npv_MediumID[1]                 = CreateH1F("h_Eff_npv_MediumID_Mu2",
						    "h_Eff_npv_MediumID_Mu2",                  45, 0, 45);
  h_Eff_npv_HWWID[0]                    = CreateH1F("h_Eff_npv_HWWID_Mu1",
						    "h_Eff_npv_HWWID_Mu1",                     45, 0, 45);
  h_Eff_npv_HWWID[1]                    = CreateH1F("h_Eff_npv_HWWID_Mu2",
						    "h_Eff_npv_HWWID_Mu2",                     45, 0, 45);
  h_Eff_npv_TightIDGoT[0]               = CreateH1F("h_Eff_npv_TightIDGoT_Mu1",
						    "h_Eff_npv_TightIDGoT_Mu1",                45, 0, 45);
  h_Eff_npv_TightIDGoT[1]               = CreateH1F("h_Eff_npv_TightIDGoT_Mu2",
						    "h_Eff_npv_TightIDGoT_Mu2",                45, 0, 45);
  h_Eff_npv_TightIDipsHWW[0]            = CreateH1F("h_Eff_npv_TightIDipsHWW_Mu1",
						    "h_Eff_npv_TightIDipsHWW_Mu1",             45, 0, 45);
  h_Eff_npv_TightIDipsHWW[1]            = CreateH1F("h_Eff_npv_TightIDipsHWW_Mu2",
						    "h_Eff_npv_TightIDipsHWW_Mu2",             45, 0, 45);
  h_Eff_npv_MediumIDipsHWW[0]           = CreateH1F("h_Eff_npv_MediumIDipsHWW_Mu1",
						    "h_Eff_npv_MediumIDipsHWW_Mu1",            45, 0, 45);
  h_Eff_npv_MediumIDipsHWW[1]           = CreateH1F("h_Eff_npv_MediumIDipsHWW_Mu2",
						    "h_Eff_npv_MediumIDipsHWW_Mu2",            45, 0, 45);
  h_Eff_npv_TightID_ISO03[0]            = CreateH1F("h_Eff_npv_TightID_ISO03_Mu1", 
						    "h_Eff_npv_TightID_ISO03_Mu1",             45, 0, 45);
  h_Eff_npv_TightID_ISO03[1]            = CreateH1F("h_Eff_npv_TightID_ISO03_Mu2", 
						    "h_Eff_npv_TightID_ISO03_Mu2",             45, 0, 45);
  h_Eff_npv_TightID_ISO04[0]            = CreateH1F("h_Eff_npv_TightID_ISO04_Mu1", 
						    "h_Eff_npv_TightID_ISO04_Mu1",             45, 0, 45);
  h_Eff_npv_TightID_ISO04[1]            = CreateH1F("h_Eff_npv_TightID_ISO04_Mu2", 
						    "h_Eff_npv_TightID_ISO04_Mu2",             45, 0, 45);
  h_Eff_npv_TightID_ISO03dBeta[0]       = CreateH1F("h_Eff_npv_TightID_ISO03dBeta_Mu1", 
						    "h_Eff_npv_TightID_ISO03dBeta_Mu1",        45, 0, 45);
  h_Eff_npv_TightID_ISO03dBeta[1]       = CreateH1F("h_Eff_npv_TightID_ISO03dBeta_Mu2", 
						    "h_Eff_npv_TightID_ISO03dBeta_Mu2",        45, 0, 45);
  h_Eff_npv_TightID_ISO04dBeta[0]       = CreateH1F("h_Eff_npv_TightID_ISO04dBeta_Mu1", 
						    "h_Eff_npv_TightID_ISO04dBeta_Mu1",        45, 0, 45);
  h_Eff_npv_TightID_ISO04dBeta[1]       = CreateH1F("h_Eff_npv_TightID_ISO04dBeta_Mu2", 
						    "h_Eff_npv_TightID_ISO04dBeta_Mu2",        45, 0, 45);
  h_Eff_npv_TightID_ISO03PFWeighted[0]  = CreateH1F("h_Eff_npv_TightID_ISO03PFWeighted_Mu1", 
						    "h_Eff_npv_TightID_ISO03PFWeighted_Mu1",   45, 0, 45);
  h_Eff_npv_TightID_ISO03PFWeighted[1]  = CreateH1F("h_Eff_npv_TightID_ISO03PFWeighted_Mu2", 
						    "h_Eff_npv_TightID_ISO03PFWeighted_Mu2",   45, 0, 45);
  h_Eff_npv_TightID_ISO04PFWeighted[0]  = CreateH1F("h_Eff_npv_TightID_ISO04PFWeighted_Mu1", 
						    "h_Eff_npv_TightID_ISO04PFWeighted_Mu1",   45, 0, 45);
  h_Eff_npv_TightID_ISO04PFWeighted[1]  = CreateH1F("h_Eff_npv_TightID_ISO04PFWeighted_Mu2", 
						    "h_Eff_npv_TightID_ISO04PFWeighted_Mu2",   45, 0, 45);
  h_Eff_npv_TightID_ISO03PUPPI[0]       = CreateH1F("h_Eff_npv_TightID_ISO03PUPPI_Mu1", 
						    "h_Eff_npv_TightID_ISO03PUPPI_Mu1",        45, 0, 45);
  h_Eff_npv_TightID_ISO03PUPPI[1]       = CreateH1F("h_Eff_npv_TightID_ISO03PUPPI_Mu2", 
						    "h_Eff_npv_TightID_ISO03PUPPI_Mu2",        45, 0, 45);
  h_Eff_npv_TightID_ISO04PUPPI[0]       = CreateH1F("h_Eff_npv_TightID_ISO04PUPPI_Mu1", 
						    "h_Eff_npv_TightID_ISO04PUPPI_Mu1",        45, 0, 45);
  h_Eff_npv_TightID_ISO04PUPPI[1]       = CreateH1F("h_Eff_npv_TightID_ISO04PUPPI_Mu2", 
						    "h_Eff_npv_TightID_ISO04PUPPI_Mu2",        45, 0, 45);
  h_Eff_npv_MediumID_ISO03[0]           = CreateH1F("h_Eff_npv_MediumID_ISO03_Mu1", 
						    "h_Eff_npv_MediumID_ISO03_Mu1",            45, 0, 45);
  h_Eff_npv_MediumID_ISO03[1]           = CreateH1F("h_Eff_npv_MediumID_ISO03_Mu2", 
						    "h_Eff_npv_MediumID_ISO03_Mu2",            45, 0, 45);
  h_Eff_npv_MediumID_ISO04[0]           = CreateH1F("h_Eff_npv_MediumID_ISO04_Mu1", 
						    "h_Eff_npv_MediumID_ISO04_Mu1",            45, 0, 45);
  h_Eff_npv_MediumID_ISO04[1]           = CreateH1F("h_Eff_npv_MediumID_ISO04_Mu2", 
						    "h_Eff_npv_MediumID_ISO04_Mu2",            45, 0, 45);
  h_Eff_npv_MediumID_ISO03dBeta[0]      = CreateH1F("h_Eff_npv_MediumID_ISO03dBeta_Mu1", 
						    "h_Eff_npv_MediumID_ISO03dBeta_Mu1",       45, 0, 45);
  h_Eff_npv_MediumID_ISO03dBeta[1]      = CreateH1F("h_Eff_npv_MediumID_ISO03dBeta_Mu2", 
						    "h_Eff_npv_MediumID_ISO03dBeta_Mu2",       45, 0, 45);
  h_Eff_npv_MediumID_ISO04dBeta[0]      = CreateH1F("h_Eff_npv_MediumID_ISO04dBeta_Mu1", 
						    "h_Eff_npv_MediumID_ISO04dBeta_Mu1",       45, 0, 45);
  h_Eff_npv_MediumID_ISO04dBeta[1]      = CreateH1F("h_Eff_npv_MediumID_ISO04dBeta_Mu2", 
						    "h_Eff_npv_MediumID_ISO04dBeta_Mu2",       45, 0, 45);
  h_Eff_npv_MediumID_ISO03PFWeighted[0] = CreateH1F("h_Eff_npv_MediumID_ISO03PFWeighted_Mu1", 
						    "h_Eff_npv_MediumID_ISO03PFWeighted_Mu1",  45, 0, 45);
  h_Eff_npv_MediumID_ISO03PFWeighted[1] = CreateH1F("h_Eff_npv_MediumID_ISO03PFWeighted_Mu2", 
						    "h_Eff_npv_MediumID_ISO03PFWeighted_Mu2",  45, 0, 45);
  h_Eff_npv_MediumID_ISO04PFWeighted[0] = CreateH1F("h_Eff_npv_MediumID_ISO04PFWeighted_Mu1", 
						    "h_Eff_npv_MediumID_ISO04PFWeighted_Mu1",  45, 0, 45);
  h_Eff_npv_MediumID_ISO04PFWeighted[1] = CreateH1F("h_Eff_npv_MediumID_ISO04PFWeighted_Mu2", 
						    "h_Eff_npv_MediumID_ISO04PFWeighted_Mu2",  45, 0, 45);
  h_Eff_npv_MediumID_ISO03PUPPI[0]      = CreateH1F("h_Eff_npv_MediumID_ISO03PUPPI_Mu1", 
						    "h_Eff_npv_MediumID_ISO03PUPPI_Mu1",       45, 0, 45);
  h_Eff_npv_MediumID_ISO03PUPPI[1]      = CreateH1F("h_Eff_npv_MediumID_ISO03PUPPI_Mu2", 
						    "h_Eff_npv_MediumID_ISO03PUPPI_Mu2",       45, 0, 45);
  h_Eff_npv_MediumID_ISO04PUPPI[0]      = CreateH1F("h_Eff_npv_MediumID_ISO04PUPPI_Mu1", 
						    "h_Eff_npv_MediumID_ISO04PUPPI_Mu1",       45, 0, 45);
  h_Eff_npv_MediumID_ISO04PUPPI[1]      = CreateH1F("h_Eff_npv_MediumID_ISO04PUPPI_Mu2", 
						    "h_Eff_npv_MediumID_ISO04PUPPI_Mu2",       45, 0, 45);

  // Dimuon ISO efficiencies vs pt, eta and npv

  h_Eff_pt_TightID_Dilep                   = CreateH1F("h_Eff_pt_TightID_Dilep", 
						       "h_Eff_pt_TightID_Dilep",                    200, 0, 200);
  h_Eff_pt_TightID_ISO03_Dilep             = CreateH1F("h_Eff_pt_TightID_ISO03_Dilep", 
						       "h_Eff_pt_TightID_ISO03_Dilep",              200, 0, 200);
  h_Eff_pt_TightID_ISO04_Dilep             = CreateH1F("h_Eff_pt_TightID_ISO04_Dilep", 
						       "h_Eff_pt_TightID_ISO04_Dilep",              200, 0, 200);
  h_Eff_pt_TightID_ISO03dBeta_Dilep        = CreateH1F("h_Eff_pt_TightID_ISO03dBeta_Dilep", 
						       "h_Eff_pt_TightID_ISO03dBeta_Dilep",         200, 0, 200);
  h_Eff_pt_TightID_ISO04dBeta_Dilep        = CreateH1F("h_Eff_pt_TightID_ISO04dBeta_Dilep", 
						       "h_Eff_pt_TightID_ISO04dBeta_Dilep",         200, 0, 200);
  h_Eff_pt_TightID_ISO03PFWeighted_Dilep   = CreateH1F("h_Eff_pt_TightID_ISO03PFWeighted_Dilep", 
						       "h_Eff_pt_TightID_ISO03PFWeighted_Dilep",    200, 0, 200);
  h_Eff_pt_TightID_ISO04PFWeighted_Dilep   = CreateH1F("h_Eff_pt_TightID_ISO04PFWeighted_Dilep", 
						       "h_Eff_pt_TightID_ISO04PFWeighted_Dilep",    200, 0, 200);
  h_Eff_pt_TightID_ISO03PUPPI_Dilep        = CreateH1F("h_Eff_pt_TightID_ISO03PUPPI_Dilep", 
						       "h_Eff_pt_TightID_ISO03PUPPI_Dilep",         200, 0, 200);
  h_Eff_pt_TightID_ISO04PUPPI_Dilep        = CreateH1F("h_Eff_pt_TightID_ISO04PUPPI_Dilep", 
						       "h_Eff_pt_TightID_ISO04PUPPI_Dilep",         200, 0, 200);
  h_Eff_pt_MediumID_Dilep                  = CreateH1F("h_Eff_pt_MediumID_Dilep", 
						       "h_Eff_pt_MediumID_Dilep",                   200, 0, 200);
  h_Eff_pt_MediumID_ISO03_Dilep            = CreateH1F("h_Eff_pt_MediumID_ISO03_Dilep", 
						       "h_Eff_pt_MediumID_ISO03_Dilep",             200, 0, 200);
  h_Eff_pt_MediumID_ISO04_Dilep            = CreateH1F("h_Eff_pt_MediumID_ISO04_Dilep", 
						       "h_Eff_pt_MediumID_ISO04_Dilep",             200, 0, 200);
  h_Eff_pt_MediumID_ISO03dBeta_Dilep       = CreateH1F("h_Eff_pt_MediumID_ISO03dBeta_Dilep", 
						       "h_Eff_pt_MediumID_ISO03dBeta_Dilep",        200, 0, 200);
  h_Eff_pt_MediumID_ISO04dBeta_Dilep       = CreateH1F("h_Eff_pt_MediumID_ISO04dBeta_Dilep", 
						       "h_Eff_pt_MediumID_ISO04dBeta_Dilep",        200, 0, 200);
  h_Eff_pt_MediumID_ISO03PFWeighted_Dilep  = CreateH1F("h_Eff_pt_MediumID_ISO03PFWeighted_Dilep", 
						       "h_Eff_pt_MediumID_ISO03PFWeighted_Dilep",   200, 0, 200);
  h_Eff_pt_MediumID_ISO04PFWeighted_Dilep  = CreateH1F("h_Eff_pt_MediumID_ISO04PFWeighted_Dilep", 
						       "h_Eff_pt_MediumID_ISO04PFWeighted_Dilep",   200, 0, 200);
  h_Eff_pt_MediumID_ISO03PUPPI_Dilep       = CreateH1F("h_Eff_pt_MediumID_ISO03PUPPI_Dilep", 
						       "h_Eff_pt_MediumID_ISO03PUPPI_Dilep",        200, 0, 200);
  h_Eff_pt_MediumID_ISO04PUPPI_Dilep       = CreateH1F("h_Eff_pt_MediumID_ISO04PUPPI_Dilep", 
						       "h_Eff_pt_MediumID_ISO04PUPPI_Dilep",        200, 0, 200);


  h_Eff_eta_TightID_Dilep                  = CreateH1F("h_Eff_eta_TightID_Dilep", 
						       "h_Eff_eta_TightID_Dilep",                   50, -2.5, 2.5);
  h_Eff_eta_TightID_ISO03_Dilep            = CreateH1F("h_Eff_eta_TightID_ISO03_Dilep", 
						       "h_Eff_eta_TightID_ISO03_Dilep",             50, -2.5, 2.5);
  h_Eff_eta_TightID_ISO04_Dilep            = CreateH1F("h_Eff_eta_TightID_ISO04_Dilep", 
						       "h_Eff_eta_TightID_ISO04_Dilep",             50, -2.5, 2.5);
  h_Eff_eta_TightID_ISO03dBeta_Dilep       = CreateH1F("h_Eff_eta_TightID_ISO03dBeta_Dilep", 
						       "h_Eff_eta_TightID_ISO03dBeta_Dilep",        50, -2.5, 2.5);
  h_Eff_eta_TightID_ISO04dBeta_Dilep       = CreateH1F("h_Eff_eta_TightID_ISO04dBeta_Dilep", 
						       "h_Eff_eta_TightID_ISO04dBeta_Dilep",        50, -2.5, 2.5);
  h_Eff_eta_TightID_ISO03PFWeighted_Dilep  = CreateH1F("h_Eff_eta_TightID_ISO03PFWeighted_Dilep", 
						       "h_Eff_eta_TightID_ISO03PFWeighted_Dilep",   50, -2.5, 2.5);
  h_Eff_eta_TightID_ISO04PFWeighted_Dilep  = CreateH1F("h_Eff_eta_TightID_ISO04PFWeighted_Dilep", 
						       "h_Eff_eta_TightID_ISO04PFWeighted_Dilep",   50, -2.5, 2.5);
  h_Eff_eta_TightID_ISO03PUPPI_Dilep       = CreateH1F("h_Eff_eta_TightID_ISO03PUPPI_Dilep", 
						       "h_Eff_eta_TightID_ISO03PUPPI_Dilep",        50, -2.5, 2.5);
  h_Eff_eta_TightID_ISO04PUPPI_Dilep       = CreateH1F("h_Eff_eta_TightID_ISO04PUPPI_Dilep", 
						       "h_Eff_eta_TightID_ISO04PUPPI_Dilep",        50, -2.5, 2.5);
  h_Eff_eta_MediumID_Dilep                 = CreateH1F("h_Eff_eta_MediumID_Dilep", 
						       "h_Eff_eta_MediumID_Dilep",                  50, -2.5, 2.5);
  h_Eff_eta_MediumID_ISO03_Dilep           = CreateH1F("h_Eff_eta_MediumID_ISO03_Dilep", 
						       "h_Eff_eta_MediumID_ISO03_Dilep",            50, -2.5, 2.5);
  h_Eff_eta_MediumID_ISO04_Dilep           = CreateH1F("h_Eff_eta_MediumID_ISO04_Dilep", 
						       "h_Eff_eta_MediumID_ISO04_Dilep",            50, -2.5, 2.5);
  h_Eff_eta_MediumID_ISO03dBeta_Dilep      = CreateH1F("h_Eff_eta_MediumID_ISO03dBeta_Dilep", 
						       "h_Eff_eta_MediumID_ISO03dBeta_Dilep",       50, -2.5, 2.5);
  h_Eff_eta_MediumID_ISO04dBeta_Dilep      = CreateH1F("h_Eff_eta_MediumID_ISO04dBeta_Dilep", 
						       "h_Eff_eta_MediumID_ISO04dBeta_Dilep",       50, -2.5, 2.5);
  h_Eff_eta_MediumID_ISO03PFWeighted_Dilep = CreateH1F("h_Eff_eta_MediumID_ISO03PFWeighted_Dilep", 
						       "h_Eff_eta_MediumID_ISO03PFWeighted_Dilep",  50, -2.5, 2.5);
  h_Eff_eta_MediumID_ISO04PFWeighted_Dilep = CreateH1F("h_Eff_eta_MediumID_ISO04PFWeighted_Dilep", 
						       "h_Eff_eta_MediumID_ISO04PFWeighted_Dilep",  50, -2.5, 2.5);
  h_Eff_eta_MediumID_ISO03PUPPI_Dilep      = CreateH1F("h_Eff_eta_MediumID_ISO03PUPPI_Dilep", 
						       "h_Eff_eta_MediumID_ISO03PUPPI_Dilep",       50, -2.5, 2.5);
  h_Eff_eta_MediumID_ISO04PUPPI_Dilep      = CreateH1F("h_Eff_eta_MediumID_ISO04PUPPI_Dilep", 
						       "h_Eff_eta_MediumID_ISO04PUPPI_Dilep",       50, -2.5, 2.5);
  

  h_Eff_npv_TightID_Dilep                  = CreateH1F("h_Eff_npv_TightID_Dilep", 
						       "h_Eff_npv_TightID_Dilep",                   45, 0, 45);
  h_Eff_npv_TightID_ISO03_Dilep            = CreateH1F("h_Eff_npv_TightID_ISO03_Dilep", 
						       "h_Eff_npv_TightID_ISO03_Dilep",             45, 0, 45);
  h_Eff_npv_TightID_ISO04_Dilep            = CreateH1F("h_Eff_npv_TightID_ISO04_Dilep", 
						       "h_Eff_npv_TightID_ISO04_Dilep",             45, 0, 45);
  h_Eff_npv_TightID_ISO03dBeta_Dilep       = CreateH1F("h_Eff_npv_TightID_ISO03dBeta_Dilep", 
						       "h_Eff_npv_TightID_ISO03dBeta_Dilep",        45, 0, 45);
  h_Eff_npv_TightID_ISO04dBeta_Dilep       = CreateH1F("h_Eff_npv_TightID_ISO04dBeta_Dilep", 
						       "h_Eff_npv_TightID_ISO04dBeta_Dilep",        45, 0, 45);
  h_Eff_npv_TightID_ISO03PFWeighted_Dilep  = CreateH1F("h_Eff_npv_TightID_ISO03PFWeighted_Dilep", 
						       "h_Eff_npv_TightID_ISO03PFWeighted_Dilep",   45, 0, 45);
  h_Eff_npv_TightID_ISO04PFWeighted_Dilep  = CreateH1F("h_Eff_npv_TightID_ISO04PFWeighted_Dilep", 
						       "h_Eff_npv_TightID_ISO04PFWeighted_Dilep",   45, 0, 45);
  h_Eff_npv_TightID_ISO03PUPPI_Dilep       = CreateH1F("h_Eff_npv_TightID_ISO03PUPPI_Dilep", 
						       "h_Eff_npv_TightID_ISO03PUPPI_Dilep",        45, 0, 45);
  h_Eff_npv_TightID_ISO04PUPPI_Dilep       = CreateH1F("h_Eff_npv_TightID_ISO04PUPPI_Dilep", 
						       "h_Eff_npv_TightID_ISO04PUPPI_Dilep",        45, 0, 45);
  h_Eff_npv_MediumID_Dilep                 = CreateH1F("h_Eff_npv_MediumID_Dilep", 
						       "h_Eff_npv_MediumID_Dilep",                  45, 0, 45);
  h_Eff_npv_MediumID_ISO03_Dilep           = CreateH1F("h_Eff_npv_MediumID_ISO03_Dilep", 
						       "h_Eff_npv_MediumID_ISO03_Dilep",            45, 0, 45);
  h_Eff_npv_MediumID_ISO04_Dilep           = CreateH1F("h_Eff_npv_MediumID_ISO04_Dilep", 
						       "h_Eff_npv_MediumID_ISO04_Dilep",            45, 0, 45);
  h_Eff_npv_MediumID_ISO03dBeta_Dilep      = CreateH1F("h_Eff_npv_MediumID_ISO03dBeta_Dilep", 
						       "h_Eff_npv_MediumID_ISO03dBeta_Dilep",       45, 0, 45);
  h_Eff_npv_MediumID_ISO04dBeta_Dilep      = CreateH1F("h_Eff_npv_MediumID_ISO04dBeta_Dilep", 
						       "h_Eff_npv_MediumID_ISO04dBeta_Dilep",       45, 0, 45);
  h_Eff_npv_MediumID_ISO03PFWeighted_Dilep = CreateH1F("h_Eff_npv_MediumID_ISO03PFWeighted_Dilep", 
						       "h_Eff_npv_MediumID_ISO03PFWeighted_Dilep",  45, 0, 45);
  h_Eff_npv_MediumID_ISO04PFWeighted_Dilep = CreateH1F("h_Eff_npv_MediumID_ISO04PFWeighted_Dilep", 
						       "h_Eff_npv_MediumID_ISO04PFWeighted_Dilep",  45, 0, 45);
  h_Eff_npv_MediumID_ISO03PUPPI_Dilep      = CreateH1F("h_Eff_npv_MediumID_ISO03PUPPI_Dilep", 
						       "h_Eff_npv_MediumID_ISO03PUPPI_Dilep",       45, 0, 45);
  h_Eff_npv_MediumID_ISO04PUPPI_Dilep      = CreateH1F("h_Eff_npv_MediumID_ISO04PUPPI_Dilep", 
						       "h_Eff_npv_MediumID_ISO04PUPPI_Dilep",       45, 0, 45);

  //ISO ROC Curves
  h_RC_TightID_ISO03_Dilep            = CreateH1F("h_RC_TightID_ISO03_Dilep", 
						  "h_RC_TightID_ISO03_Dilep",             50, 0, 50);
  h_RC_TightID_ISO04_Dilep            = CreateH1F("h_RC_TightID_ISO04_Dilep", 
						  "h_RC_TightID_ISO04_Dilep",             50, 0, 50);
  h_RC_TightID_ISO03dBeta_Dilep       = CreateH1F("h_RC_TightID_ISO03dBeta_Dilep", 
						  "h_RC_TightID_ISO03dBeta_Dilep",        50, 0, 50);
  h_RC_TightID_ISO04dBeta_Dilep       = CreateH1F("h_RC_TightID_ISO04dBeta_Dilep", 
						  "h_RC_TightID_ISO04dBeta_Dilep",        50, 0, 50);
  h_RC_TightID_ISO03PFWeighted_Dilep  = CreateH1F("h_RC_TightID_ISO03PFWeighted_Dilep", 
						  "h_RC_TightID_ISO03PFWeighted_Dilep",   50, 0, 50);
  h_RC_TightID_ISO04PFWeighted_Dilep  = CreateH1F("h_RC_TightID_ISO04PFWeighted_Dilep", 
						  "h_RC_TightID_ISO04PFWeighted_Dilep",   50, 0, 50);
  h_RC_TightID_ISO03PUPPI_Dilep       = CreateH1F("h_RC_TightID_ISO03PUPPI_Dilep", 
						  "h_RC_TightID_ISO03PUPPI_Dilep",        50, 0, 50);
  h_RC_TightID_ISO04PUPPI_Dilep       = CreateH1F("h_RC_TightID_ISO04PUPPI_Dilep", 
						  "h_RC_TightID_ISO04PUPPI_Dilep",        50, 0, 50);
  h_RC_MediumID_ISO03_Dilep           = CreateH1F("h_RC_MediumID_ISO03_Dilep", 
						  "h_RC_MediumID_ISO03_Dilep",            50, 0, 50);
  h_RC_MediumID_ISO04_Dilep           = CreateH1F("h_RC_MediumID_ISO04_Dilep", 
						  "h_RC_MediumID_ISO04_Dilep",            50, 0, 50);
  h_RC_MediumID_ISO03dBeta_Dilep      = CreateH1F("h_RC_MediumID_ISO03dBeta_Dilep", 
						  "h_RC_MediumID_ISO03dBeta_Dilep",       50, 0, 50);
  h_RC_MediumID_ISO04dBeta_Dilep      = CreateH1F("h_RC_MediumID_ISO04dBeta_Dilep", 
						  "h_RC_MediumID_ISO04dBeta_Dilep",       50, 0, 50);
  h_RC_MediumID_ISO03PFWeighted_Dilep = CreateH1F("h_RC_MediumID_ISO03PFWeighted_Dilep", 
						  "h_RC_MediumID_ISO03PFWeighted_Dilep",  50, 0, 50);
  h_RC_MediumID_ISO04PFWeighted_Dilep = CreateH1F("h_RC_MediumID_ISO04PFWeighted_Dilep", 
						  "h_RC_MediumID_ISO04PFWeighted_Dilep",  50, 0, 50);
  h_RC_MediumID_ISO03PUPPI_Dilep      = CreateH1F("h_RC_MediumID_ISO03PUPPI_Dilep", 
						  "h_RC_MediumID_ISO03PUPPI_Dilep",       50, 0, 50);
  h_RC_MediumID_ISO04PUPPI_Dilep      = CreateH1F("h_RC_MediumID_ISO04PUPPI_Dilep", 
						  "h_RC_MediumID_ISO04PUPPI_Dilep",       50, 0, 50); 

}


void CoreMuonSelector::InsideLoop() {
 
 // The InsideLoop() function is called for each entry in the tree to be processed  

  if (_Debug) std::cout << "[DEBUG][Event "<< Get<int>("T_Event_EventNumber") <<"]" << std::endl;

  //------------------------------------------------------------------------------
  // Initialise data members for each event
  //------------------------------------------------------------------------------

  // GEN Info
  G_GEN_PromptMuon_4vec.clear();
  G_GEN_Muon_4vec.clear();
  G_GEN_isMuMu      = false;
  G_GEN_isMuTau     = false;
  G_GEN_isTauMu     = false;
  G_GEN_isTauTau    = false;
  G_GEN_isNonPrompt = false;

  // RECO muons
  G_Muon_4vec.clear();

  G_MuonID_Tight.clear();
  G_MuonID_Medium.clear();
  G_MuonID_HWW.clear();
  G_MuonID_Tight_GoT.clear();
  G_MuonID_IPs_HWW.clear();
  G_MuonID_GLBorTRKArb.clear();
  G_MuonID_Fiducial.clear();

  G_MuonISO03.clear();
  G_MuonISO03_dBeta.clear();
  G_MuonISO03_PFWeighted.clear();
  G_MuonISO03_PUPPI.clear();
  G_MuonISO04.clear();
  G_MuonISO04_dBeta.clear();
  G_MuonISO04_PFWeighted.clear();
  G_MuonISO04_PUPPI.clear();

  G_Muon_ChCompatible.clear();
  G_Muon_Matching.clear();
  
  // Sizes
  G_RecoMuSize  = 0;
  G_NPV         = 0;

  G_GEN_Pass     = false;
  G_PassMatching = false;

  // Event Flags
  EvtFlag_Fiducial = false;
  EvtFlag_Gen      = false;
  EvtFlag_Matching = false;


  //------------------------------------------------------------------------------
  // Get all RECO muons
  //------------------------------------------------------------------------------

  //G_RecoMuSize = Get<std::vector<float>*>("T_Muon_Px")->size();
  G_RecoMuSize = GetSizeOf("T_Muon_Px");
  
  for (unsigned int i = 0; i < G_RecoMuSize; ++i) {
    
    //-->define the 4D momentum for each RECO muon. 
    G_Muon_4vec.push_back(TLorentzVector(Get<float>("T_Muon_Px",i), Get<float>("T_Muon_Py",i), 
					 Get<float>("T_Muon_Pz",i), Get<float>("T_Muon_Energy",i)));
    
  }

  //------------------------------------------------------------------------------
  // Check if RECO muons pass IDs and ISOs 
  //------------------------------------------------------------------------------
  
  CheckMuons();

  //------------------------------------------------------------------------------
  // Get all GEN prompt muons
  //------------------------------------------------------------------------------
  
  if (!_IsDATA) SetGenInfo();

  //------------------------------------------------------------------------------
  // Get RECO-GEN matching
  //------------------------------------------------------------------------------
  
  GetMatching();

  //------------------------------------------------------------------------------
  // Get number of good vertex per event
  //------------------------------------------------------------------------------

  //G_NPV = Get<std::vector<float>*>("T_Vertex_z")->size();
  G_NPV = GetSizeOf("T_Vertex_z");
  h_N_PV->Fill(G_NPV, _factN);

  //------------------------------------------------------------------------------
  // Set Event Flags
  //------------------------------------------------------------------------------

  SetEventFlags();

  //------------------------------------------------------------------------------
  // Do Efficiencies in function of Pt, Eta and NPV for muons that passed matching
  //------------------------------------------------------------------------------

  if (EvtFlag_Matching) {
    doEffsRECO(0,0);
    doEffsRECO(1,1);
    //doEffsRECODilep();
  }

  doEffsRECODilep();
  //ISORocCurve();

  //------------------------------------------------------------------------------
  // Set Parameters for other selectors. This is the main point of this selector
  //------------------------------------------------------------------------------

  // SetParam("Muon_4vec",            G_Muon_4vec);
  // SetParam("GEN_PromptMuon_4vec",  G_GEN_PromptMuon_4vec);
  // SetParam("GEN_Muon_4vec",        G_GEN_Muon_4vec);

  // SetParam("MuonID_Tight",         G_MuonID_Tight);
  // SetParam("MuonID_Medium",        G_MuonID_Medium);
  // SetParam("MuonID_HWW",           G_MuonID_HWW);
  // SetParam("MuonID_IPs_HWW",       G_MuonID_IPs_HWW);
  // SetParam("MuonID_GLBorTRKArb",   G_MuonID_GLBorTRKArb);
  // SetParam("MuonID_Fiducial",      G_MuonID_Fiducial);

  // SetParam("MuonISO03",            G_MuonISO03);
  // SetParam("MuonISO03_dBeta",      G_MuonISO03_dBeta);
  // SetParam("MuonISO03_PFWeighted", G_MuonISO03_PFWeighted);
  // SetParam("MuonISO03_PUPPI",      G_MuonISO03_PUPPI);
  // SetParam("MuonISO04",            G_MuonISO04);
  // SetParam("MuonISO04_dBeta",      G_MuonISO04_dBeta);
  // SetParam("MuonISO04_PFWeighted", G_MuonISO04_PFWeighted);
  // SetParam("MuonISO04_PUPPI",      G_MuonISO04_PUPPI);

  // SetParam("Muon_Matching",        G_Muon_Matching);

  // SetParam("RecoMuSize",           G_RecoMuSize);
  // SetParam("NPV",                  G_NPV);

  // SetParam("FLAG_Fiducial",        EvtFlag_Fiducial);
  // SetParam("FLAG_Gen",             EvtFlag_Gen);
  // SetParam("FLAG_Matching",        EvtFlag_Matching);



  //------------------------------------------------------------------------------
  // Do Event Counting
  //------------------------------------------------------------------------------

  if (_Report) Counting();
    
  
} // end inside Loop

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// MEMBER FUNCTIONS
//

//---------------------------------------------------------------------------------------------------------------------
// CheckMuons: Fill bool std::vectors of size G_RecoMuSize to indicate if the RECO muons pass several IDs and ISOs.
//---------------------------------------------------------------------------------------------------------------------
void CoreMuonSelector::CheckMuons() {

  int count = 0; // for counting how many muons pass the fiducial selection
  
  if ( G_RecoMuSize > 0 ) {  // asking for at least one muon in the event 

    for (unsigned int i = 0; i < G_RecoMuSize; ++i) {
      
      // Define selection for several Muon IDs (or components of IDs) and PF Rel. ISOs

      const int NFLAGS = 19;
      bool muon_sel[NFLAGS];   
      for (int j=0; j<NFLAGS; ++j) muon_sel[j] = false;

      muon_sel[0] = (i==0) ? 
	(G_Muon_4vec[i].Pt() > 20. && fabs(G_Muon_4vec[i].Eta()) < 2.4) :
	(G_Muon_4vec[i].Pt() > 10. && fabs(G_Muon_4vec[i].Eta()) < 2.4);
      muon_sel[1] = Get<bool>("T_Muon_IsPFMuon",i);
      muon_sel[2] = Get<bool>("T_Muon_IsGlobalMuon",i);
      muon_sel[3] = Get<bool>("T_Muon_IsTrackerMuon",i) && Get<bool>("T_Muon_IsTrackerMuonArbitrated",i);
      muon_sel[4] = Get<float>("T_Muon_NormChi2GTrk",i) < 10.;
      muon_sel[5] = Get<int>("T_Muon_NValidHitsSATrk",i) > 0;
      muon_sel[6] = Get<int>("T_Muon_NumOfMatchedStations",i) > 1;
      muon_sel[7] = Get<int>("T_Muon_NValidPixelHitsInTrk",i) > 0;
      muon_sel[8] = Get<int>("T_Muon_NLayers",i) > 5;
      muon_sel[9] = Get<float>("T_Muon_IPwrtAveBSInTrack",i) < 0.2;
      muon_sel[10] = fabs(Get<float>("T_Muon_BestTrack_dz",i)) < 0.5;      
      muon_sel[11] = (G_Muon_4vec[i].Pt() < 20   &&  Get<float>("T_Muon_IPwrtAveBSInTrack",i)  < 0.01) ||
	             (G_Muon_4vec[i].Pt() >= 20  &&  Get<float>("T_Muon_IPwrtAveBSInTrack",i)  < 0.02);
      muon_sel[12] = fabs(Get<float>("T_Muon_BestTrack_dz",i)) < 0.1;
      muon_sel[13] = Get<bool>("T_Muon_IsTightMuon",i);
      muon_sel[14] = passMediumID(i);
      muon_sel[14] = Exists("T_Muon_IsMediumMuon") ? Get<bool>("T_Muon_IsMediumMuon",i) : passMediumID(i);
      
      muon_sel[15] = (muon_sel[2] || muon_sel[3]);
      muon_sel[16] = (((muon_sel[2] && muon_sel[4] && muon_sel[5] && muon_sel[6]) || muon_sel[3]) &&
		     muon_sel[7] && muon_sel[8] && muon_sel[11] && muon_sel[12]);
      muon_sel[17] = (((muon_sel[2] && muon_sel[4] && muon_sel[5] && muon_sel[6]) || muon_sel[3]) &&
		     muon_sel[7] && muon_sel[8] && muon_sel[9] && muon_sel[10]);
      muon_sel[18] = (i==0) ? true : (Get<int>("T_Muon_Charge",0)*Get<int>("T_Muon_Charge",i)) < 0;

      
      G_MuonID_Tight.push_back(muon_sel[13]);
      G_MuonID_Medium.push_back(muon_sel[14]);
      G_MuonID_HWW.push_back(muon_sel[16]);
      G_MuonID_Tight_GoT.push_back(muon_sel[17]);
      G_MuonID_IPs_HWW.push_back(muon_sel[11] * muon_sel[12]);
      G_MuonID_GLBorTRKArb.push_back(muon_sel[15]);
      G_MuonID_Fiducial.push_back(muon_sel[0]);

      G_Muon_ChCompatible.push_back(muon_sel[18]);

      G_MuonISO03.push_back(           passISO(i, "R03",           0.12));
      G_MuonISO03_dBeta.push_back(     passISO(i, "dBetaR03",      0.12));
      G_MuonISO03_PFWeighted.push_back(passISO(i, "PFWeightedR03", 0.12));
      G_MuonISO03_PUPPI.push_back(     passISO(i, "PUPPIR03",      0.12));
      G_MuonISO04.push_back(           passISO(i, "R04",           0.12));
      G_MuonISO04_dBeta.push_back(     passISO(i, "dBetaR04",      0.12));
      G_MuonISO04_PFWeighted.push_back(passISO(i, "PFWeightedR04", 0.12));
      G_MuonISO04_PUPPI.push_back(     passISO(i, "PUPPIR04",      0.12));

      if (muon_sel[0]) count++;

      
    } // end loop on muons

    if (_Debug) std::cout << "[DEBUG] Got " << G_RecoMuSize <<" RECO muon(s)" << std::endl;

    if (_Debug) std::cout << "[DEBUG] Got " << count <<" RECO fiducial muon(s)" << std::endl;

  } // end loop on at least one muon
  
}

//---------------------------------------------------------------------------------------------------------------------
// passMediumID: return true if the RECO muon of index 'iMu' passes the Medium ID
//---------------------------------------------------------------------------------------------------------------------
bool CoreMuonSelector::passMediumID(int iMu) {

  bool isMuonID = false;
  bool goodGLB = false;

  goodGLB = Get<bool>("T_Muon_IsGlobalMuon",iMu)      && 
    Get<float>("T_Muon_NormChi2GTrk",iMu) < 3.        &&
    Get<float>("T_Muon_StaTrkChi2LocalPos",iMu) < 12. &&
    Get<float>("T_Muon_trkKink",iMu) < 20.;

  isMuonID = Get<float>("T_Muon_ValidFractionInTrk",iMu) >= 0.8 && Get<bool>("T_Muon_IsPFMuon",iMu) &&
    Get<float>("T_Muon_SegmentCompatibility",iMu) >= (goodGLB ? 0.303 : 0.451);

  return isMuonID;

}

//---------------------------------------------------------------------------------------------------------------------
// passISO: return true if the RECO muon of index 'iMu' passes the PF Relative Isolation indicated by 'typeIso'
//          with a working point 'wp', with the help of getISO() member function
//---------------------------------------------------------------------------------------------------------------------
bool CoreMuonSelector::passISO(int iMu, string typeIso, float wp) {

  bool passIso = false;

  float PFRelIsoBeta = getISO(iMu, typeIso);
     
  if (PFRelIsoBeta <=  wp)  passIso = true;	  

  return passIso;

}

//---------------------------------------------------------------------------------------------------------------------
// getISO: return, for the RECO muon of index 'iMu', the PF Relative Isolation indicated by 'typeIso' 
//---------------------------------------------------------------------------------------------------------------------
float CoreMuonSelector::getISO(int iMu, string typeIso) {
  
  float PFRelIso = 999.9;
  float pt = Get<float>("T_Muon_Pt",iMu);

  if (typeIso == "R03") // PF Rel. ISO, dR=0.3
    PFRelIso = ( Get<float>("T_Muon_chargedHadronIsoR03",iMu) + 
		 Get<float>("T_Muon_neutralHadronIsoR03",iMu) + 
		 Get<float>("T_Muon_photonIsoR03",iMu) )
      / pt;

  else if (typeIso == "R04") // PF Rel. ISO, dR=0.4
    PFRelIso = ( Get<float>("T_Muon_chargedHadronIsoR04",iMu) + 
		 Get<float>("T_Muon_neutralHadronIsoR04",iMu) + 
		 Get<float>("T_Muon_photonIsoR04",iMu) )
      / pt;

  else if (typeIso == "dBetaR03") // PF Rel. ISO, dR=0.3 and dBeta corrections
    PFRelIso = ( Get<float>("T_Muon_chargedHadronIsoR03",iMu) + max(0., Get<float>("T_Muon_neutralHadronIsoR03",iMu) + 
								        Get<float>("T_Muon_photonIsoR03",iMu) - 
								        0.5 * Get<float>("T_Muon_sumPUPtR03",iMu)) )
      / pt;

  else if (typeIso == "dBetaR04") // PF Rel. ISO, dR=0.4 and dBeta corrections 
    PFRelIso = ( Get<float>("T_Muon_chargedHadronIsoR04",iMu) + max(0., Get<float>("T_Muon_neutralHadronIsoR04",iMu) + 
								        Get<float>("T_Muon_photonIsoR04",iMu) - 
								        0.5 * Get<float>("T_Muon_sumPUPtR04",iMu)) )
      / pt;

  else if (typeIso == "PFWeightedR03") // PF Rel. ISO, dR=0.3 and PF-weighted corrections
    PFRelIso = ( Get<float>("T_Muon_chargedHadronIsoR03",iMu) + Get<float>("T_Muon_neutralIsoPFweightR03",iMu) )
      / pt;

  else if (typeIso == "PFWeightedR04") // PF Rel. ISO, dR=0.4 and PF-weighted corrections
    PFRelIso = ( Get<float>("T_Muon_chargedHadronIsoR04",iMu) + Get<float>("T_Muon_neutralIsoPFweightR04",iMu) )
      / pt;

  else if (Exists("T_Muon_neutralIsoPUPPIR03") && typeIso == "PUPPIR03") // PF Rel. ISO, dR=0.3 and PUPPI corr.
    PFRelIso = ( Get<float>("T_Muon_chargedHadronIsoR03",iMu) + Get<float>("T_Muon_neutralIsoPUPPIR03",iMu) )
      / pt;

  else if (Exists("T_Muon_neutralIsoPUPPIR04") && typeIso == "PUPPIR04") // PF Rel. ISO, dR=0.4 and PUPPI corr.
    PFRelIso = ( Get<float>("T_Muon_chargedHadronIsoR04",iMu) + Get<float>("T_Muon_neutralIsoPUPPIR04",iMu) )
      / pt;

  return PFRelIso;

}


// --------------------------------------------------------------------------------------------------------------------
// SetGenInfo: Retrieve all GEN prompt muons, or muons coming from a prompt tau decay, ordered by Pt.
//             Also indicate if there are 2 GEN prompt muons found.
// --------------------------------------------------------------------------------------------------------------------
void CoreMuonSelector::SetGenInfo() {
  
  UInt_t genPromptMuSize = 0;
  //genPromptMuSize = Get<std::vector<float>*>("T_Gen_PromptMuon_Px")->size();
  genPromptMuSize = GetSizeOf("T_Gen_PromptMuon_Px");
  
  UInt_t genPromptTauSize = 0;
  //genPromptTauSize = Get<std::vector<float>*>("T_Gen_PromptTau_Px")->size();
  genPromptTauSize = GetSizeOf("T_Gen_PromptTau_Px");
  
  TLorentzVector p1 = TLorentzVector(0,0,0,0);
  TLorentzVector p2 = TLorentzVector(0,0,0,0);


  if ( _Signal.Contains("Wjets"))
  {
     
    UInt_t genNonPromptMuSize = 0;
    //genNonPromptMuSize = Get<std::vector<float>*>("T_Gen_Muon_Px")->size();
    genNonPromptMuSize = GetSizeOf("T_Gen_Muon_Px");
    
    if ( genPromptMuSize == 1 && fabs(Get<int>("T_Gen_PromptMuon_MpdgId",0)) == 24) 
      
      G_GEN_isMuMu = true; 
    
    if ( genPromptMuSize < 1 && genPromptTauSize == 1 && 
  	 fabs(Get<int>("T_Gen_PromptTau_MpdgId",0))== 24 && 
  	 fabs(Get<int>("T_Gen_PromptTau_LepDec_pdgId",0)) == 13) 

      G_GEN_isTauMu = true; 
    
    
    if ( G_GEN_isMuMu ) {
      G_GEN_PromptMuon_4vec.push_back(TLorentzVector(Get<float>("T_Gen_PromptMuon_Px",0), 
  						     Get<float>("T_Gen_PromptMuon_Py",0),
  						     Get<float>("T_Gen_PromptMuon_Pz",0), 
  						     Get<float>("T_Gen_PromptMuon_Energy",0)));
      if ( genNonPromptMuSize > 0) {
  	G_GEN_isNonPrompt = true;
  	G_GEN_Muon_4vec.push_back(TLorentzVector(Get<float>("T_Gen_Muon_Px",0), 
  						 Get<float>("T_Gen_Muon_Py",0),
  						 Get<float>("T_Gen_Muon_Pz",0), 
  						 Get<float>("T_Gen_Muon_Energy",0)));
      }
    } 
    
    if ( G_GEN_isTauMu ) {
      G_GEN_PromptMuon_4vec.push_back(TLorentzVector(Get<float>("T_Gen_PromptTau_LepDec_Px",0),
  						     Get<float>("T_Gen_PromptTau_LepDec_Py",0),
  						     Get<float>("T_Gen_PromptTau_LepDec_Pz",0), 
  						     Get<float>("T_Gen_PromptTau_LepDec_Energy",0)));
      if ( genNonPromptMuSize > 1) {
  	G_GEN_isNonPrompt = true;
  	G_GEN_Muon_4vec.push_back(TLorentzVector(Get<float>("T_Gen_Muon_Px",1), 
  						 Get<float>("T_Gen_Muon_Py",1),
  						 Get<float>("T_Gen_Muon_Pz",1), 
  						 Get<float>("T_Gen_Muon_Energy",1)));
      }
    }
    
  }
    
  if (_Signal.Contains("GGHWW") || _Signal.Contains("TTbar") || _Signal.Contains("DY"))
    {

      int bosonPdgId = 0;
      if      (_Signal.Contains("GGHWW") || _Signal.Contains("TTbar"))  bosonPdgId = 24;
      else if (_Signal.Contains("DY"))                                  bosonPdgId = 23;        
      
      if ( genPromptMuSize == 2 && fabs(Get<int>("T_Gen_PromptMuon_MpdgId",0)) == bosonPdgId && 
  	   fabs(Get<int>("T_Gen_PromptMuon_MpdgId",1)) == bosonPdgId &&
  	   (Get<int>("T_Gen_PromptMuon_pdgId",0)*Get<int>("T_Gen_PromptMuon_pdgId",1)) < 0) {

  	p1 = TLorentzVector(Get<float>("T_Gen_PromptMuon_Px",0), 
  			    Get<float>("T_Gen_PromptMuon_Py",0),
  			    Get<float>("T_Gen_PromptMuon_Pz",0), 
  			    Get<float>("T_Gen_PromptMuon_Energy",0));

  	p2 = TLorentzVector(Get<float>("T_Gen_PromptMuon_Px",1), 
  			    Get<float>("T_Gen_PromptMuon_Py",1),
  			    Get<float>("T_Gen_PromptMuon_Pz",1), 
  			    Get<float>("T_Gen_PromptMuon_Energy",1));

  	G_GEN_isMuMu = true; 

      }
      
      if ( genPromptMuSize == 1 && fabs(Get<int>("T_Gen_PromptMuon_MpdgId",0)) == bosonPdgId && 
  	   genPromptTauSize == 1 && fabs(Get<int>("T_Gen_PromptTau_MpdgId",0))== bosonPdgId && 
  	   fabs(Get<int>("T_Gen_PromptTau_LepDec_pdgId",0)) == 13 &&
  	   (Get<int>("T_Gen_PromptMuon_pdgId",0)*Get<int>("T_Gen_PromptTau_LepDec_pdgId",0)) < 0) {

  	p1 = TLorentzVector(Get<float>("T_Gen_PromptMuon_Px",0), 
  			    Get<float>("T_Gen_PromptMuon_Py",0),
  			    Get<float>("T_Gen_PromptMuon_Pz",0), 
  			    Get<float>("T_Gen_PromptMuon_Energy",0));

  	p2 = TLorentzVector(Get<float>("T_Gen_PromptTau_LepDec_Px",0), 
  			    Get<float>("T_Gen_PromptTau_LepDec_Py",0),
  			    Get<float>("T_Gen_PromptTau_LepDec_Pz",0), 
  			    Get<float>("T_Gen_PromptTau_LepDec_Energy",0));

  	if (p1.Pt() >= p2.Pt()) G_GEN_isMuTau = true;
  	else                    G_GEN_isTauMu = true;

      }
      
      if ( genPromptMuSize < 1 && genPromptTauSize == 2 && 
  	   fabs(Get<int>("T_Gen_PromptTau_MpdgId",0))== bosonPdgId && 
  	   fabs(Get<int>("T_Gen_PromptTau_MpdgId",1))== bosonPdgId && 
  	   fabs(Get<int>("T_Gen_PromptTau_LepDec_pdgId",0)) == 13 && 
  	   fabs(Get<int>("T_Gen_PromptTau_LepDec_pdgId",1)) == 13 &&
  	   (Get<int>("T_Gen_PromptTau_LepDec_pdgId",0)*Get<int>("T_Gen_PromptTau_LepDec_pdgId",1)) < 0) {

  	p1 = TLorentzVector(Get<float>("T_Gen_PromptTau_LepDec_Px",0), 
  			    Get<float>("T_Gen_PromptTau_LepDec_Py",0),
  			    Get<float>("T_Gen_PromptTau_LepDec_Pz",0), 
  			    Get<float>("T_Gen_PromptTau_LepDec_Energy",0));

  	p2 = TLorentzVector(Get<float>("T_Gen_PromptTau_LepDec_Px",1), 
  			    Get<float>("T_Gen_PromptTau_LepDec_Py",1),
  			    Get<float>("T_Gen_PromptTau_LepDec_Pz",1), 
  			    Get<float>("T_Gen_PromptTau_LepDec_Energy",1));

  	G_GEN_isTauTau = true; 

      }
      
    }


  G_GEN_Pass = G_GEN_isMuMu || G_GEN_isMuTau || G_GEN_isTauMu || G_GEN_isTauTau;

  if ((_Signal.Contains("DY") || _Signal.Contains("GGHWW") || _Signal.Contains("TTbar")) && G_GEN_Pass) {

    if ( p1.Pt() >= p2.Pt() ) {

      G_GEN_PromptMuon_4vec.push_back(p1);
      G_GEN_PromptMuon_4vec.push_back(p2);

    }

    else {

      G_GEN_PromptMuon_4vec.push_back(p2);
      G_GEN_PromptMuon_4vec.push_back(p1);

    }

    if (_Debug) std::cout << "[DEBUG] Got 2 GEN muons" << std::endl;

  }
        

}


//---------------------------------------------------------------------------------------------------------------------
// GetMatching: fill an int std::vector indicating to which GEN prompt muon the RECO muons are matched
//               * 1: matched to the 1st GEN prompt muon
//               * 2: matched to the 2nd GEN prompt muon
//               * 0: not matched to any GEN prompt muon
//---------------------------------------------------------------------------------------------------------------------
void CoreMuonSelector::GetMatching() {

  UInt_t GenSize = 0;
  GenSize = G_GEN_PromptMuon_4vec.size();

  for (unsigned int i = 0; i < G_RecoMuSize; ++i) {

    int isMatchedTo = 0;

    if (G_MuonID_Fiducial[i] && G_MuonID_GLBorTRKArb[i]) { 

      for (UInt_t j = 0; j < GenSize; ++j) {
      
	Double_t dR  = 999.;
	//Double_t dPt = 999.;

	dR  = G_Muon_4vec[i].DeltaR(G_GEN_PromptMuon_4vec[j]);
	//dPt = fabs(G_Muon_4vec[i].Pt() - G_GEN_PromptMuon_4vec[j].Pt());
      
	if (dR <= 0.01) { // && dPt <= 10
	  isMatchedTo = j+1;
	  break;
	}
      
      }

    }

    G_Muon_Matching.push_back(isMatchedTo);

  }

  int count_matched = 0;
  int count_1 = 0;
  int count_2 = 0;

  for (unsigned int i = 0; i < G_RecoMuSize; ++i) {

    if (G_Muon_Matching[i]) {
      ++count_matched;
      if (G_Muon_Matching[i] == 1) ++count_1;
      else if (G_Muon_Matching[i] == 2) ++count_2;
    }

  }

  if (G_Muon_Matching[0]) {
    if (count_matched == 2) {
      if (G_Muon_Matching[1] && G_Muon_ChCompatible[1]) {
	if (count_1 == 1 && count_2 == 1) {
	  G_PassMatching = true;
	  if (_Debug) std::cout << "[DEBUG] Matching successful" << std::endl;
	}
      }
    }
  }
  
}

//---------------------------------------------------------------------------------------------------------------------
// SetEventFlags: Assign some event flags to get ready to pass them as parameters to other selectors
//---------------------------------------------------------------------------------------------------------------------
void CoreMuonSelector::SetEventFlags() {

  if (G_RecoMuSize >= 2) {
    EvtFlag_Fiducial = G_MuonID_Fiducial[0] && G_MuonID_Fiducial[1];
  }

  EvtFlag_Gen = G_GEN_Pass;

  EvtFlag_Matching = G_PassMatching;
  if (_IsDATA || _Signal.Contains("QCD")) EvtFlag_Matching = EvtFlag_Fiducial;

}

//---------------------------------------------------------------------------------------------------------------------
// Counting: Increments several counting variables to display the number of events fulfilling certain criteria during
//           the Summary method
//---------------------------------------------------------------------------------------------------------------------
void CoreMuonSelector::Counting() {

  ++GCount_AllEvents;
  if (G_GEN_Pass) ++GCount_GenEvents;

  UInt_t n_fiducial = 0;
  UInt_t n_match    = 0;
  UInt_t n_tight    = 0;
  UInt_t n_iso      = 0;
  UInt_t n_mtight   = 0;
  UInt_t n_miso     = 0;

  for (unsigned int i = 0; i < G_RecoMuSize; i++) {

    if (G_MuonID_Fiducial[i]) ++n_fiducial;
    if (G_MuonID_Fiducial[i] && G_Muon_Matching[i])   ++n_match;
    if (G_MuonID_Fiducial[i] && G_MuonID_Tight[i])    ++n_tight;
    if (G_MuonID_Fiducial[i] && G_MuonID_Tight[i] && G_MuonISO04_dBeta[i]) ++n_iso;
    if (G_MuonID_Fiducial[i] && G_Muon_Matching[i] && G_MuonID_Tight[i])    ++n_mtight;
    if (G_MuonID_Fiducial[i] && G_Muon_Matching[i] && G_MuonID_Tight[i] && G_MuonISO04_dBeta[i]) ++n_miso;

  }

  if (G_GEN_Pass) {
    if (G_MuonID_Fiducial[0] && n_fiducial >= 2) {
      GCount_Fiducial_AtLeast2++;
      if (n_fiducial == 2) {
	GCount_Fiducial_2++;
	if (G_MuonID_Fiducial[1]) {
	  ++GCount_Fiducial_1st2nd;
	  if (n_match == 2 && (G_Muon_Matching[0] != G_Muon_Matching[1])) {
	    ++GCount_Match_1st2nd;
	    if (n_mtight == 2) {
	      ++GCount_MatchTight_1st2nd;
	      if (n_miso == 2) ++GCount_MatchTightIso_1st2nd;
	      else if (n_miso < 2) {
		if (n_miso == 1 && G_MuonISO04_dBeta[0])      ++GCount_MatchTightIso_Only1st;
		else if (n_miso == 1 && G_MuonISO04_dBeta[1]) ++GCount_MatchTightIso_Only2nd;
		else if (n_miso == 0)                         ++GCount_MatchTightIso_None;
	      }
	    }
	    else if (n_mtight < 2) {
	      if (n_mtight == 1 && G_MuonID_Tight[0])      ++GCount_MatchTight_Only1st;
	      else if (n_mtight == 1 && G_MuonID_Tight[1]) ++GCount_MatchTight_Only2nd;
	      else if (n_mtight == 0)                      ++GCount_MatchTight_None;
	    }
	  }
	  else if (n_match < 2) {
	    ++GCount_NoMatch_1st2nd;
	    if (n_tight == 2) {
	      ++GCount_Tight_1st2nd;
	      if (n_iso == 2) ++GCount_TightIso_1st2nd;
	      else if (n_iso < 2) {
		if (n_iso == 1 && G_MuonISO04_dBeta[0])      ++GCount_TightIso_Only1st;
		else if (n_iso == 1 && G_MuonISO04_dBeta[1]) ++GCount_TightIso_Only2nd;
		else if (n_iso == 0)                         ++GCount_TightIso_None;
	      }
	    }
	    else if (n_tight < 2) {
	      if (n_tight == 1 && G_MuonID_Tight[0])      ++GCount_Tight_Only1st;
	      else if (n_tight == 1 && G_MuonID_Tight[1]) ++GCount_Tight_Only2nd;
	      else if (n_tight == 0)                      ++GCount_Tight_None;
	    }
	  }
	}
	else if (G_MuonID_Fiducial[2]) {
	  ++GCount_Fiducial_1st3rd;
	  if (n_match == 2 && (G_Muon_Matching[0] != G_Muon_Matching[2])) {
	    ++GCount_Match_1st3rd;
	    if (n_mtight == 2) {
	      ++GCount_MatchTight_1st3rd;
	      if (n_miso == 2) ++GCount_MatchTightIso_1st3rd;
	      else if (n_miso < 2) {
		if (n_miso == 1 && G_MuonISO04_dBeta[0])      ++GCount_MatchTightIso_Only1st;
		else if (n_miso == 1 && G_MuonISO04_dBeta[2]) ++GCount_MatchTightIso_Only3rd;
		else if (n_miso == 0)                         ++GCount_MatchTightIso_None;
	      }
	    }
	    else if (n_mtight < 2) {
	      if (n_mtight == 1 && G_MuonID_Tight[0])      ++GCount_MatchTight_Only1st;
	      else if (n_mtight == 1 && G_MuonID_Tight[2]) ++GCount_MatchTight_Only3rd;
	      else if (n_mtight == 0)                      ++GCount_MatchTight_None;
	    }
	  }
	  if (n_match < 2) {
	    ++GCount_NoMatch_1st3rd;
	    if (n_tight == 2) {
	      ++GCount_Tight_1st3rd;
	      if (n_iso == 2) ++GCount_TightIso_1st3rd;
	      else if (n_iso < 2) {
		if (n_iso == 1 && G_MuonISO04_dBeta[0])      ++GCount_TightIso_Only1st;
		else if (n_iso == 1 && G_MuonISO04_dBeta[2]) ++GCount_TightIso_Only3rd;
		else if (n_iso == 0)                         ++GCount_TightIso_None;
	      }
	    }
	    else if (n_tight < 2) {
	      if (n_tight == 1 && G_MuonID_Tight[0])      ++GCount_Tight_Only1st;
	      else if (n_tight == 1 && G_MuonID_Tight[2]) ++GCount_Tight_Only3rd;
	      else if (n_tight == 0)                      ++GCount_Tight_None;
	    }
	  }
	}
	else {
	  ++GCount_Fiducial_1stOther;
	}
      }
      else if (n_fiducial > 2) {
	++GCount_Fiducial_MoreThan2;
	if (G_Muon_Matching[0] && n_match == 2) ++GCount_Match_MoreThan2_OK;
	if (G_Muon_Matching[0] && n_match > 2) ++GCount_Match_MoreThan2;
	if (G_Muon_Matching[0] && n_match < 2) ++GCount_NoMatch_MoreThan2;
	if      (n_fiducial == 3) ++GCount_Fiducial_3;
	else if (n_fiducial > 3)  ++GCount_Fiducial_MoreThan3;
      }
    }
    else if (G_MuonID_Fiducial[0] && n_fiducial == 1) {
      ++GCount_Fiducial_Only1st;
    }
    else if (n_fiducial == 0) {
      ++GCount_Fiducial_None;
    }
  }	  
      
}

void CoreMuonSelector::doEffsRECO(int iMu, int indexMuon) {

  float pt  = Get<float>("T_Muon_Pt", iMu); 
  float eta = Get<float>("T_Muon_Eta",iMu);
  float npv = G_NPV;

  
  h_Eff_pt_Matched[indexMuon] ->Fill(pt);
  h_Eff_eta_Matched[indexMuon]->Fill(eta);
  h_Eff_npv_Matched[indexMuon]->Fill(npv);

  if (G_MuonID_Tight[iMu]) {
    h_Eff_pt_TightID[indexMuon] ->Fill(pt);
    h_Eff_eta_TightID[indexMuon]->Fill(eta);
    h_Eff_npv_TightID[indexMuon]->Fill(npv);
  }

  if (G_MuonID_Medium[iMu]) {
    h_Eff_pt_MediumID[indexMuon] ->Fill(pt);
    h_Eff_eta_MediumID[indexMuon]->Fill(eta);
    h_Eff_npv_MediumID[indexMuon]->Fill(npv);
  }
  
  if (G_MuonID_HWW[iMu]) {
    h_Eff_pt_HWWID[indexMuon] ->Fill(pt);
    h_Eff_eta_HWWID[indexMuon]->Fill(eta);
    h_Eff_npv_HWWID[indexMuon]->Fill(npv);
  }

  if (G_MuonID_Tight_GoT[iMu]) {
    h_Eff_pt_TightIDGoT[indexMuon] ->Fill(pt);
    h_Eff_eta_TightIDGoT[indexMuon]->Fill(eta);
    h_Eff_npv_TightIDGoT[indexMuon]->Fill(npv);
  }
      
  if (G_MuonID_Tight[iMu] && G_MuonID_IPs_HWW[iMu]) {
    h_Eff_pt_TightIDipsHWW[indexMuon] ->Fill(pt);
    h_Eff_eta_TightIDipsHWW[indexMuon]->Fill(eta);
    h_Eff_npv_TightIDipsHWW[indexMuon]->Fill(npv);
  }

  if (G_MuonID_Medium[iMu] && G_MuonID_IPs_HWW[iMu]) {
    h_Eff_pt_MediumIDipsHWW[indexMuon] ->Fill(pt);
    h_Eff_eta_MediumIDipsHWW[indexMuon]->Fill(eta);
    h_Eff_npv_MediumIDipsHWW[indexMuon]->Fill(npv);
  }

  if (G_MuonID_Tight[iMu] && G_MuonISO03[iMu]) {
    h_Eff_pt_TightID_ISO03[indexMuon] ->Fill(pt);
    h_Eff_eta_TightID_ISO03[indexMuon]->Fill(eta);
    h_Eff_npv_TightID_ISO03[indexMuon]->Fill(npv);
  }
  
  if (G_MuonID_Tight[iMu] && G_MuonISO04[iMu]) {
    h_Eff_pt_TightID_ISO04[indexMuon] ->Fill(pt);
    h_Eff_eta_TightID_ISO04[indexMuon]->Fill(eta);
    h_Eff_npv_TightID_ISO04[indexMuon]->Fill(npv);
  }
  
  if (G_MuonID_Tight[iMu] && G_MuonISO03_dBeta[iMu]) {
    h_Eff_pt_TightID_ISO03dBeta[indexMuon] ->Fill(pt);
    h_Eff_eta_TightID_ISO03dBeta[indexMuon]->Fill(eta);
    h_Eff_npv_TightID_ISO03dBeta[indexMuon]->Fill(npv);
  }
  
  if (G_MuonID_Tight[iMu] && G_MuonISO04_dBeta[iMu]) {
    h_Eff_pt_TightID_ISO04dBeta[indexMuon] ->Fill(pt);
    h_Eff_eta_TightID_ISO04dBeta[indexMuon]->Fill(eta);
    h_Eff_npv_TightID_ISO04dBeta[indexMuon]->Fill(npv);
  }
  
  if (G_MuonID_Tight[iMu] && G_MuonISO03_PFWeighted[iMu]) {
    h_Eff_pt_TightID_ISO03PFWeighted[indexMuon] ->Fill(pt);
    h_Eff_eta_TightID_ISO03PFWeighted[indexMuon]->Fill(eta);
    h_Eff_npv_TightID_ISO03PFWeighted[indexMuon]->Fill(npv);
  }
  
  if (G_MuonID_Tight[iMu] && G_MuonISO04_PFWeighted[iMu]) {
    h_Eff_pt_TightID_ISO04PFWeighted[indexMuon] ->Fill(pt);
    h_Eff_eta_TightID_ISO04PFWeighted[indexMuon]->Fill(eta);
    h_Eff_npv_TightID_ISO04PFWeighted[indexMuon]->Fill(npv);
  } 

  if (G_MuonID_Tight[iMu] && G_MuonISO03_PUPPI[iMu]) {
    h_Eff_pt_TightID_ISO03PUPPI[indexMuon] ->Fill(pt);
    h_Eff_eta_TightID_ISO03PUPPI[indexMuon]->Fill(eta);
    h_Eff_npv_TightID_ISO03PUPPI[indexMuon]->Fill(npv);
  }
  
  if (G_MuonID_Tight[iMu] && G_MuonISO04_PUPPI[iMu]) {
    h_Eff_pt_TightID_ISO04PUPPI[indexMuon] ->Fill(pt);
    h_Eff_eta_TightID_ISO04PUPPI[indexMuon]->Fill(eta);
    h_Eff_npv_TightID_ISO04PUPPI[indexMuon]->Fill(npv);
  } 

  if (G_MuonID_Medium[iMu] && G_MuonISO03[iMu]) {
    h_Eff_pt_MediumID_ISO03[indexMuon] ->Fill(pt);
    h_Eff_eta_MediumID_ISO03[indexMuon]->Fill(eta);
    h_Eff_npv_MediumID_ISO03[indexMuon]->Fill(npv);
  }
  
  if (G_MuonID_Medium[iMu] && G_MuonISO04[iMu]) {
    h_Eff_pt_MediumID_ISO04[indexMuon] ->Fill(pt);
    h_Eff_eta_MediumID_ISO04[indexMuon]->Fill(eta);
    h_Eff_npv_MediumID_ISO04[indexMuon]->Fill(npv);
  }
  
  if (G_MuonID_Medium[iMu] && G_MuonISO03_dBeta[iMu]) {
    h_Eff_pt_MediumID_ISO03dBeta[indexMuon] ->Fill(pt);
    h_Eff_eta_MediumID_ISO03dBeta[indexMuon]->Fill(eta);
    h_Eff_npv_MediumID_ISO03dBeta[indexMuon]->Fill(npv);
  }
  
  if (G_MuonID_Medium[iMu] && G_MuonISO04_dBeta[iMu]) {
    h_Eff_pt_MediumID_ISO04dBeta[indexMuon] ->Fill(pt);
    h_Eff_eta_MediumID_ISO04dBeta[indexMuon]->Fill(eta);
    h_Eff_npv_MediumID_ISO04dBeta[indexMuon]->Fill(npv);
  }
  
  if (G_MuonID_Medium[iMu] && G_MuonISO03_PFWeighted[iMu]) {
    h_Eff_pt_MediumID_ISO03PFWeighted[indexMuon] ->Fill(pt);
    h_Eff_eta_MediumID_ISO03PFWeighted[indexMuon]->Fill(eta);
    h_Eff_npv_MediumID_ISO03PFWeighted[indexMuon]->Fill(npv);
  }
  
  if (G_MuonID_Medium[iMu] && G_MuonISO04_PFWeighted[iMu]) {
    h_Eff_pt_MediumID_ISO04PFWeighted[indexMuon] ->Fill(pt);
    h_Eff_eta_MediumID_ISO04PFWeighted[indexMuon]->Fill(eta);
    h_Eff_npv_MediumID_ISO04PFWeighted[indexMuon]->Fill(npv);
  } 

  if (G_MuonID_Medium[iMu] && G_MuonISO03_PUPPI[iMu]) {
    h_Eff_pt_MediumID_ISO03PUPPI[indexMuon] ->Fill(pt);
    h_Eff_eta_MediumID_ISO03PUPPI[indexMuon]->Fill(eta);
    h_Eff_npv_MediumID_ISO03PUPPI[indexMuon]->Fill(npv);
  }
  
  if (G_MuonID_Medium[iMu] && G_MuonISO04_PUPPI[iMu]) {
    h_Eff_pt_MediumID_ISO04PUPPI[indexMuon] ->Fill(pt);
    h_Eff_eta_MediumID_ISO04PUPPI[indexMuon]->Fill(eta);
    h_Eff_npv_MediumID_ISO04PUPPI[indexMuon]->Fill(npv);
  } 
  
}


void CoreMuonSelector::doEffsRECODilep() {

  for (UInt_t i=0; i<G_RecoMuSize; i++) {

    float pt  = Get<float>("T_Muon_Pt", i); 
    float eta = Get<float>("T_Muon_Eta",i);
    float npv = G_NPV;

    bool tight  = G_MuonID_Tight_GoT[i];
    bool medium = G_MuonID_Medium[i] && G_MuonID_IPs_HWW[i];

    if (pt < 20 && fabs(eta) > 2.4) continue;
    if (_Signal.Contains("DY") && !G_Muon_Matching[i]) continue;
    
    if (tight) {
       h_Eff_pt_TightID_Dilep->Fill(pt);
      h_Eff_eta_TightID_Dilep->Fill(eta);
      h_Eff_npv_TightID_Dilep->Fill(npv);
    }
    
    if (medium) {
       h_Eff_pt_MediumID_Dilep->Fill(pt);
      h_Eff_eta_MediumID_Dilep->Fill(eta);
      h_Eff_npv_MediumID_Dilep->Fill(npv);
    }
    
    if (tight && G_MuonISO03[i]) {
       h_Eff_pt_TightID_ISO03_Dilep->Fill(pt);
      h_Eff_eta_TightID_ISO03_Dilep->Fill(eta);
      h_Eff_npv_TightID_ISO03_Dilep->Fill(npv);
    }
  
    if (tight && G_MuonISO04[i]) {
       h_Eff_pt_TightID_ISO04_Dilep->Fill(pt);
      h_Eff_eta_TightID_ISO04_Dilep->Fill(eta);
      h_Eff_npv_TightID_ISO04_Dilep->Fill(npv);
    }
  
    if (tight && G_MuonISO03_dBeta[i]) {
       h_Eff_pt_TightID_ISO03dBeta_Dilep->Fill(pt);
      h_Eff_eta_TightID_ISO03dBeta_Dilep->Fill(eta);
      h_Eff_npv_TightID_ISO03dBeta_Dilep->Fill(npv);
    }
  
    if (tight && G_MuonISO04_dBeta[i]) {
       h_Eff_pt_TightID_ISO04dBeta_Dilep->Fill(pt);
      h_Eff_eta_TightID_ISO04dBeta_Dilep->Fill(eta);
      h_Eff_npv_TightID_ISO04dBeta_Dilep->Fill(npv);
    }
  
    if (tight && G_MuonISO03_PFWeighted[i]) {
       h_Eff_pt_TightID_ISO03PFWeighted_Dilep->Fill(pt);
      h_Eff_eta_TightID_ISO03PFWeighted_Dilep->Fill(eta);
      h_Eff_npv_TightID_ISO03PFWeighted_Dilep->Fill(npv);
    }
  
    if (tight && G_MuonISO04_PFWeighted[i]) {
       h_Eff_pt_TightID_ISO04PFWeighted_Dilep->Fill(pt);
      h_Eff_eta_TightID_ISO04PFWeighted_Dilep->Fill(eta);
      h_Eff_npv_TightID_ISO04PFWeighted_Dilep->Fill(npv);
    } 

    if (tight && G_MuonISO03_PUPPI[i]) {
       h_Eff_pt_TightID_ISO03PUPPI_Dilep->Fill(pt);
      h_Eff_eta_TightID_ISO03PUPPI_Dilep->Fill(eta);
      h_Eff_npv_TightID_ISO03PUPPI_Dilep->Fill(npv);
    }
  
    if (tight && G_MuonISO04_PUPPI[i]) {
       h_Eff_pt_TightID_ISO04PUPPI_Dilep->Fill(pt);
      h_Eff_eta_TightID_ISO04PUPPI_Dilep->Fill(eta);
      h_Eff_npv_TightID_ISO04PUPPI_Dilep->Fill(npv);
    } 

    if (medium && G_MuonISO03[i]) {
       h_Eff_pt_MediumID_ISO03_Dilep->Fill(pt);
      h_Eff_eta_MediumID_ISO03_Dilep->Fill(eta);
      h_Eff_npv_MediumID_ISO03_Dilep->Fill(npv);
    }
  
    if (medium && G_MuonISO04[i]) {
       h_Eff_pt_MediumID_ISO04_Dilep->Fill(pt);
      h_Eff_eta_MediumID_ISO04_Dilep->Fill(eta);
      h_Eff_npv_MediumID_ISO04_Dilep->Fill(npv);
    }
  
    if (medium && G_MuonISO03_dBeta[i]) {
       h_Eff_pt_MediumID_ISO03dBeta_Dilep->Fill(pt);
      h_Eff_eta_MediumID_ISO03dBeta_Dilep->Fill(eta);
      h_Eff_npv_MediumID_ISO03dBeta_Dilep->Fill(npv);
    }
  
    if (medium && G_MuonISO04_dBeta[i]) {
       h_Eff_pt_MediumID_ISO04dBeta_Dilep->Fill(pt);
      h_Eff_eta_MediumID_ISO04dBeta_Dilep->Fill(eta);
      h_Eff_npv_MediumID_ISO04dBeta_Dilep->Fill(npv);
    }
  
    if (medium && G_MuonISO03_PFWeighted[i]) {
       h_Eff_pt_MediumID_ISO03PFWeighted_Dilep->Fill(pt);
      h_Eff_eta_MediumID_ISO03PFWeighted_Dilep->Fill(eta);
      h_Eff_npv_MediumID_ISO03PFWeighted_Dilep->Fill(npv);
    }
  
    if (medium && G_MuonISO04_PFWeighted[i]) {
       h_Eff_pt_MediumID_ISO04PFWeighted_Dilep->Fill(pt);
      h_Eff_eta_MediumID_ISO04PFWeighted_Dilep->Fill(eta);
      h_Eff_npv_MediumID_ISO04PFWeighted_Dilep->Fill(npv);
    } 

    if (medium && G_MuonISO03_PUPPI[i]) {
       h_Eff_pt_MediumID_ISO03PUPPI_Dilep->Fill(pt);
      h_Eff_eta_MediumID_ISO03PUPPI_Dilep->Fill(eta);
      h_Eff_npv_MediumID_ISO03PUPPI_Dilep->Fill(npv);
    }
  
    if (medium && G_MuonISO04_PUPPI[i]) {
       h_Eff_pt_MediumID_ISO04PUPPI_Dilep->Fill(pt);
      h_Eff_eta_MediumID_ISO04PUPPI_Dilep->Fill(eta);
      h_Eff_npv_MediumID_ISO04PUPPI_Dilep->Fill(npv);
    } 

  }
  
}

void CoreMuonSelector::ISORocCurve() {

  for (UInt_t i=0; i<G_RecoMuSize; i++) {

    float pt  = Get<float>("T_Muon_Pt", i); 
    float eta = Get<float>("T_Muon_Eta",i);

    if (pt < 20 && fabs(eta) > 2.4) continue;
    if (_Signal.Contains("DY") && !G_Muon_Matching[i]) continue;

    if (G_MuonID_Tight[i]) {

      h_RC_TightID_ISO03_Dilep          ->Fill(-1);
      h_RC_TightID_ISO04_Dilep          ->Fill(-1);
      h_RC_TightID_ISO03dBeta_Dilep     ->Fill(-1);
      h_RC_TightID_ISO04dBeta_Dilep     ->Fill(-1);
      h_RC_TightID_ISO03PFWeighted_Dilep->Fill(-1);
      h_RC_TightID_ISO04PFWeighted_Dilep->Fill(-1);
      h_RC_TightID_ISO03PUPPI_Dilep     ->Fill(-1);
      h_RC_TightID_ISO04PUPPI_Dilep     ->Fill(-1);

      for (int j = 0; j < 50; ++j) {

	if (passISO(i, "R03",           j*1.0/100.0)) h_RC_TightID_ISO03_Dilep          ->Fill(j);
	if (passISO(i, "dBetaR03",      j*1.0/100.0)) h_RC_TightID_ISO03dBeta_Dilep     ->Fill(j);
	if (passISO(i, "PFWeightedR03", j*1.0/100.0)) h_RC_TightID_ISO03PFWeighted_Dilep->Fill(j);
	if (passISO(i, "PUPPIR03",      j*1.0/100.0)) h_RC_TightID_ISO03PUPPI_Dilep     ->Fill(j);
	if (passISO(i, "R04",           j*1.0/100.0)) h_RC_TightID_ISO04_Dilep          ->Fill(j);
	if (passISO(i, "dBetaR04",      j*1.0/100.0)) h_RC_TightID_ISO04dBeta_Dilep     ->Fill(j);
	if (passISO(i, "PFWeightedR04", j*1.0/100.0)) h_RC_TightID_ISO04PFWeighted_Dilep->Fill(j);
	if (passISO(i, "PUPPIR04",      j*1.0/100.0)) h_RC_TightID_ISO04PUPPI_Dilep     ->Fill(j);

      }

    }

    if (G_MuonID_Medium[i]) {

      h_RC_MediumID_ISO03_Dilep          ->Fill(-1);
      h_RC_MediumID_ISO04_Dilep          ->Fill(-1);
      h_RC_MediumID_ISO03dBeta_Dilep     ->Fill(-1);
      h_RC_MediumID_ISO04dBeta_Dilep     ->Fill(-1);
      h_RC_MediumID_ISO03PFWeighted_Dilep->Fill(-1);
      h_RC_MediumID_ISO04PFWeighted_Dilep->Fill(-1);
      h_RC_MediumID_ISO03PUPPI_Dilep     ->Fill(-1);
      h_RC_MediumID_ISO04PUPPI_Dilep     ->Fill(-1);

      for (int j = 0; j < 50; ++j) {

	if (passISO(i, "R03",           j*1.0/100.0)) h_RC_MediumID_ISO03_Dilep          ->Fill(j);
	if (passISO(i, "dBetaR03",      j*1.0/100.0)) h_RC_MediumID_ISO03dBeta_Dilep     ->Fill(j);
	if (passISO(i, "PFWeightedR03", j*1.0/100.0)) h_RC_MediumID_ISO03PFWeighted_Dilep->Fill(j);
	if (passISO(i, "PUPPIR03",      j*1.0/100.0)) h_RC_MediumID_ISO03PUPPI_Dilep     ->Fill(j);
	if (passISO(i, "R04",           j*1.0/100.0)) h_RC_MediumID_ISO04_Dilep          ->Fill(j);
	if (passISO(i, "dBetaR04",      j*1.0/100.0)) h_RC_MediumID_ISO04dBeta_Dilep     ->Fill(j);
	if (passISO(i, "PFWeightedR04", j*1.0/100.0)) h_RC_MediumID_ISO04PFWeighted_Dilep->Fill(j);
	if (passISO(i, "PUPPIR04",      j*1.0/100.0)) h_RC_MediumID_ISO04PUPPI_Dilep     ->Fill(j);

      }

    }

  }
  
}

void CoreMuonSelector::Summary() {
  // Get Data Members at the client-master (after finishing the analysis at the workers nodes)
  // Only data members set here will be accesible at the client-master


  // 1D histos  

  h_N_PV = FindOutput<TH1F*>("h_N_PV");

  // Single Muon ID and ISO efficiencies vs pt, eta and npv
  h_Eff_pt_Matched[0]                   = FindOutput<TH1F*>("h_Eff_pt_Matched_Mu1");
  h_Eff_pt_Matched[1]                   = FindOutput<TH1F*>("h_Eff_pt_Matched_Mu2");
  h_Eff_pt_TightID[0]                   = FindOutput<TH1F*>("h_Eff_pt_TightID_Mu1");
  h_Eff_pt_TightID[1]                   = FindOutput<TH1F*>("h_Eff_pt_TightID_Mu2");
  h_Eff_pt_MediumID[0]                  = FindOutput<TH1F*>("h_Eff_pt_MediumID_Mu1");
  h_Eff_pt_MediumID[1]                  = FindOutput<TH1F*>("h_Eff_pt_MediumID_Mu2");
  h_Eff_pt_HWWID[0]                     = FindOutput<TH1F*>("h_Eff_pt_HWWID_Mu1");
  h_Eff_pt_HWWID[1]                     = FindOutput<TH1F*>("h_Eff_pt_HWWID_Mu2");
  h_Eff_pt_TightIDGoT[0]                = FindOutput<TH1F*>("h_Eff_pt_TightIDGoT_Mu1");
  h_Eff_pt_TightIDGoT[1]                = FindOutput<TH1F*>("h_Eff_pt_TightIDGoT_Mu2");
  h_Eff_pt_TightIDipsHWW[0]             = FindOutput<TH1F*>("h_Eff_pt_TightIDipsHWW_Mu1");
  h_Eff_pt_TightIDipsHWW[1]             = FindOutput<TH1F*>("h_Eff_pt_TightIDipsHWW_Mu2");
  h_Eff_pt_MediumIDipsHWW[0]            = FindOutput<TH1F*>("h_Eff_pt_MediumIDipsHWW_Mu1");
  h_Eff_pt_MediumIDipsHWW[1]            = FindOutput<TH1F*>("h_Eff_pt_MediumIDipsHWW_Mu2");
  h_Eff_pt_TightID_ISO03[0]             = FindOutput<TH1F*>("h_Eff_pt_TightID_ISO03_Mu1");
  h_Eff_pt_TightID_ISO03[1]             = FindOutput<TH1F*>("h_Eff_pt_TightID_ISO03_Mu2");
  h_Eff_pt_TightID_ISO04[0]             = FindOutput<TH1F*>("h_Eff_pt_TightID_ISO04_Mu1");
  h_Eff_pt_TightID_ISO04[1]             = FindOutput<TH1F*>("h_Eff_pt_TightID_ISO04_Mu2");
  h_Eff_pt_TightID_ISO03dBeta[0]        = FindOutput<TH1F*>("h_Eff_pt_TightID_ISO03dBeta_Mu1");
  h_Eff_pt_TightID_ISO03dBeta[1]        = FindOutput<TH1F*>("h_Eff_pt_TightID_ISO03dBeta_Mu2");
  h_Eff_pt_TightID_ISO04dBeta[0]        = FindOutput<TH1F*>("h_Eff_pt_TightID_ISO04dBeta_Mu1");
  h_Eff_pt_TightID_ISO04dBeta[1]        = FindOutput<TH1F*>("h_Eff_pt_TightID_ISO04dBeta_Mu2");
  h_Eff_pt_TightID_ISO03PFWeighted[0]   = FindOutput<TH1F*>("h_Eff_pt_TightID_ISO03PFWeighted_Mu1");
  h_Eff_pt_TightID_ISO03PFWeighted[1]   = FindOutput<TH1F*>("h_Eff_pt_TightID_ISO03PFWeighted_Mu2");
  h_Eff_pt_TightID_ISO04PFWeighted[0]   = FindOutput<TH1F*>("h_Eff_pt_TightID_ISO04PFWeighted_Mu1");
  h_Eff_pt_TightID_ISO04PFWeighted[1]   = FindOutput<TH1F*>("h_Eff_pt_TightID_ISO04PFWeighted_Mu2"); 
  h_Eff_pt_TightID_ISO03PUPPI[0]        = FindOutput<TH1F*>("h_Eff_pt_TightID_ISO03PUPPI_Mu1");
  h_Eff_pt_TightID_ISO03PUPPI[1]        = FindOutput<TH1F*>("h_Eff_pt_TightID_ISO03PUPPI_Mu2");
  h_Eff_pt_TightID_ISO04PUPPI[0]        = FindOutput<TH1F*>("h_Eff_pt_TightID_ISO04PUPPI_Mu1");
  h_Eff_pt_TightID_ISO04PUPPI[1]        = FindOutput<TH1F*>("h_Eff_pt_TightID_ISO04PUPPI_Mu2"); 
  h_Eff_pt_MediumID_ISO03[0]            = FindOutput<TH1F*>("h_Eff_pt_MediumID_ISO03_Mu1");
  h_Eff_pt_MediumID_ISO03[1]            = FindOutput<TH1F*>("h_Eff_pt_MediumID_ISO03_Mu2");
  h_Eff_pt_MediumID_ISO04[0]            = FindOutput<TH1F*>("h_Eff_pt_MediumID_ISO04_Mu1");
  h_Eff_pt_MediumID_ISO04[1]            = FindOutput<TH1F*>("h_Eff_pt_MediumID_ISO04_Mu2");
  h_Eff_pt_MediumID_ISO03dBeta[0]       = FindOutput<TH1F*>("h_Eff_pt_MediumID_ISO03dBeta_Mu1");
  h_Eff_pt_MediumID_ISO03dBeta[1]       = FindOutput<TH1F*>("h_Eff_pt_MediumID_ISO03dBeta_Mu2");
  h_Eff_pt_MediumID_ISO04dBeta[0]       = FindOutput<TH1F*>("h_Eff_pt_MediumID_ISO04dBeta_Mu1");
  h_Eff_pt_MediumID_ISO04dBeta[1]       = FindOutput<TH1F*>("h_Eff_pt_MediumID_ISO04dBeta_Mu2");
  h_Eff_pt_MediumID_ISO03PFWeighted[0]  = FindOutput<TH1F*>("h_Eff_pt_MediumID_ISO03PFWeighted_Mu1");
  h_Eff_pt_MediumID_ISO03PFWeighted[1]  = FindOutput<TH1F*>("h_Eff_pt_MediumID_ISO03PFWeighted_Mu2");
  h_Eff_pt_MediumID_ISO04PFWeighted[0]  = FindOutput<TH1F*>("h_Eff_pt_MediumID_ISO04PFWeighted_Mu1");
  h_Eff_pt_MediumID_ISO04PFWeighted[1]  = FindOutput<TH1F*>("h_Eff_pt_MediumID_ISO04PFWeighted_Mu2"); 
  h_Eff_pt_MediumID_ISO03PUPPI[0]       = FindOutput<TH1F*>("h_Eff_pt_MediumID_ISO03PUPPI_Mu1");
  h_Eff_pt_MediumID_ISO03PUPPI[1]       = FindOutput<TH1F*>("h_Eff_pt_MediumID_ISO03PUPPI_Mu2");
  h_Eff_pt_MediumID_ISO04PUPPI[0]       = FindOutput<TH1F*>("h_Eff_pt_MediumID_ISO04PUPPI_Mu1");
  h_Eff_pt_MediumID_ISO04PUPPI[1]       = FindOutput<TH1F*>("h_Eff_pt_MediumID_ISO04PUPPI_Mu2"); 

  h_Eff_eta_Matched[0]                  = FindOutput<TH1F*>("h_Eff_eta_Matched_Mu1");
  h_Eff_eta_Matched[1]                  = FindOutput<TH1F*>("h_Eff_eta_Matched_Mu2");
  h_Eff_eta_TightID[0]                  = FindOutput<TH1F*>("h_Eff_eta_TightID_Mu1");
  h_Eff_eta_TightID[1]                  = FindOutput<TH1F*>("h_Eff_eta_TightID_Mu2");
  h_Eff_eta_MediumID[0]                 = FindOutput<TH1F*>("h_Eff_eta_MediumID_Mu1");
  h_Eff_eta_MediumID[1]                 = FindOutput<TH1F*>("h_Eff_eta_MediumID_Mu2");
  h_Eff_eta_HWWID[0]                    = FindOutput<TH1F*>("h_Eff_eta_HWWID_Mu1");
  h_Eff_eta_HWWID[1]                    = FindOutput<TH1F*>("h_Eff_eta_HWWID_Mu2");
  h_Eff_eta_TightIDGoT[0]               = FindOutput<TH1F*>("h_Eff_eta_TightIDGoT_Mu1");
  h_Eff_eta_TightIDGoT[1]               = FindOutput<TH1F*>("h_Eff_eta_TightIDGoT_Mu2");
  h_Eff_eta_TightIDipsHWW[0]            = FindOutput<TH1F*>("h_Eff_eta_TightIDipsHWW_Mu1");
  h_Eff_eta_TightIDipsHWW[1]            = FindOutput<TH1F*>("h_Eff_eta_TightIDipsHWW_Mu2");
  h_Eff_eta_MediumIDipsHWW[0]           = FindOutput<TH1F*>("h_Eff_eta_MediumIDipsHWW_Mu1");
  h_Eff_eta_MediumIDipsHWW[1]           = FindOutput<TH1F*>("h_Eff_eta_MediumIDipsHWW_Mu2");
  h_Eff_eta_TightID_ISO03[0]            = FindOutput<TH1F*>("h_Eff_eta_TightID_ISO03_Mu1");
  h_Eff_eta_TightID_ISO03[1]            = FindOutput<TH1F*>("h_Eff_eta_TightID_ISO03_Mu2");
  h_Eff_eta_TightID_ISO04[0]            = FindOutput<TH1F*>("h_Eff_eta_TightID_ISO04_Mu1");
  h_Eff_eta_TightID_ISO04[1]            = FindOutput<TH1F*>("h_Eff_eta_TightID_ISO04_Mu2");
  h_Eff_eta_TightID_ISO03dBeta[0]       = FindOutput<TH1F*>("h_Eff_eta_TightID_ISO03dBeta_Mu1");
  h_Eff_eta_TightID_ISO03dBeta[1]       = FindOutput<TH1F*>("h_Eff_eta_TightID_ISO03dBeta_Mu2");
  h_Eff_eta_TightID_ISO04dBeta[0]       = FindOutput<TH1F*>("h_Eff_eta_TightID_ISO04dBeta_Mu1");
  h_Eff_eta_TightID_ISO04dBeta[1]       = FindOutput<TH1F*>("h_Eff_eta_TightID_ISO04dBeta_Mu2");
  h_Eff_eta_TightID_ISO03PFWeighted[0]  = FindOutput<TH1F*>("h_Eff_eta_TightID_ISO03PFWeighted_Mu1");
  h_Eff_eta_TightID_ISO03PFWeighted[1]  = FindOutput<TH1F*>("h_Eff_eta_TightID_ISO03PFWeighted_Mu2");
  h_Eff_eta_TightID_ISO04PFWeighted[0]  = FindOutput<TH1F*>("h_Eff_eta_TightID_ISO04PFWeighted_Mu1");
  h_Eff_eta_TightID_ISO04PFWeighted[1]  = FindOutput<TH1F*>("h_Eff_eta_TightID_ISO04PFWeighted_Mu2"); 
  h_Eff_eta_TightID_ISO03PUPPI[0]       = FindOutput<TH1F*>("h_Eff_eta_TightID_ISO03PUPPI_Mu1");
  h_Eff_eta_TightID_ISO03PUPPI[1]       = FindOutput<TH1F*>("h_Eff_eta_TightID_ISO03PUPPI_Mu2");
  h_Eff_eta_TightID_ISO04PUPPI[0]       = FindOutput<TH1F*>("h_Eff_eta_TightID_ISO04PUPPI_Mu1");
  h_Eff_eta_TightID_ISO04PUPPI[1]       = FindOutput<TH1F*>("h_Eff_eta_TightID_ISO04PUPPI_Mu2");
  h_Eff_eta_MediumID_ISO03[0]           = FindOutput<TH1F*>("h_Eff_eta_MediumID_ISO03_Mu1");
  h_Eff_eta_MediumID_ISO03[1]           = FindOutput<TH1F*>("h_Eff_eta_MediumID_ISO03_Mu2");
  h_Eff_eta_MediumID_ISO04[0]           = FindOutput<TH1F*>("h_Eff_eta_MediumID_ISO04_Mu1");
  h_Eff_eta_MediumID_ISO04[1]           = FindOutput<TH1F*>("h_Eff_eta_MediumID_ISO04_Mu2");
  h_Eff_eta_MediumID_ISO03dBeta[0]      = FindOutput<TH1F*>("h_Eff_eta_MediumID_ISO03dBeta_Mu1");
  h_Eff_eta_MediumID_ISO03dBeta[1]      = FindOutput<TH1F*>("h_Eff_eta_MediumID_ISO03dBeta_Mu2");
  h_Eff_eta_MediumID_ISO04dBeta[0]      = FindOutput<TH1F*>("h_Eff_eta_MediumID_ISO04dBeta_Mu1");
  h_Eff_eta_MediumID_ISO04dBeta[1]      = FindOutput<TH1F*>("h_Eff_eta_MediumID_ISO04dBeta_Mu2");
  h_Eff_eta_MediumID_ISO03PFWeighted[0] = FindOutput<TH1F*>("h_Eff_eta_MediumID_ISO03PFWeighted_Mu1");
  h_Eff_eta_MediumID_ISO03PFWeighted[1] = FindOutput<TH1F*>("h_Eff_eta_MediumID_ISO03PFWeighted_Mu2");
  h_Eff_eta_MediumID_ISO04PFWeighted[0] = FindOutput<TH1F*>("h_Eff_eta_MediumID_ISO04PFWeighted_Mu1");
  h_Eff_eta_MediumID_ISO04PFWeighted[1] = FindOutput<TH1F*>("h_Eff_eta_MediumID_ISO04PFWeighted_Mu2"); 
  h_Eff_eta_MediumID_ISO03PUPPI[0]      = FindOutput<TH1F*>("h_Eff_eta_MediumID_ISO03PUPPI_Mu1");
  h_Eff_eta_MediumID_ISO03PUPPI[1]      = FindOutput<TH1F*>("h_Eff_eta_MediumID_ISO03PUPPI_Mu2");
  h_Eff_eta_MediumID_ISO04PUPPI[0]      = FindOutput<TH1F*>("h_Eff_eta_MediumID_ISO04PUPPI_Mu1");
  h_Eff_eta_MediumID_ISO04PUPPI[1]      = FindOutput<TH1F*>("h_Eff_eta_MediumID_ISO04PUPPI_Mu2");

  h_Eff_npv_Matched[0]                  = FindOutput<TH1F*>("h_Eff_npv_Matched_Mu1");
  h_Eff_npv_Matched[1]                  = FindOutput<TH1F*>("h_Eff_npv_Matched_Mu2");
  h_Eff_npv_TightID[0]                  = FindOutput<TH1F*>("h_Eff_npv_TightID_Mu1");
  h_Eff_npv_TightID[1]                  = FindOutput<TH1F*>("h_Eff_npv_TightID_Mu2");
  h_Eff_npv_MediumID[0]                 = FindOutput<TH1F*>("h_Eff_npv_MediumID_Mu1");
  h_Eff_npv_MediumID[1]                 = FindOutput<TH1F*>("h_Eff_npv_MediumID_Mu2");
  h_Eff_npv_HWWID[0]                    = FindOutput<TH1F*>("h_Eff_npv_HWWID_Mu1");
  h_Eff_npv_HWWID[1]                    = FindOutput<TH1F*>("h_Eff_npv_HWWID_Mu2");
  h_Eff_npv_TightIDGoT[0]               = FindOutput<TH1F*>("h_Eff_npv_TightIDGoT_Mu1");
  h_Eff_npv_TightIDGoT[1]               = FindOutput<TH1F*>("h_Eff_npv_TightIDGoT_Mu2");
  h_Eff_npv_TightIDipsHWW[0]            = FindOutput<TH1F*>("h_Eff_npv_TightIDipsHWW_Mu1");
  h_Eff_npv_TightIDipsHWW[1]            = FindOutput<TH1F*>("h_Eff_npv_TightIDipsHWW_Mu2");
  h_Eff_npv_MediumIDipsHWW[0]           = FindOutput<TH1F*>("h_Eff_npv_MediumIDipsHWW_Mu1");
  h_Eff_npv_MediumIDipsHWW[1]           = FindOutput<TH1F*>("h_Eff_npv_MediumIDipsHWW_Mu2");
  h_Eff_npv_TightID_ISO03[0]            = FindOutput<TH1F*>("h_Eff_npv_TightID_ISO03_Mu1");
  h_Eff_npv_TightID_ISO03[1]            = FindOutput<TH1F*>("h_Eff_npv_TightID_ISO03_Mu2");
  h_Eff_npv_TightID_ISO04[0]            = FindOutput<TH1F*>("h_Eff_npv_TightID_ISO04_Mu1");
  h_Eff_npv_TightID_ISO04[1]            = FindOutput<TH1F*>("h_Eff_npv_TightID_ISO04_Mu2");
  h_Eff_npv_TightID_ISO03dBeta[0]       = FindOutput<TH1F*>("h_Eff_npv_TightID_ISO03dBeta_Mu1");
  h_Eff_npv_TightID_ISO03dBeta[1]       = FindOutput<TH1F*>("h_Eff_npv_TightID_ISO03dBeta_Mu2");
  h_Eff_npv_TightID_ISO04dBeta[0]       = FindOutput<TH1F*>("h_Eff_npv_TightID_ISO04dBeta_Mu1");
  h_Eff_npv_TightID_ISO04dBeta[1]       = FindOutput<TH1F*>("h_Eff_npv_TightID_ISO04dBeta_Mu2");
  h_Eff_npv_TightID_ISO03PFWeighted[0]  = FindOutput<TH1F*>("h_Eff_npv_TightID_ISO03PFWeighted_Mu1");
  h_Eff_npv_TightID_ISO03PFWeighted[1]  = FindOutput<TH1F*>("h_Eff_npv_TightID_ISO03PFWeighted_Mu2");
  h_Eff_npv_TightID_ISO04PFWeighted[0]  = FindOutput<TH1F*>("h_Eff_npv_TightID_ISO04PFWeighted_Mu1");
  h_Eff_npv_TightID_ISO04PFWeighted[1]  = FindOutput<TH1F*>("h_Eff_npv_TightID_ISO04PFWeighted_Mu2"); 
  h_Eff_npv_TightID_ISO03PUPPI[0]       = FindOutput<TH1F*>("h_Eff_npv_TightID_ISO03PUPPI_Mu1");
  h_Eff_npv_TightID_ISO03PUPPI[1]       = FindOutput<TH1F*>("h_Eff_npv_TightID_ISO03PUPPI_Mu2");
  h_Eff_npv_TightID_ISO04PUPPI[0]       = FindOutput<TH1F*>("h_Eff_npv_TightID_ISO04PUPPI_Mu1");
  h_Eff_npv_TightID_ISO04PUPPI[1]       = FindOutput<TH1F*>("h_Eff_npv_TightID_ISO04PUPPI_Mu2");
  h_Eff_npv_MediumID_ISO03[0]           = FindOutput<TH1F*>("h_Eff_npv_MediumID_ISO03_Mu1");
  h_Eff_npv_MediumID_ISO03[1]           = FindOutput<TH1F*>("h_Eff_npv_MediumID_ISO03_Mu2");
  h_Eff_npv_MediumID_ISO04[0]           = FindOutput<TH1F*>("h_Eff_npv_MediumID_ISO04_Mu1");
  h_Eff_npv_MediumID_ISO04[1]           = FindOutput<TH1F*>("h_Eff_npv_MediumID_ISO04_Mu2");
  h_Eff_npv_MediumID_ISO03dBeta[0]      = FindOutput<TH1F*>("h_Eff_npv_MediumID_ISO03dBeta_Mu1");
  h_Eff_npv_MediumID_ISO03dBeta[1]      = FindOutput<TH1F*>("h_Eff_npv_MediumID_ISO03dBeta_Mu2");
  h_Eff_npv_MediumID_ISO04dBeta[0]      = FindOutput<TH1F*>("h_Eff_npv_MediumID_ISO04dBeta_Mu1");
  h_Eff_npv_MediumID_ISO04dBeta[1]      = FindOutput<TH1F*>("h_Eff_npv_MediumID_ISO04dBeta_Mu2");
  h_Eff_npv_MediumID_ISO03PFWeighted[0] = FindOutput<TH1F*>("h_Eff_npv_MediumID_ISO03PFWeighted_Mu1");
  h_Eff_npv_MediumID_ISO03PFWeighted[1] = FindOutput<TH1F*>("h_Eff_npv_MediumID_ISO03PFWeighted_Mu2");
  h_Eff_npv_MediumID_ISO04PFWeighted[0] = FindOutput<TH1F*>("h_Eff_npv_MediumID_ISO04PFWeighted_Mu1");
  h_Eff_npv_MediumID_ISO04PFWeighted[1] = FindOutput<TH1F*>("h_Eff_npv_MediumID_ISO04PFWeighted_Mu2"); 
  h_Eff_npv_MediumID_ISO03PUPPI[0]      = FindOutput<TH1F*>("h_Eff_npv_MediumID_ISO03PUPPI_Mu1");
  h_Eff_npv_MediumID_ISO03PUPPI[1]      = FindOutput<TH1F*>("h_Eff_npv_MediumID_ISO03PUPPI_Mu2");
  h_Eff_npv_MediumID_ISO04PUPPI[0]      = FindOutput<TH1F*>("h_Eff_npv_MediumID_ISO04PUPPI_Mu1");
  h_Eff_npv_MediumID_ISO04PUPPI[1]      = FindOutput<TH1F*>("h_Eff_npv_MediumID_ISO04PUPPI_Mu2");

  // Dimuon ISO efficiencies vs pt, eta and npv
  h_Eff_pt_TightID_Dilep                   = FindOutput<TH1F*>("h_Eff_pt_TightID_Dilep");
  h_Eff_pt_TightID_ISO03_Dilep             = FindOutput<TH1F*>("h_Eff_pt_TightID_ISO03_Dilep");
  h_Eff_pt_TightID_ISO04_Dilep             = FindOutput<TH1F*>("h_Eff_pt_TightID_ISO04_Dilep");
  h_Eff_pt_TightID_ISO03dBeta_Dilep        = FindOutput<TH1F*>("h_Eff_pt_TightID_ISO03dBeta_Dilep");
  h_Eff_pt_TightID_ISO04dBeta_Dilep        = FindOutput<TH1F*>("h_Eff_pt_TightID_ISO04dBeta_Dilep");
  h_Eff_pt_TightID_ISO03PFWeighted_Dilep   = FindOutput<TH1F*>("h_Eff_pt_TightID_ISO03PFWeighted_Dilep");
  h_Eff_pt_TightID_ISO04PFWeighted_Dilep   = FindOutput<TH1F*>("h_Eff_pt_TightID_ISO04PFWeighted_Dilep");
  h_Eff_pt_TightID_ISO03PUPPI_Dilep        = FindOutput<TH1F*>("h_Eff_pt_TightID_ISO03PUPPI_Dilep");
  h_Eff_pt_TightID_ISO04PUPPI_Dilep        = FindOutput<TH1F*>("h_Eff_pt_TightID_ISO04PUPPI_Dilep");
  h_Eff_pt_MediumID_Dilep                  = FindOutput<TH1F*>("h_Eff_pt_MediumID_Dilep");
  h_Eff_pt_MediumID_ISO03_Dilep            = FindOutput<TH1F*>("h_Eff_pt_MediumID_ISO03_Dilep");
  h_Eff_pt_MediumID_ISO04_Dilep            = FindOutput<TH1F*>("h_Eff_pt_MediumID_ISO04_Dilep");
  h_Eff_pt_MediumID_ISO03dBeta_Dilep       = FindOutput<TH1F*>("h_Eff_pt_MediumID_ISO03dBeta_Dilep");
  h_Eff_pt_MediumID_ISO04dBeta_Dilep       = FindOutput<TH1F*>("h_Eff_pt_MediumID_ISO04dBeta_Dilep");
  h_Eff_pt_MediumID_ISO03PFWeighted_Dilep  = FindOutput<TH1F*>("h_Eff_pt_MediumID_ISO03PFWeighted_Dilep");
  h_Eff_pt_MediumID_ISO04PFWeighted_Dilep  = FindOutput<TH1F*>("h_Eff_pt_MediumID_ISO04PFWeighted_Dilep");
  h_Eff_pt_MediumID_ISO03PUPPI_Dilep       = FindOutput<TH1F*>("h_Eff_pt_MediumID_ISO03PUPPI_Dilep");
  h_Eff_pt_MediumID_ISO04PUPPI_Dilep       = FindOutput<TH1F*>("h_Eff_pt_MediumID_ISO04PUPPI_Dilep");

  h_Eff_eta_TightID_Dilep                  = FindOutput<TH1F*>("h_Eff_eta_TightID_Dilep");
  h_Eff_eta_TightID_ISO03_Dilep            = FindOutput<TH1F*>("h_Eff_eta_TightID_ISO03_Dilep");
  h_Eff_eta_TightID_ISO04_Dilep            = FindOutput<TH1F*>("h_Eff_eta_TightID_ISO04_Dilep");
  h_Eff_eta_TightID_ISO03dBeta_Dilep       = FindOutput<TH1F*>("h_Eff_eta_TightID_ISO03dBeta_Dilep");
  h_Eff_eta_TightID_ISO04dBeta_Dilep       = FindOutput<TH1F*>("h_Eff_eta_TightID_ISO04dBeta_Dilep");
  h_Eff_eta_TightID_ISO03PFWeighted_Dilep  = FindOutput<TH1F*>("h_Eff_eta_TightID_ISO03PFWeighted_Dilep");
  h_Eff_eta_TightID_ISO04PFWeighted_Dilep  = FindOutput<TH1F*>("h_Eff_eta_TightID_ISO04PFWeighted_Dilep");
  h_Eff_eta_TightID_ISO03PUPPI_Dilep       = FindOutput<TH1F*>("h_Eff_eta_TightID_ISO03PUPPI_Dilep");
  h_Eff_eta_TightID_ISO04PUPPI_Dilep       = FindOutput<TH1F*>("h_Eff_eta_TightID_ISO04PUPPI_Dilep");
  h_Eff_eta_MediumID_Dilep                 = FindOutput<TH1F*>("h_Eff_eta_MediumID_Dilep");
  h_Eff_eta_MediumID_ISO03_Dilep           = FindOutput<TH1F*>("h_Eff_eta_MediumID_ISO03_Dilep");
  h_Eff_eta_MediumID_ISO04_Dilep           = FindOutput<TH1F*>("h_Eff_eta_MediumID_ISO04_Dilep");
  h_Eff_eta_MediumID_ISO03dBeta_Dilep      = FindOutput<TH1F*>("h_Eff_eta_MediumID_ISO03dBeta_Dilep");
  h_Eff_eta_MediumID_ISO04dBeta_Dilep      = FindOutput<TH1F*>("h_Eff_eta_MediumID_ISO04dBeta_Dilep");
  h_Eff_eta_MediumID_ISO03PFWeighted_Dilep = FindOutput<TH1F*>("h_Eff_eta_MediumID_ISO03PFWeighted_Dilep");
  h_Eff_eta_MediumID_ISO04PFWeighted_Dilep = FindOutput<TH1F*>("h_Eff_eta_MediumID_ISO04PFWeighted_Dilep");
  h_Eff_eta_MediumID_ISO03PUPPI_Dilep      = FindOutput<TH1F*>("h_Eff_eta_MediumID_ISO03PUPPI_Dilep");
  h_Eff_eta_MediumID_ISO04PUPPI_Dilep      = FindOutput<TH1F*>("h_Eff_eta_MediumID_ISO04PUPPI_Dilep");

  h_Eff_npv_TightID_Dilep                  = FindOutput<TH1F*>("h_Eff_npv_TightID_Dilep");
  h_Eff_npv_TightID_ISO03_Dilep            = FindOutput<TH1F*>("h_Eff_npv_TightID_ISO03_Dilep");
  h_Eff_npv_TightID_ISO04_Dilep            = FindOutput<TH1F*>("h_Eff_npv_TightID_ISO04_Dilep");
  h_Eff_npv_TightID_ISO03dBeta_Dilep       = FindOutput<TH1F*>("h_Eff_npv_TightID_ISO03dBeta_Dilep");
  h_Eff_npv_TightID_ISO04dBeta_Dilep       = FindOutput<TH1F*>("h_Eff_npv_TightID_ISO04dBeta_Dilep");
  h_Eff_npv_TightID_ISO03PFWeighted_Dilep  = FindOutput<TH1F*>("h_Eff_npv_TightID_ISO03PFWeighted_Dilep");
  h_Eff_npv_TightID_ISO04PFWeighted_Dilep  = FindOutput<TH1F*>("h_Eff_npv_TightID_ISO04PFWeighted_Dilep");
  h_Eff_npv_TightID_ISO03PUPPI_Dilep       = FindOutput<TH1F*>("h_Eff_npv_TightID_ISO03PUPPI_Dilep");
  h_Eff_npv_TightID_ISO04PUPPI_Dilep       = FindOutput<TH1F*>("h_Eff_npv_TightID_ISO04PUPPI_Dilep");
  h_Eff_npv_MediumID_Dilep                 = FindOutput<TH1F*>("h_Eff_npv_MediumID_Dilep");
  h_Eff_npv_MediumID_ISO03_Dilep           = FindOutput<TH1F*>("h_Eff_npv_MediumID_ISO03_Dilep");
  h_Eff_npv_MediumID_ISO04_Dilep           = FindOutput<TH1F*>("h_Eff_npv_MediumID_ISO04_Dilep");
  h_Eff_npv_MediumID_ISO03dBeta_Dilep      = FindOutput<TH1F*>("h_Eff_npv_MediumID_ISO03dBeta_Dilep");
  h_Eff_npv_MediumID_ISO04dBeta_Dilep      = FindOutput<TH1F*>("h_Eff_npv_MediumID_ISO04dBeta_Dilep");
  h_Eff_npv_MediumID_ISO03PFWeighted_Dilep = FindOutput<TH1F*>("h_Eff_npv_MediumID_ISO03PFWeighted_Dilep");
  h_Eff_npv_MediumID_ISO04PFWeighted_Dilep = FindOutput<TH1F*>("h_Eff_npv_MediumID_ISO04PFWeighted_Dilep");
  h_Eff_npv_MediumID_ISO03PUPPI_Dilep      = FindOutput<TH1F*>("h_Eff_npv_MediumID_ISO03PUPPI_Dilep");
  h_Eff_npv_MediumID_ISO04PUPPI_Dilep      = FindOutput<TH1F*>("h_Eff_npv_MediumID_ISO04PUPPI_Dilep");

  //ISO ROC Curves
  h_RC_TightID_ISO03_Dilep            = FindOutput<TH1F*>("h_RC_TightID_ISO03_Dilep");
  h_RC_TightID_ISO04_Dilep            = FindOutput<TH1F*>("h_RC_TightID_ISO04_Dilep");
  h_RC_TightID_ISO03dBeta_Dilep       = FindOutput<TH1F*>("h_RC_TightID_ISO03dBeta_Dilep");
  h_RC_TightID_ISO04dBeta_Dilep       = FindOutput<TH1F*>("h_RC_TightID_ISO04dBeta_Dilep");
  h_RC_TightID_ISO03PFWeighted_Dilep  = FindOutput<TH1F*>("h_RC_TightID_ISO03PFWeighted_Dilep");
  h_RC_TightID_ISO04PFWeighted_Dilep  = FindOutput<TH1F*>("h_RC_TightID_ISO04PFWeighted_Dilep");
  h_RC_TightID_ISO03PUPPI_Dilep       = FindOutput<TH1F*>("h_RC_TightID_ISO03PUPPI_Dilep");
  h_RC_TightID_ISO04PUPPI_Dilep       = FindOutput<TH1F*>("h_RC_TightID_ISO04PUPPI_Dilep");
  h_RC_MediumID_ISO03_Dilep           = FindOutput<TH1F*>("h_RC_MediumID_ISO03_Dilep");
  h_RC_MediumID_ISO04_Dilep           = FindOutput<TH1F*>("h_RC_MediumID_ISO04_Dilep");
  h_RC_MediumID_ISO03dBeta_Dilep      = FindOutput<TH1F*>("h_RC_MediumID_ISO03dBeta_Dilep");
  h_RC_MediumID_ISO04dBeta_Dilep      = FindOutput<TH1F*>("h_RC_MediumID_ISO04dBeta_Dilep");
  h_RC_MediumID_ISO03PFWeighted_Dilep = FindOutput<TH1F*>("h_RC_MediumID_ISO03PFWeighted_Dilep");
  h_RC_MediumID_ISO04PFWeighted_Dilep = FindOutput<TH1F*>("h_RC_MediumID_ISO04PFWeighted_Dilep");
  h_RC_MediumID_ISO03PUPPI_Dilep      = FindOutput<TH1F*>("h_RC_MediumID_ISO03PUPPI_Dilep");
  h_RC_MediumID_ISO04PUPPI_Dilep      = FindOutput<TH1F*>("h_RC_MediumID_ISO04PUPPI_Dilep");

  // Final Report

  cout << " ---------------------------------------------------" << endl;
  cout << " " << endl;  
  cout << " Number of Events in the sample:  " << _NEvents  << endl;
  cout << " Normalization factor: " << _factN << endl;
  cout << endl;

  if (_Report) {
    cout << " ---------------------------------------------------" << endl;
    cout << "Counting Report" << endl;
    cout << " ---------------------------------------------------" << endl;
    cout << "" << endl;
    cout << "Total NoE (Number of Events): " << GCount_AllEvents << endl;
    cout << " NoE with 2 GEN muons: " << GCount_GenEvents << endl;
    cout << " * NoE with at least 2 valid muons: " << GCount_Fiducial_AtLeast2 << 
      " (" << 100*GCount_Fiducial_AtLeast2/GCount_GenEvents << "%)" << endl;
    cout << "  * NoE with exactly 2 valid muons: " << GCount_Fiducial_2 << 
      " (" << 100*GCount_Fiducial_2/GCount_GenEvents << "%)" << endl;
    
    cout << "   + NoE where the valid muons are the 1st and the 2nd: " << GCount_Fiducial_1st2nd <<
      " (" << 100*GCount_Fiducial_1st2nd/GCount_GenEvents << "%)" << endl;
    
    cout << "     - NoE where both are matched: " << GCount_Match_1st2nd <<
      " (" << 100*GCount_Match_1st2nd/GCount_Fiducial_1st2nd << "%)" << endl;
    cout << "      ~ NoE where both pass Tight ID: " << GCount_MatchTight_1st2nd <<
      " (" << 100*GCount_MatchTight_1st2nd/GCount_Match_1st2nd << "%)" << endl;
    cout << "      ~ NoE where both pass Tight ID and ISO: " << GCount_MatchTightIso_1st2nd <<
      " (" << 100*GCount_MatchTightIso_1st2nd/GCount_Match_1st2nd << "%)" << endl;
    
    cout << "     - NoE where both are not matched: " << GCount_NoMatch_1st2nd <<
      " (" << 100*GCount_NoMatch_1st2nd/GCount_Fiducial_1st2nd << "%)" << endl;
    cout << "      ~ NoE where both pass Tight ID: " << GCount_Tight_1st2nd <<
      " (" << 100*GCount_Tight_1st2nd/GCount_NoMatch_1st2nd << "%)" << endl;
    cout << "      ~ NoE where both pass Tight ID and ISO: " << GCount_TightIso_1st2nd <<
      " (" << 100*GCount_TightIso_1st2nd/GCount_NoMatch_1st2nd << "%)" << endl;
    
    cout << "   + NoE where the valid muons are the 1st and the 3rd: " << GCount_Fiducial_1st3rd << 
      " (" << 100*GCount_Fiducial_1st3rd/GCount_GenEvents << "%)" << endl;
    
    cout << "     - NoE where both are matched: " << GCount_Match_1st3rd <<
      " (" << 100*GCount_Match_1st3rd/GCount_Fiducial_1st3rd << "%)" << endl;
    cout << "      ~ NoE where both pass Tight ID: " << GCount_MatchTight_1st3rd <<
      " (" << 100*GCount_MatchTight_1st3rd/GCount_Match_1st3rd << "%)" << endl;
    cout << "      ~ NoE where both pass Tight ID and ISO: " << GCount_MatchTightIso_1st3rd <<
      " (" << 100*GCount_MatchTightIso_1st3rd/GCount_Match_1st3rd << "%)" << endl;
    
    cout << "     - NoE where both are not matched: " << GCount_NoMatch_1st3rd <<
      " (" << 100*GCount_NoMatch_1st3rd/GCount_Fiducial_1st3rd << "%)" << endl;
    cout << "      ~ NoE where both pass Tight ID: " << GCount_Tight_1st3rd <<
      " (" << 100*GCount_Tight_1st3rd/GCount_NoMatch_1st3rd << "%)" << endl;
    cout << "      ~ NoE where both pass Tight ID and ISO: " << GCount_TightIso_1st3rd <<
      " (" << 100*GCount_TightIso_1st3rd/GCount_NoMatch_1st3rd << "%)" << endl;
    
    cout << "   + NoE where the valid muons are the 1st and other one: " << GCount_Fiducial_1stOther <<
      " (" << 100*GCount_Fiducial_1stOther/GCount_GenEvents << "%)" << endl; 
    
    cout << "  * NoE with more than 2 valid muons: " << GCount_Fiducial_MoreThan2 <<
      " (" << 100*GCount_Fiducial_MoreThan2/GCount_GenEvents << "%)" << endl;
    
    cout << "   + NoE with 2 muons matched: " << GCount_Match_MoreThan2_OK << 
      " (" << 100*GCount_Match_MoreThan2_OK/GCount_GenEvents << "%)" << endl;
    cout << "   + NoE with more than 2 muons matched: " << GCount_Match_MoreThan2 << 
      " (" << 100*GCount_Match_MoreThan2/GCount_GenEvents << "%)" << endl;
    cout << "   + NoE with less than 2 muons matched: " << GCount_NoMatch_MoreThan2 << 
      " (" << 100*GCount_NoMatch_MoreThan2/GCount_GenEvents << "%)" << endl; 
    cout << "   + NoE with exactly 3 valid muons: " << GCount_Fiducial_3 << 
      " (" << 100*GCount_Fiducial_3/GCount_GenEvents << "%)" << endl;
    cout << "   + NoE with more than 3 valid muons: " << GCount_Fiducial_MoreThan3 << 
      " (" << 100*GCount_Fiducial_MoreThan3/GCount_GenEvents << "%)" << endl;
    cout << " * NoE where only the 1st muon is valid: " << GCount_Fiducial_Only1st << 
      
      " (" << 100*GCount_Fiducial_Only1st/GCount_GenEvents << "%)" << endl;
    cout << " * NoE without valid muons: " << GCount_Fiducial_None << 
      
      " (" << 100*GCount_Fiducial_None/GCount_GenEvents << "%)" << endl;
    cout << "" << endl;
    cout << " ---------------------------------------------------" << endl;
  }

}
