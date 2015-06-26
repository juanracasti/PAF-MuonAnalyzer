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

#include "MuonEfficienciesSelector.h"

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

ClassImp(MuonEfficienciesSelector)

// Initialise input parameters and data members for all events
void MuonEfficienciesSelector::Initialise() {

  _Signal     = GetParam<TString>("Signal");
  _IsDATA     = GetParam<bool>("IsDATA");
  _NEvents    = GetParam<int>("NEvents");
  _WhichRun   = GetParam<int>("WhichRun"); 
  _Debug      = GetParam<bool>("Debug");

 
//------------------------------------------------------------------------------
// Create histos
//------------------------------------------------------------------------------

  // Single Muon ID and ISO efficiencies vs pt, eta and npv

  h_Eff_pt_NoID[0]                      = CreateH1F("h_Eff_pt_NoID_Mu1",
						    "h_Eff_pt_NoID_Mu1",                       200, 0, 200);
  h_Eff_pt_NoID[0]->TH1::SetDefaultSumw2();
  h_Eff_pt_NoID[1]                      = CreateH1F("h_Eff_pt_NoID_Mu2",
						    "h_Eff_pt_NoID_Mu2",                       200, 0, 200);
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
  
  h_Eff_eta_NoID[0]                     = CreateH1F("h_Eff_eta_NoID_Mu1",
						    "h_Eff_eta_NoID_Mu1",                      50, -2.5, 2.5);
  h_Eff_eta_NoID[1]                     = CreateH1F("h_Eff_eta_NoID_Mu2",
						    "h_Eff_eta_NoID_Mu2",                      50, -2.5, 2.5);
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
  
  h_Eff_npv_NoID[0]                     = CreateH1F("h_Eff_npv_NoID_Mu1",
						    "h_Eff_npv_NoID_Mu1",                      45, 0, 45);
  h_Eff_npv_NoID[1]                     = CreateH1F("h_Eff_npv_NoID_Mu2",
						    "h_Eff_npv_NoID_Mu2",                      45, 0, 45);
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

  // ISO efficiencies vs pt, eta and npv for all muons

  h_Eff_pt_NoID_AllMu                      = CreateH1F("h_Eff_pt_NoID_AllMu",
						       "h_Eff_pt_NoID_AllMu",                       200, 0, 200);
  h_Eff_pt_HWWID_AllMu                     = CreateH1F("h_Eff_pt_HWWID_AllMu",
						       "h_Eff_pt_HWWID_AllMu",                      200, 0, 200);
  h_Eff_pt_TightIDGoT_AllMu                = CreateH1F("h_Eff_pt_TightIDGoT_AllMu",
						       "h_Eff_pt_TightIDGoT_AllMu",                 200, 0, 200);
  h_Eff_pt_TightIDipsHWW_AllMu             = CreateH1F("h_Eff_pt_TightIDipsHWW_AllMu",
						       "h_Eff_pt_TightIDipsHWW_AllMu",              200, 0, 200);
  h_Eff_pt_MediumIDipsHWW_AllMu            = CreateH1F("h_Eff_pt_MediumIDipsHWW_AllMu",
						       "h_Eff_pt_MediumIDipsHWW_AllMu",             200, 0, 200);
  h_Eff_pt_TightID_AllMu                   = CreateH1F("h_Eff_pt_TightID_AllMu", 
						       "h_Eff_pt_TightID_AllMu",                    200, 0, 200);
  h_Eff_pt_TightID_ISO03_AllMu             = CreateH1F("h_Eff_pt_TightID_ISO03_AllMu", 
						       "h_Eff_pt_TightID_ISO03_AllMu",              200, 0, 200);
  h_Eff_pt_TightID_ISO04_AllMu             = CreateH1F("h_Eff_pt_TightID_ISO04_AllMu", 
						       "h_Eff_pt_TightID_ISO04_AllMu",              200, 0, 200);
  h_Eff_pt_TightID_ISO03dBeta_AllMu        = CreateH1F("h_Eff_pt_TightID_ISO03dBeta_AllMu", 
						       "h_Eff_pt_TightID_ISO03dBeta_AllMu",         200, 0, 200);
  h_Eff_pt_TightID_ISO04dBeta_AllMu        = CreateH1F("h_Eff_pt_TightID_ISO04dBeta_AllMu", 
						       "h_Eff_pt_TightID_ISO04dBeta_AllMu",         200, 0, 200);
  h_Eff_pt_TightID_ISO03PFWeighted_AllMu   = CreateH1F("h_Eff_pt_TightID_ISO03PFWeighted_AllMu", 
						       "h_Eff_pt_TightID_ISO03PFWeighted_AllMu",    200, 0, 200);
  h_Eff_pt_TightID_ISO04PFWeighted_AllMu   = CreateH1F("h_Eff_pt_TightID_ISO04PFWeighted_AllMu", 
						       "h_Eff_pt_TightID_ISO04PFWeighted_AllMu",    200, 0, 200);
  h_Eff_pt_TightID_ISO03PUPPI_AllMu        = CreateH1F("h_Eff_pt_TightID_ISO03PUPPI_AllMu", 
						       "h_Eff_pt_TightID_ISO03PUPPI_AllMu",         200, 0, 200);
  h_Eff_pt_TightID_ISO04PUPPI_AllMu        = CreateH1F("h_Eff_pt_TightID_ISO04PUPPI_AllMu", 
						       "h_Eff_pt_TightID_ISO04PUPPI_AllMu",         200, 0, 200);
  h_Eff_pt_MediumID_AllMu                  = CreateH1F("h_Eff_pt_MediumID_AllMu", 
						       "h_Eff_pt_MediumID_AllMu",                   200, 0, 200);
  h_Eff_pt_MediumID_ISO03_AllMu            = CreateH1F("h_Eff_pt_MediumID_ISO03_AllMu", 
						       "h_Eff_pt_MediumID_ISO03_AllMu",             200, 0, 200);
  h_Eff_pt_MediumID_ISO04_AllMu            = CreateH1F("h_Eff_pt_MediumID_ISO04_AllMu", 
						       "h_Eff_pt_MediumID_ISO04_AllMu",             200, 0, 200);
  h_Eff_pt_MediumID_ISO03dBeta_AllMu       = CreateH1F("h_Eff_pt_MediumID_ISO03dBeta_AllMu", 
						       "h_Eff_pt_MediumID_ISO03dBeta_AllMu",        200, 0, 200);
  h_Eff_pt_MediumID_ISO04dBeta_AllMu       = CreateH1F("h_Eff_pt_MediumID_ISO04dBeta_AllMu", 
						       "h_Eff_pt_MediumID_ISO04dBeta_AllMu",        200, 0, 200);
  h_Eff_pt_MediumID_ISO03PFWeighted_AllMu  = CreateH1F("h_Eff_pt_MediumID_ISO03PFWeighted_AllMu", 
						       "h_Eff_pt_MediumID_ISO03PFWeighted_AllMu",   200, 0, 200);
  h_Eff_pt_MediumID_ISO04PFWeighted_AllMu  = CreateH1F("h_Eff_pt_MediumID_ISO04PFWeighted_AllMu", 
						       "h_Eff_pt_MediumID_ISO04PFWeighted_AllMu",   200, 0, 200);
  h_Eff_pt_MediumID_ISO03PUPPI_AllMu       = CreateH1F("h_Eff_pt_MediumID_ISO03PUPPI_AllMu", 
						       "h_Eff_pt_MediumID_ISO03PUPPI_AllMu",        200, 0, 200);
  h_Eff_pt_MediumID_ISO04PUPPI_AllMu       = CreateH1F("h_Eff_pt_MediumID_ISO04PUPPI_AllMu", 
						       "h_Eff_pt_MediumID_ISO04PUPPI_AllMu",        200, 0, 200);


  h_Eff_eta_NoID_AllMu                     = CreateH1F("h_Eff_eta_NoID_AllMu",
						       "h_Eff_eta_NoID_AllMu",                      50, -2.5, 2.5);
  h_Eff_eta_HWWID_AllMu                    = CreateH1F("h_Eff_eta_HWWID_AllMu",
						       "h_Eff_eta_HWWID_AllMu",                     50, -2.5, 2.5);
  h_Eff_eta_TightIDGoT_AllMu               = CreateH1F("h_Eff_eta_TightIDGoT_AllMu",
						       "h_Eff_eta_TightIDGoT_AllMu",                50, -2.5, 2.5);
  h_Eff_eta_TightIDipsHWW_AllMu            = CreateH1F("h_Eff_eta_TightIDipsHWW_AllMu",
						       "h_Eff_eta_TightIDipsHWW_AllMu",             50, -2.5, 2.5);
  h_Eff_eta_MediumIDipsHWW_AllMu           = CreateH1F("h_Eff_eta_MediumIDipsHWW_AllMu",
						       "h_Eff_eta_MediumIDipsHWW_AllMu",            50, -2.5, 2.5);
  h_Eff_eta_TightID_AllMu                  = CreateH1F("h_Eff_eta_TightID_AllMu", 
						       "h_Eff_eta_TightID_AllMu",                   50, -2.5, 2.5);
  h_Eff_eta_TightID_ISO03_AllMu            = CreateH1F("h_Eff_eta_TightID_ISO03_AllMu", 
						       "h_Eff_eta_TightID_ISO03_AllMu",             50, -2.5, 2.5);
  h_Eff_eta_TightID_ISO04_AllMu            = CreateH1F("h_Eff_eta_TightID_ISO04_AllMu", 
						       "h_Eff_eta_TightID_ISO04_AllMu",             50, -2.5, 2.5);
  h_Eff_eta_TightID_ISO03dBeta_AllMu       = CreateH1F("h_Eff_eta_TightID_ISO03dBeta_AllMu", 
						       "h_Eff_eta_TightID_ISO03dBeta_AllMu",        50, -2.5, 2.5);
  h_Eff_eta_TightID_ISO04dBeta_AllMu       = CreateH1F("h_Eff_eta_TightID_ISO04dBeta_AllMu", 
						       "h_Eff_eta_TightID_ISO04dBeta_AllMu",        50, -2.5, 2.5);
  h_Eff_eta_TightID_ISO03PFWeighted_AllMu  = CreateH1F("h_Eff_eta_TightID_ISO03PFWeighted_AllMu", 
						       "h_Eff_eta_TightID_ISO03PFWeighted_AllMu",   50, -2.5, 2.5);
  h_Eff_eta_TightID_ISO04PFWeighted_AllMu  = CreateH1F("h_Eff_eta_TightID_ISO04PFWeighted_AllMu", 
						       "h_Eff_eta_TightID_ISO04PFWeighted_AllMu",   50, -2.5, 2.5);
  h_Eff_eta_TightID_ISO03PUPPI_AllMu       = CreateH1F("h_Eff_eta_TightID_ISO03PUPPI_AllMu", 
						       "h_Eff_eta_TightID_ISO03PUPPI_AllMu",        50, -2.5, 2.5);
  h_Eff_eta_TightID_ISO04PUPPI_AllMu       = CreateH1F("h_Eff_eta_TightID_ISO04PUPPI_AllMu", 
						       "h_Eff_eta_TightID_ISO04PUPPI_AllMu",        50, -2.5, 2.5);
  h_Eff_eta_MediumID_AllMu                 = CreateH1F("h_Eff_eta_MediumID_AllMu", 
						       "h_Eff_eta_MediumID_AllMu",                  50, -2.5, 2.5);
  h_Eff_eta_MediumID_ISO03_AllMu           = CreateH1F("h_Eff_eta_MediumID_ISO03_AllMu", 
						       "h_Eff_eta_MediumID_ISO03_AllMu",            50, -2.5, 2.5);
  h_Eff_eta_MediumID_ISO04_AllMu           = CreateH1F("h_Eff_eta_MediumID_ISO04_AllMu", 
						       "h_Eff_eta_MediumID_ISO04_AllMu",            50, -2.5, 2.5);
  h_Eff_eta_MediumID_ISO03dBeta_AllMu      = CreateH1F("h_Eff_eta_MediumID_ISO03dBeta_AllMu", 
						       "h_Eff_eta_MediumID_ISO03dBeta_AllMu",       50, -2.5, 2.5);
  h_Eff_eta_MediumID_ISO04dBeta_AllMu      = CreateH1F("h_Eff_eta_MediumID_ISO04dBeta_AllMu", 
						       "h_Eff_eta_MediumID_ISO04dBeta_AllMu",       50, -2.5, 2.5);
  h_Eff_eta_MediumID_ISO03PFWeighted_AllMu = CreateH1F("h_Eff_eta_MediumID_ISO03PFWeighted_AllMu", 
						       "h_Eff_eta_MediumID_ISO03PFWeighted_AllMu",  50, -2.5, 2.5);
  h_Eff_eta_MediumID_ISO04PFWeighted_AllMu = CreateH1F("h_Eff_eta_MediumID_ISO04PFWeighted_AllMu", 
						       "h_Eff_eta_MediumID_ISO04PFWeighted_AllMu",  50, -2.5, 2.5);
  h_Eff_eta_MediumID_ISO03PUPPI_AllMu      = CreateH1F("h_Eff_eta_MediumID_ISO03PUPPI_AllMu", 
						       "h_Eff_eta_MediumID_ISO03PUPPI_AllMu",       50, -2.5, 2.5);
  h_Eff_eta_MediumID_ISO04PUPPI_AllMu      = CreateH1F("h_Eff_eta_MediumID_ISO04PUPPI_AllMu", 
						       "h_Eff_eta_MediumID_ISO04PUPPI_AllMu",       50, -2.5, 2.5);
  

  h_Eff_npv_NoID_AllMu                     = CreateH1F("h_Eff_npv_NoID_AllMu",
						       "h_Eff_npv_NoID_AllMu",                      45, 0, 45);
  h_Eff_npv_HWWID_AllMu                    = CreateH1F("h_Eff_npv_HWWID_AllMu",
						       "h_Eff_npv_HWWID_AllMu",                     45, 0, 45);
  h_Eff_npv_TightIDGoT_AllMu               = CreateH1F("h_Eff_npv_TightIDGoT_AllMu",
						       "h_Eff_npv_TightIDGoT_AllMu",                45, 0, 45);
  h_Eff_npv_TightIDipsHWW_AllMu            = CreateH1F("h_Eff_npv_TightIDipsHWW_AllMu",
						       "h_Eff_npv_TightIDipsHWW_AllMu",             45, 0, 45);
  h_Eff_npv_MediumIDipsHWW_AllMu           = CreateH1F("h_Eff_npv_MediumIDipsHWW_AllMu",
						       "h_Eff_npv_MediumIDipsHWW_AllMu",            45, 0, 45);
  h_Eff_npv_TightID_AllMu                  = CreateH1F("h_Eff_npv_TightID_AllMu", 
						       "h_Eff_npv_TightID_AllMu",                   45, 0, 45);
  h_Eff_npv_TightID_ISO03_AllMu            = CreateH1F("h_Eff_npv_TightID_ISO03_AllMu", 
						       "h_Eff_npv_TightID_ISO03_AllMu",             45, 0, 45);
  h_Eff_npv_TightID_ISO04_AllMu            = CreateH1F("h_Eff_npv_TightID_ISO04_AllMu", 
						       "h_Eff_npv_TightID_ISO04_AllMu",             45, 0, 45);
  h_Eff_npv_TightID_ISO03dBeta_AllMu       = CreateH1F("h_Eff_npv_TightID_ISO03dBeta_AllMu", 
						       "h_Eff_npv_TightID_ISO03dBeta_AllMu",        45, 0, 45);
  h_Eff_npv_TightID_ISO04dBeta_AllMu       = CreateH1F("h_Eff_npv_TightID_ISO04dBeta_AllMu", 
						       "h_Eff_npv_TightID_ISO04dBeta_AllMu",        45, 0, 45);
  h_Eff_npv_TightID_ISO03PFWeighted_AllMu  = CreateH1F("h_Eff_npv_TightID_ISO03PFWeighted_AllMu", 
						       "h_Eff_npv_TightID_ISO03PFWeighted_AllMu",   45, 0, 45);
  h_Eff_npv_TightID_ISO04PFWeighted_AllMu  = CreateH1F("h_Eff_npv_TightID_ISO04PFWeighted_AllMu", 
						       "h_Eff_npv_TightID_ISO04PFWeighted_AllMu",   45, 0, 45);
  h_Eff_npv_TightID_ISO03PUPPI_AllMu       = CreateH1F("h_Eff_npv_TightID_ISO03PUPPI_AllMu", 
						       "h_Eff_npv_TightID_ISO03PUPPI_AllMu",        45, 0, 45);
  h_Eff_npv_TightID_ISO04PUPPI_AllMu       = CreateH1F("h_Eff_npv_TightID_ISO04PUPPI_AllMu", 
						       "h_Eff_npv_TightID_ISO04PUPPI_AllMu",        45, 0, 45);
  h_Eff_npv_MediumID_AllMu                 = CreateH1F("h_Eff_npv_MediumID_AllMu", 
						       "h_Eff_npv_MediumID_AllMu",                  45, 0, 45);
  h_Eff_npv_MediumID_ISO03_AllMu           = CreateH1F("h_Eff_npv_MediumID_ISO03_AllMu", 
						       "h_Eff_npv_MediumID_ISO03_AllMu",            45, 0, 45);
  h_Eff_npv_MediumID_ISO04_AllMu           = CreateH1F("h_Eff_npv_MediumID_ISO04_AllMu", 
						       "h_Eff_npv_MediumID_ISO04_AllMu",            45, 0, 45);
  h_Eff_npv_MediumID_ISO03dBeta_AllMu      = CreateH1F("h_Eff_npv_MediumID_ISO03dBeta_AllMu", 
						       "h_Eff_npv_MediumID_ISO03dBeta_AllMu",       45, 0, 45);
  h_Eff_npv_MediumID_ISO04dBeta_AllMu      = CreateH1F("h_Eff_npv_MediumID_ISO04dBeta_AllMu", 
						       "h_Eff_npv_MediumID_ISO04dBeta_AllMu",       45, 0, 45);
  h_Eff_npv_MediumID_ISO03PFWeighted_AllMu = CreateH1F("h_Eff_npv_MediumID_ISO03PFWeighted_AllMu", 
						       "h_Eff_npv_MediumID_ISO03PFWeighted_AllMu",  45, 0, 45);
  h_Eff_npv_MediumID_ISO04PFWeighted_AllMu = CreateH1F("h_Eff_npv_MediumID_ISO04PFWeighted_AllMu", 
						       "h_Eff_npv_MediumID_ISO04PFWeighted_AllMu",  45, 0, 45);
  h_Eff_npv_MediumID_ISO03PUPPI_AllMu      = CreateH1F("h_Eff_npv_MediumID_ISO03PUPPI_AllMu", 
						       "h_Eff_npv_MediumID_ISO03PUPPI_AllMu",       45, 0, 45);
  h_Eff_npv_MediumID_ISO04PUPPI_AllMu      = CreateH1F("h_Eff_npv_MediumID_ISO04PUPPI_AllMu", 
						       "h_Eff_npv_MediumID_ISO04PUPPI_AllMu",       45, 0, 45);

}


void MuonEfficienciesSelector::InsideLoop() {
 
 // The InsideLoop() function is called for each entry in the tree to be processed  

  if (_Debug) std::cout << "[DEBUG][Event "<< Get<int>("T_Event_EventNumber") <<"]" << std::endl;

  //------------------------------------------------------------------------------
  // Initialise data members for each event
  //------------------------------------------------------------------------------

  // // GEN Info
  // G_GEN_PromptMuon_4vec  = *(GetParam<std::vector<TLorentzVector>*>("GEN_PromptMuon_4vec"));
  // G_GEN_isMuMu           = *(GetParam<bool*>("GEN_isMuMu"));
  // G_GEN_isMuTau          = *(GetParam<bool*>("GEN_isMuTau"));
  // G_GEN_isTauMu          = *(GetParam<bool*>("GEN_isTauMu")); 
  // G_GEN_isTauTau         = *(GetParam<bool*>("GEN_isTauTau")); 

  // // RECO muons
  // G_Muon_4vec            = *(GetParam<std::vector<TLorentzVector>*>("Muon_4vec"));

  // G_MuonID_Tight         = *(GetParam<std::vector<bool>*>("MuonID_Tight"));
  // G_MuonID_Medium        = *(GetParam<std::vector<bool>*>("MuonID_Medium"));
  // G_MuonID_HWW           = *(GetParam<std::vector<bool>*>("MuonID_HWW"));
  // G_MuonID_Tight_GoT     = *(GetParam<std::vector<bool>*>("MuonID_Tight_GoT"));
  // G_MuonID_IPs_HWW       = *(GetParam<std::vector<bool>*>("MuonID_IPs_HWW"));
  // G_MuonID_GLBorTRKArb   = *(GetParam<std::vector<bool>*>("MuonID_GLBorTRKArb"));
  // G_MuonID_Fiducial      = *(GetParam<std::vector<bool>*>("MuonID_Fiducial"));

  // G_MuonISO03            = *(GetParam<std::vector<bool>*>("MuonISO03"));
  // G_MuonISO03_dBeta      = *(GetParam<std::vector<bool>*>("MuonISO03_dBeta"));
  // G_MuonISO03_PFWeighted = *(GetParam<std::vector<bool>*>("MuonISO03_PFWeighted"));
  // G_MuonISO03_PUPPI      = *(GetParam<std::vector<bool>*>("MuonISO03_PUPPI"));
  // G_MuonISO04            = *(GetParam<std::vector<bool>*>("MuonISO04"));
  // G_MuonISO04_dBeta      = *(GetParam<std::vector<bool>*>("MuonISO04_dBeta"));
  // G_MuonISO04_PFWeighted = *(GetParam<std::vector<bool>*>("MuonISO04_PFWeighted"));
  // G_MuonISO04_PUPPI      = *(GetParam<std::vector<bool>*>("MuonISO04_PUPPI"));

  // G_Muon_ChCompatible    = *(GetParam<std::vector<bool>*>("Muon_ChCompatible"));
  // G_Muon_Matching        = *(GetParam<std::vector<int>*>("Muon_Matching"));
  
  // // Sizes
  // G_RecoMuSize  = *(GetParam<UInt_t*>("RecoMuSize"));
  // G_NPV         = *(GetParam<UInt_t*>("NPV"));

  // // Event Flags
  // EvtFlag_Fiducial = *(GetParam<bool*>("FLAG_Fiducial"));
  // EvtFlag_Gen      = *(GetParam<bool*>("FLAG_Gen"));
  // EvtFlag_Matching = *(GetParam<bool*>("FLAG_Matching"));


  //------------------------------------------------------------------------------
  // Do Efficiencies in function of Pt, Eta and NPV for muons that passed matching
  //------------------------------------------------------------------------------

  if (EvtFlag_Matching) {
    EffsSingleMu(0);
    EffsSingleMu(1);
  }

  EffsAllMu();

  
} // end inside Loop

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// MEMBER FUNCTIONS
//

//---------------------------------------------------------------------------------------------------------------------
// EffsSingleMu: Calculate ID and ISO efficiencies individually for the 2 muons passing Evt_Flag_Matching
//---------------------------------------------------------------------------------------------------------------------
void MuonEfficienciesSelector::EffsSingleMu(UInt_t iMu) {

  float pt  = Get<float>("T_Muon_Pt", iMu); 
  float eta = Get<float>("T_Muon_Eta",iMu);
  float npv = G_NPV;

  bool tight  = G_MuonID_Tight[iMu];
  bool medium = G_MuonID_Medium[iMu];
  
   h_Eff_pt_NoID[iMu]->Fill(pt);
  h_Eff_eta_NoID[iMu]->Fill(eta);
  h_Eff_npv_NoID[iMu]->Fill(npv);
  
  if (G_MuonID_HWW[iMu]) {
     h_Eff_pt_HWWID[iMu]->Fill(pt);
    h_Eff_eta_HWWID[iMu]->Fill(eta);
    h_Eff_npv_HWWID[iMu]->Fill(npv);
  }

  if (G_MuonID_Tight_GoT[iMu]) {
     h_Eff_pt_TightIDGoT[iMu]->Fill(pt);
    h_Eff_eta_TightIDGoT[iMu]->Fill(eta);
    h_Eff_npv_TightIDGoT[iMu]->Fill(npv);
  }
      
  if (tight && G_MuonID_IPs_HWW[iMu]) {
     h_Eff_pt_TightIDipsHWW[iMu]->Fill(pt);
    h_Eff_eta_TightIDipsHWW[iMu]->Fill(eta);
    h_Eff_npv_TightIDipsHWW[iMu]->Fill(npv);
  }

  if (G_MuonID_Medium[iMu] && G_MuonID_IPs_HWW[iMu]) {
     h_Eff_pt_MediumIDipsHWW[iMu]->Fill(pt);
    h_Eff_eta_MediumIDipsHWW[iMu]->Fill(eta);
    h_Eff_npv_MediumIDipsHWW[iMu]->Fill(npv);
  }

  if (tight) {
     h_Eff_pt_TightID[iMu]->Fill(pt);
    h_Eff_eta_TightID[iMu]->Fill(eta);
    h_Eff_npv_TightID[iMu]->Fill(npv);
  }

  if (tight && G_MuonISO03[iMu]) {
     h_Eff_pt_TightID_ISO03[iMu]->Fill(pt);
    h_Eff_eta_TightID_ISO03[iMu]->Fill(eta);
    h_Eff_npv_TightID_ISO03[iMu]->Fill(npv);
  }
  
  if (tight && G_MuonISO04[iMu]) {
     h_Eff_pt_TightID_ISO04[iMu]->Fill(pt);
    h_Eff_eta_TightID_ISO04[iMu]->Fill(eta);
    h_Eff_npv_TightID_ISO04[iMu]->Fill(npv);
  }
  
  if (tight && G_MuonISO03_dBeta[iMu]) {
     h_Eff_pt_TightID_ISO03dBeta[iMu]->Fill(pt);
    h_Eff_eta_TightID_ISO03dBeta[iMu]->Fill(eta);
    h_Eff_npv_TightID_ISO03dBeta[iMu]->Fill(npv);
  }
  
  if (tight && G_MuonISO04_dBeta[iMu]) {
     h_Eff_pt_TightID_ISO04dBeta[iMu]->Fill(pt);
    h_Eff_eta_TightID_ISO04dBeta[iMu]->Fill(eta);
    h_Eff_npv_TightID_ISO04dBeta[iMu]->Fill(npv);
  }
  
  if (tight && G_MuonISO03_PFWeighted[iMu]) {
     h_Eff_pt_TightID_ISO03PFWeighted[iMu]->Fill(pt);
    h_Eff_eta_TightID_ISO03PFWeighted[iMu]->Fill(eta);
    h_Eff_npv_TightID_ISO03PFWeighted[iMu]->Fill(npv);
  }
  
  if (tight && G_MuonISO04_PFWeighted[iMu]) {
     h_Eff_pt_TightID_ISO04PFWeighted[iMu]->Fill(pt);
    h_Eff_eta_TightID_ISO04PFWeighted[iMu]->Fill(eta);
    h_Eff_npv_TightID_ISO04PFWeighted[iMu]->Fill(npv);
  } 

  if (tight && G_MuonISO03_PUPPI[iMu]) {
     h_Eff_pt_TightID_ISO03PUPPI[iMu]->Fill(pt);
    h_Eff_eta_TightID_ISO03PUPPI[iMu]->Fill(eta);
    h_Eff_npv_TightID_ISO03PUPPI[iMu]->Fill(npv);
  }
  
  if (tight && G_MuonISO04_PUPPI[iMu]) {
     h_Eff_pt_TightID_ISO04PUPPI[iMu]->Fill(pt);
    h_Eff_eta_TightID_ISO04PUPPI[iMu]->Fill(eta);
    h_Eff_npv_TightID_ISO04PUPPI[iMu]->Fill(npv);
  }

  if (medium) {
     h_Eff_pt_MediumID[iMu]->Fill(pt);
    h_Eff_eta_MediumID[iMu]->Fill(eta);
    h_Eff_npv_MediumID[iMu]->Fill(npv);
  } 

  if (medium && G_MuonISO03[iMu]) {
     h_Eff_pt_MediumID_ISO03[iMu]->Fill(pt);
    h_Eff_eta_MediumID_ISO03[iMu]->Fill(eta);
    h_Eff_npv_MediumID_ISO03[iMu]->Fill(npv);
  }
  
  if (medium && G_MuonISO04[iMu]) {
     h_Eff_pt_MediumID_ISO04[iMu]->Fill(pt);
    h_Eff_eta_MediumID_ISO04[iMu]->Fill(eta);
    h_Eff_npv_MediumID_ISO04[iMu]->Fill(npv);
  }
  
  if (medium && G_MuonISO03_dBeta[iMu]) {
     h_Eff_pt_MediumID_ISO03dBeta[iMu]->Fill(pt);
    h_Eff_eta_MediumID_ISO03dBeta[iMu]->Fill(eta);
    h_Eff_npv_MediumID_ISO03dBeta[iMu]->Fill(npv);
  }
  
  if (medium && G_MuonISO04_dBeta[iMu]) {
     h_Eff_pt_MediumID_ISO04dBeta[iMu]->Fill(pt);
    h_Eff_eta_MediumID_ISO04dBeta[iMu]->Fill(eta);
    h_Eff_npv_MediumID_ISO04dBeta[iMu]->Fill(npv);
  }
  
  if (medium && G_MuonISO03_PFWeighted[iMu]) {
     h_Eff_pt_MediumID_ISO03PFWeighted[iMu]->Fill(pt);
    h_Eff_eta_MediumID_ISO03PFWeighted[iMu]->Fill(eta);
    h_Eff_npv_MediumID_ISO03PFWeighted[iMu]->Fill(npv);
  }
  
  if (medium && G_MuonISO04_PFWeighted[iMu]) {
     h_Eff_pt_MediumID_ISO04PFWeighted[iMu]->Fill(pt);
    h_Eff_eta_MediumID_ISO04PFWeighted[iMu]->Fill(eta);
    h_Eff_npv_MediumID_ISO04PFWeighted[iMu]->Fill(npv);
  } 

  if (medium && G_MuonISO03_PUPPI[iMu]) {
     h_Eff_pt_MediumID_ISO03PUPPI[iMu]->Fill(pt);
    h_Eff_eta_MediumID_ISO03PUPPI[iMu]->Fill(eta);
    h_Eff_npv_MediumID_ISO03PUPPI[iMu]->Fill(npv);
  }
  
  if (medium && G_MuonISO04_PUPPI[iMu]) {
     h_Eff_pt_MediumID_ISO04PUPPI[iMu]->Fill(pt);
    h_Eff_eta_MediumID_ISO04PUPPI[iMu]->Fill(eta);
    h_Eff_npv_MediumID_ISO04PUPPI[iMu]->Fill(npv);
  } 
  
}

//---------------------------------------------------------------------------------------------------------------------
// EffsAllMu: Calculate ID and ISO efficiencies for all muons (matched in case of DY, TTbar, WW, HWW, etc.)
//---------------------------------------------------------------------------------------------------------------------
void MuonEfficienciesSelector::EffsAllMu() {

  for (UInt_t i=0; i<G_RecoMuSize; i++) {

    float pt  = Get<float>("T_Muon_Pt", i); 
    float eta = Get<float>("T_Muon_Eta",i);
    float npv = G_NPV;

    bool tight  = G_MuonID_Tight[i];
    bool medium = G_MuonID_Medium[i];

    if (pt < 20) continue;
    if (fabs(eta) > 2.4) continue;
    if (_Signal.Contains("DY") && !G_Muon_Matching[i]) continue;

     h_Eff_pt_NoID_AllMu->Fill(pt);
    h_Eff_eta_NoID_AllMu->Fill(eta);
    h_Eff_npv_NoID_AllMu->Fill(npv);
  
    if (G_MuonID_HWW[i]) {
       h_Eff_pt_HWWID_AllMu->Fill(pt);
      h_Eff_eta_HWWID_AllMu->Fill(eta);
      h_Eff_npv_HWWID_AllMu->Fill(npv);
    }

    if (G_MuonID_Tight_GoT[i]) {
       h_Eff_pt_TightIDGoT_AllMu->Fill(pt);
      h_Eff_eta_TightIDGoT_AllMu->Fill(eta);
      h_Eff_npv_TightIDGoT_AllMu->Fill(npv);
    }
      
    if (tight && G_MuonID_IPs_HWW[i]) {
       h_Eff_pt_TightIDipsHWW_AllMu->Fill(pt);
      h_Eff_eta_TightIDipsHWW_AllMu->Fill(eta);
      h_Eff_npv_TightIDipsHWW_AllMu->Fill(npv);
    }

    if (G_MuonID_Medium[i] && G_MuonID_IPs_HWW[i]) {
       h_Eff_pt_MediumIDipsHWW_AllMu->Fill(pt);
      h_Eff_eta_MediumIDipsHWW_AllMu->Fill(eta);
      h_Eff_npv_MediumIDipsHWW_AllMu->Fill(npv);
    }
    
    if (tight) {
       h_Eff_pt_TightID_AllMu->Fill(pt);
      h_Eff_eta_TightID_AllMu->Fill(eta);
      h_Eff_npv_TightID_AllMu->Fill(npv);
    }
    
    if (medium) {
       h_Eff_pt_MediumID_AllMu->Fill(pt);
      h_Eff_eta_MediumID_AllMu->Fill(eta);
      h_Eff_npv_MediumID_AllMu->Fill(npv);
    }
    
    if (tight && G_MuonISO03[i]) {
       h_Eff_pt_TightID_ISO03_AllMu->Fill(pt);
      h_Eff_eta_TightID_ISO03_AllMu->Fill(eta);
      h_Eff_npv_TightID_ISO03_AllMu->Fill(npv);
    }
  
    if (tight && G_MuonISO04[i]) {
       h_Eff_pt_TightID_ISO04_AllMu->Fill(pt);
      h_Eff_eta_TightID_ISO04_AllMu->Fill(eta);
      h_Eff_npv_TightID_ISO04_AllMu->Fill(npv);
    }
  
    if (tight && G_MuonISO03_dBeta[i]) {
       h_Eff_pt_TightID_ISO03dBeta_AllMu->Fill(pt);
      h_Eff_eta_TightID_ISO03dBeta_AllMu->Fill(eta);
      h_Eff_npv_TightID_ISO03dBeta_AllMu->Fill(npv);
    }
  
    if (tight && G_MuonISO04_dBeta[i]) {
       h_Eff_pt_TightID_ISO04dBeta_AllMu->Fill(pt);
      h_Eff_eta_TightID_ISO04dBeta_AllMu->Fill(eta);
      h_Eff_npv_TightID_ISO04dBeta_AllMu->Fill(npv);
    }
  
    if (tight && G_MuonISO03_PFWeighted[i]) {
       h_Eff_pt_TightID_ISO03PFWeighted_AllMu->Fill(pt);
      h_Eff_eta_TightID_ISO03PFWeighted_AllMu->Fill(eta);
      h_Eff_npv_TightID_ISO03PFWeighted_AllMu->Fill(npv);
    }
  
    if (tight && G_MuonISO04_PFWeighted[i]) {
       h_Eff_pt_TightID_ISO04PFWeighted_AllMu->Fill(pt);
      h_Eff_eta_TightID_ISO04PFWeighted_AllMu->Fill(eta);
      h_Eff_npv_TightID_ISO04PFWeighted_AllMu->Fill(npv);
    } 

    if (tight && G_MuonISO03_PUPPI[i]) {
       h_Eff_pt_TightID_ISO03PUPPI_AllMu->Fill(pt);
      h_Eff_eta_TightID_ISO03PUPPI_AllMu->Fill(eta);
      h_Eff_npv_TightID_ISO03PUPPI_AllMu->Fill(npv);
    }
  
    if (tight && G_MuonISO04_PUPPI[i]) {
       h_Eff_pt_TightID_ISO04PUPPI_AllMu->Fill(pt);
      h_Eff_eta_TightID_ISO04PUPPI_AllMu->Fill(eta);
      h_Eff_npv_TightID_ISO04PUPPI_AllMu->Fill(npv);
    } 

    if (medium && G_MuonISO03[i]) {
       h_Eff_pt_MediumID_ISO03_AllMu->Fill(pt);
      h_Eff_eta_MediumID_ISO03_AllMu->Fill(eta);
      h_Eff_npv_MediumID_ISO03_AllMu->Fill(npv);
    }
  
    if (medium && G_MuonISO04[i]) {
       h_Eff_pt_MediumID_ISO04_AllMu->Fill(pt);
      h_Eff_eta_MediumID_ISO04_AllMu->Fill(eta);
      h_Eff_npv_MediumID_ISO04_AllMu->Fill(npv);
    }
  
    if (medium && G_MuonISO03_dBeta[i]) {
       h_Eff_pt_MediumID_ISO03dBeta_AllMu->Fill(pt);
      h_Eff_eta_MediumID_ISO03dBeta_AllMu->Fill(eta);
      h_Eff_npv_MediumID_ISO03dBeta_AllMu->Fill(npv);
    }
  
    if (medium && G_MuonISO04_dBeta[i]) {
       h_Eff_pt_MediumID_ISO04dBeta_AllMu->Fill(pt);
      h_Eff_eta_MediumID_ISO04dBeta_AllMu->Fill(eta);
      h_Eff_npv_MediumID_ISO04dBeta_AllMu->Fill(npv);
    }
  
    if (medium && G_MuonISO03_PFWeighted[i]) {
       h_Eff_pt_MediumID_ISO03PFWeighted_AllMu->Fill(pt);
      h_Eff_eta_MediumID_ISO03PFWeighted_AllMu->Fill(eta);
      h_Eff_npv_MediumID_ISO03PFWeighted_AllMu->Fill(npv);
    }
  
    if (medium && G_MuonISO04_PFWeighted[i]) {
       h_Eff_pt_MediumID_ISO04PFWeighted_AllMu->Fill(pt);
      h_Eff_eta_MediumID_ISO04PFWeighted_AllMu->Fill(eta);
      h_Eff_npv_MediumID_ISO04PFWeighted_AllMu->Fill(npv);
    } 

    if (medium && G_MuonISO03_PUPPI[i]) {
       h_Eff_pt_MediumID_ISO03PUPPI_AllMu->Fill(pt);
      h_Eff_eta_MediumID_ISO03PUPPI_AllMu->Fill(eta);
      h_Eff_npv_MediumID_ISO03PUPPI_AllMu->Fill(npv);
    }
  
    if (medium && G_MuonISO04_PUPPI[i]) {
       h_Eff_pt_MediumID_ISO04PUPPI_AllMu->Fill(pt);
      h_Eff_eta_MediumID_ISO04PUPPI_AllMu->Fill(eta);
      h_Eff_npv_MediumID_ISO04PUPPI_AllMu->Fill(npv);
    } 

  }
  
}


void MuonEfficienciesSelector::Summary() {
  // Get Data Members at the client-master (after finishing the analysis at the workers nodes)
  // Only data members set here will be accesible at the client-master

  // Single Muon ID and ISO efficiencies vs pt, eta and npv
  h_Eff_pt_NoID[0]                      = FindOutput<TH1F*>("h_Eff_pt_NoID_Mu1");
  h_Eff_pt_NoID[1]                      = FindOutput<TH1F*>("h_Eff_pt_NoID_Mu2");
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

  h_Eff_eta_NoID[0]                     = FindOutput<TH1F*>("h_Eff_eta_NoID_Mu1");
  h_Eff_eta_NoID[1]                     = FindOutput<TH1F*>("h_Eff_eta_NoID_Mu2");
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

  h_Eff_npv_NoID[0]                     = FindOutput<TH1F*>("h_Eff_npv_NoID_Mu1");
  h_Eff_npv_NoID[1]                     = FindOutput<TH1F*>("h_Eff_npv_NoID_Mu2");
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

  // ISO efficiencies vs pt, eta and npv for all muons
  h_Eff_pt_NoID_AllMu                      = FindOutput<TH1F*>("h_Eff_pt_NoID_AllMu");
  h_Eff_pt_HWWID_AllMu                     = FindOutput<TH1F*>("h_Eff_pt_HWWID_AllMu");
  h_Eff_pt_TightIDGoT_AllMu                = FindOutput<TH1F*>("h_Eff_pt_TightIDGoT_AllMu");
  h_Eff_pt_TightIDipsHWW_AllMu             = FindOutput<TH1F*>("h_Eff_pt_TightIDipsHWW_AllMu");
  h_Eff_pt_MediumIDipsHWW_AllMu            = FindOutput<TH1F*>("h_Eff_pt_MediumIDipsHWW_AllMu");
  h_Eff_pt_TightID_AllMu                   = FindOutput<TH1F*>("h_Eff_pt_TightID_AllMu");
  h_Eff_pt_TightID_ISO03_AllMu             = FindOutput<TH1F*>("h_Eff_pt_TightID_ISO03_AllMu");
  h_Eff_pt_TightID_ISO04_AllMu             = FindOutput<TH1F*>("h_Eff_pt_TightID_ISO04_AllMu");
  h_Eff_pt_TightID_ISO03dBeta_AllMu        = FindOutput<TH1F*>("h_Eff_pt_TightID_ISO03dBeta_AllMu");
  h_Eff_pt_TightID_ISO04dBeta_AllMu        = FindOutput<TH1F*>("h_Eff_pt_TightID_ISO04dBeta_AllMu");
  h_Eff_pt_TightID_ISO03PFWeighted_AllMu   = FindOutput<TH1F*>("h_Eff_pt_TightID_ISO03PFWeighted_AllMu");
  h_Eff_pt_TightID_ISO04PFWeighted_AllMu   = FindOutput<TH1F*>("h_Eff_pt_TightID_ISO04PFWeighted_AllMu");
  h_Eff_pt_TightID_ISO03PUPPI_AllMu        = FindOutput<TH1F*>("h_Eff_pt_TightID_ISO03PUPPI_AllMu");
  h_Eff_pt_TightID_ISO04PUPPI_AllMu        = FindOutput<TH1F*>("h_Eff_pt_TightID_ISO04PUPPI_AllMu");
  h_Eff_pt_MediumID_AllMu                  = FindOutput<TH1F*>("h_Eff_pt_MediumID_AllMu");
  h_Eff_pt_MediumID_ISO03_AllMu            = FindOutput<TH1F*>("h_Eff_pt_MediumID_ISO03_AllMu");
  h_Eff_pt_MediumID_ISO04_AllMu            = FindOutput<TH1F*>("h_Eff_pt_MediumID_ISO04_AllMu");
  h_Eff_pt_MediumID_ISO03dBeta_AllMu       = FindOutput<TH1F*>("h_Eff_pt_MediumID_ISO03dBeta_AllMu");
  h_Eff_pt_MediumID_ISO04dBeta_AllMu       = FindOutput<TH1F*>("h_Eff_pt_MediumID_ISO04dBeta_AllMu");
  h_Eff_pt_MediumID_ISO03PFWeighted_AllMu  = FindOutput<TH1F*>("h_Eff_pt_MediumID_ISO03PFWeighted_AllMu");
  h_Eff_pt_MediumID_ISO04PFWeighted_AllMu  = FindOutput<TH1F*>("h_Eff_pt_MediumID_ISO04PFWeighted_AllMu");
  h_Eff_pt_MediumID_ISO03PUPPI_AllMu       = FindOutput<TH1F*>("h_Eff_pt_MediumID_ISO03PUPPI_AllMu");
  h_Eff_pt_MediumID_ISO04PUPPI_AllMu       = FindOutput<TH1F*>("h_Eff_pt_MediumID_ISO04PUPPI_AllMu");

  h_Eff_eta_NoID_AllMu                     = FindOutput<TH1F*>("h_Eff_eta_NoID_AllMu");
  h_Eff_eta_HWWID_AllMu                    = FindOutput<TH1F*>("h_Eff_eta_HWWID_AllMu");
  h_Eff_eta_TightIDGoT_AllMu               = FindOutput<TH1F*>("h_Eff_eta_TightIDGoT_AllMu");
  h_Eff_eta_TightIDipsHWW_AllMu            = FindOutput<TH1F*>("h_Eff_eta_TightIDipsHWW_AllMu");
  h_Eff_eta_MediumIDipsHWW_AllMu           = FindOutput<TH1F*>("h_Eff_eta_MediumIDipsHWW_AllMu");
  h_Eff_eta_TightID_AllMu                  = FindOutput<TH1F*>("h_Eff_eta_TightID_AllMu");
  h_Eff_eta_TightID_ISO03_AllMu            = FindOutput<TH1F*>("h_Eff_eta_TightID_ISO03_AllMu");
  h_Eff_eta_TightID_ISO04_AllMu            = FindOutput<TH1F*>("h_Eff_eta_TightID_ISO04_AllMu");
  h_Eff_eta_TightID_ISO03dBeta_AllMu       = FindOutput<TH1F*>("h_Eff_eta_TightID_ISO03dBeta_AllMu");
  h_Eff_eta_TightID_ISO04dBeta_AllMu       = FindOutput<TH1F*>("h_Eff_eta_TightID_ISO04dBeta_AllMu");
  h_Eff_eta_TightID_ISO03PFWeighted_AllMu  = FindOutput<TH1F*>("h_Eff_eta_TightID_ISO03PFWeighted_AllMu");
  h_Eff_eta_TightID_ISO04PFWeighted_AllMu  = FindOutput<TH1F*>("h_Eff_eta_TightID_ISO04PFWeighted_AllMu");
  h_Eff_eta_TightID_ISO03PUPPI_AllMu       = FindOutput<TH1F*>("h_Eff_eta_TightID_ISO03PUPPI_AllMu");
  h_Eff_eta_TightID_ISO04PUPPI_AllMu       = FindOutput<TH1F*>("h_Eff_eta_TightID_ISO04PUPPI_AllMu");
  h_Eff_eta_MediumID_AllMu                 = FindOutput<TH1F*>("h_Eff_eta_MediumID_AllMu");
  h_Eff_eta_MediumID_ISO03_AllMu           = FindOutput<TH1F*>("h_Eff_eta_MediumID_ISO03_AllMu");
  h_Eff_eta_MediumID_ISO04_AllMu           = FindOutput<TH1F*>("h_Eff_eta_MediumID_ISO04_AllMu");
  h_Eff_eta_MediumID_ISO03dBeta_AllMu      = FindOutput<TH1F*>("h_Eff_eta_MediumID_ISO03dBeta_AllMu");
  h_Eff_eta_MediumID_ISO04dBeta_AllMu      = FindOutput<TH1F*>("h_Eff_eta_MediumID_ISO04dBeta_AllMu");
  h_Eff_eta_MediumID_ISO03PFWeighted_AllMu = FindOutput<TH1F*>("h_Eff_eta_MediumID_ISO03PFWeighted_AllMu");
  h_Eff_eta_MediumID_ISO04PFWeighted_AllMu = FindOutput<TH1F*>("h_Eff_eta_MediumID_ISO04PFWeighted_AllMu");
  h_Eff_eta_MediumID_ISO03PUPPI_AllMu      = FindOutput<TH1F*>("h_Eff_eta_MediumID_ISO03PUPPI_AllMu");
  h_Eff_eta_MediumID_ISO04PUPPI_AllMu      = FindOutput<TH1F*>("h_Eff_eta_MediumID_ISO04PUPPI_AllMu");

  h_Eff_npv_NoID_AllMu                     = FindOutput<TH1F*>("h_Eff_npv_NoID_AllMu");
  h_Eff_npv_HWWID_AllMu                    = FindOutput<TH1F*>("h_Eff_npv_HWWID_AllMu");
  h_Eff_npv_TightIDGoT_AllMu               = FindOutput<TH1F*>("h_Eff_npv_TightIDGoT_AllMu");
  h_Eff_npv_TightIDipsHWW_AllMu            = FindOutput<TH1F*>("h_Eff_npv_TightIDipsHWW_AllMu");
  h_Eff_npv_MediumIDipsHWW_AllMu           = FindOutput<TH1F*>("h_Eff_npv_MediumIDipsHWW_AllMu");
  h_Eff_npv_TightID_AllMu                  = FindOutput<TH1F*>("h_Eff_npv_TightID_AllMu");
  h_Eff_npv_TightID_ISO03_AllMu            = FindOutput<TH1F*>("h_Eff_npv_TightID_ISO03_AllMu");
  h_Eff_npv_TightID_ISO04_AllMu            = FindOutput<TH1F*>("h_Eff_npv_TightID_ISO04_AllMu");
  h_Eff_npv_TightID_ISO03dBeta_AllMu       = FindOutput<TH1F*>("h_Eff_npv_TightID_ISO03dBeta_AllMu");
  h_Eff_npv_TightID_ISO04dBeta_AllMu       = FindOutput<TH1F*>("h_Eff_npv_TightID_ISO04dBeta_AllMu");
  h_Eff_npv_TightID_ISO03PFWeighted_AllMu  = FindOutput<TH1F*>("h_Eff_npv_TightID_ISO03PFWeighted_AllMu");
  h_Eff_npv_TightID_ISO04PFWeighted_AllMu  = FindOutput<TH1F*>("h_Eff_npv_TightID_ISO04PFWeighted_AllMu");
  h_Eff_npv_TightID_ISO03PUPPI_AllMu       = FindOutput<TH1F*>("h_Eff_npv_TightID_ISO03PUPPI_AllMu");
  h_Eff_npv_TightID_ISO04PUPPI_AllMu       = FindOutput<TH1F*>("h_Eff_npv_TightID_ISO04PUPPI_AllMu");
  h_Eff_npv_MediumID_AllMu                 = FindOutput<TH1F*>("h_Eff_npv_MediumID_AllMu");
  h_Eff_npv_MediumID_ISO03_AllMu           = FindOutput<TH1F*>("h_Eff_npv_MediumID_ISO03_AllMu");
  h_Eff_npv_MediumID_ISO04_AllMu           = FindOutput<TH1F*>("h_Eff_npv_MediumID_ISO04_AllMu");
  h_Eff_npv_MediumID_ISO03dBeta_AllMu      = FindOutput<TH1F*>("h_Eff_npv_MediumID_ISO03dBeta_AllMu");
  h_Eff_npv_MediumID_ISO04dBeta_AllMu      = FindOutput<TH1F*>("h_Eff_npv_MediumID_ISO04dBeta_AllMu");
  h_Eff_npv_MediumID_ISO03PFWeighted_AllMu = FindOutput<TH1F*>("h_Eff_npv_MediumID_ISO03PFWeighted_AllMu");
  h_Eff_npv_MediumID_ISO04PFWeighted_AllMu = FindOutput<TH1F*>("h_Eff_npv_MediumID_ISO04PFWeighted_AllMu");
  h_Eff_npv_MediumID_ISO03PUPPI_AllMu      = FindOutput<TH1F*>("h_Eff_npv_MediumID_ISO03PUPPI_AllMu");
  h_Eff_npv_MediumID_ISO04PUPPI_AllMu      = FindOutput<TH1F*>("h_Eff_npv_MediumID_ISO04PUPPI_AllMu");


}
