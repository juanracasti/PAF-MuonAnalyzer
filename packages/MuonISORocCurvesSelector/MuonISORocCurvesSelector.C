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

#include "MuonISORocCurvesSelector.h"

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

ClassImp(MuonISORocCurvesSelector)

// Initialise input parameters and data members for all events
void MuonISORocCurvesSelector::Initialise() {

  _Signal     = GetParam<TString>("Signal");
  _IsDATA     = GetParam<bool>("IsDATA");
  _NEvents    = GetParam<int>("NEvents");
  _WhichRun   = GetParam<int>("WhichRun"); 
  _Debug      = GetParam<bool>("Debug");

 
//------------------------------------------------------------------------------
// Create histos
//------------------------------------------------------------------------------

  //ISO ROC Curves
  h_RC_TightID_ISO03_AllMu            = CreateH1F("h_RC_TightID_ISO03_AllMu", 
						  "h_RC_TightID_ISO03_AllMu",             25, 0, 25);
  h_RC_TightID_ISO03_AllMu->TH1::SetDefaultSumw2();
  h_RC_TightID_ISO04_AllMu            = CreateH1F("h_RC_TightID_ISO04_AllMu", 
						  "h_RC_TightID_ISO04_AllMu",             25, 0, 25);
  h_RC_TightID_ISO03dBeta_AllMu       = CreateH1F("h_RC_TightID_ISO03dBeta_AllMu", 
						  "h_RC_TightID_ISO03dBeta_AllMu",        25, 0, 25);
  h_RC_TightID_ISO04dBeta_AllMu       = CreateH1F("h_RC_TightID_ISO04dBeta_AllMu", 
						  "h_RC_TightID_ISO04dBeta_AllMu",        25, 0, 25);
  h_RC_TightID_ISO03PFWeighted_AllMu  = CreateH1F("h_RC_TightID_ISO03PFWeighted_AllMu", 
						  "h_RC_TightID_ISO03PFWeighted_AllMu",   25, 0, 25);
  h_RC_TightID_ISO04PFWeighted_AllMu  = CreateH1F("h_RC_TightID_ISO04PFWeighted_AllMu", 
						  "h_RC_TightID_ISO04PFWeighted_AllMu",   25, 0, 25);
  h_RC_TightID_ISO03PUPPI_AllMu       = CreateH1F("h_RC_TightID_ISO03PUPPI_AllMu", 
						  "h_RC_TightID_ISO03PUPPI_AllMu",        25, 0, 25);
  h_RC_TightID_ISO04PUPPI_AllMu       = CreateH1F("h_RC_TightID_ISO04PUPPI_AllMu", 
						  "h_RC_TightID_ISO04PUPPI_AllMu",        25, 0, 25);
  h_RC_MediumID_ISO03_AllMu           = CreateH1F("h_RC_MediumID_ISO03_AllMu", 
						  "h_RC_MediumID_ISO03_AllMu",            25, 0, 25);
  h_RC_MediumID_ISO04_AllMu           = CreateH1F("h_RC_MediumID_ISO04_AllMu", 
						  "h_RC_MediumID_ISO04_AllMu",            25, 0, 25);
  h_RC_MediumID_ISO03dBeta_AllMu      = CreateH1F("h_RC_MediumID_ISO03dBeta_AllMu", 
						  "h_RC_MediumID_ISO03dBeta_AllMu",       25, 0, 25);
  h_RC_MediumID_ISO04dBeta_AllMu      = CreateH1F("h_RC_MediumID_ISO04dBeta_AllMu", 
						  "h_RC_MediumID_ISO04dBeta_AllMu",       25, 0, 25);
  h_RC_MediumID_ISO03PFWeighted_AllMu = CreateH1F("h_RC_MediumID_ISO03PFWeighted_AllMu", 
						  "h_RC_MediumID_ISO03PFWeighted_AllMu",  25, 0, 25);
  h_RC_MediumID_ISO04PFWeighted_AllMu = CreateH1F("h_RC_MediumID_ISO04PFWeighted_AllMu", 
						  "h_RC_MediumID_ISO04PFWeighted_AllMu",  25, 0, 25);
  h_RC_MediumID_ISO03PUPPI_AllMu      = CreateH1F("h_RC_MediumID_ISO03PUPPI_AllMu", 
						  "h_RC_MediumID_ISO03PUPPI_AllMu",       25, 0, 25);
  h_RC_MediumID_ISO04PUPPI_AllMu      = CreateH1F("h_RC_MediumID_ISO04PUPPI_AllMu", 
						  "h_RC_MediumID_ISO04PUPPI_AllMu",       25, 0, 25); 

}


void MuonISORocCurvesSelector::InsideLoop() {
 
 // The InsideLoop() function is called for each entry in the tree to be processed  

  if (_Debug) std::cout << "[DEBUG][Event "<< Get<int>("T_Event_EventNumber") <<"]" << std::endl;

  //------------------------------------------------------------------------------
  // Initialise data members for each event
  //------------------------------------------------------------------------------

  // // RECO muons
  // G_Muon_4vec            = *(GetParam<std::vector<TLorentzVector>*>("Muon_4vec"));

  // G_MuonID_Tight         = *(GetParam<std::vector<bool>*>("MuonID_Tight"));
  // G_MuonID_Medium        = *(GetParam<std::vector<bool>*>("MuonID_Medium"));
  // G_MuonID_HWW           = *(GetParam<std::vector<bool>*>("MuonID_HWW"));
  // G_MuonID_Tight_GoT     = *(GetParam<std::vector<bool>*>("MuonID_Tight_GoT"));
  // G_MuonID_IPs_HWW       = *(GetParam<std::vector<bool>*>("MuonID_IPs_HWW"));
  // G_MuonID_GLBorTRKArb   = *(GetParam<std::vector<bool>*>("MuonID_GLBorTRKArb"));
  // G_MuonID_Fiducial      = *(GetParam<std::vector<bool>*>("MuonID_Fiducial"));

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
  // Do Isolation Roc Curves
  //------------------------------------------------------------------------------

  ISORocCurve(); //Warning! Long calculation


} // end inside Loop

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// MEMBER FUNCTIONS
//

//---------------------------------------------------------------------------------------------------------------------
// passISO: return true if the RECO muon of index 'iMu' passes the PF Relative Isolation indicated by 'typeIso'
//          with a working point 'wp', with the help of getISO() member function
//---------------------------------------------------------------------------------------------------------------------
bool MuonISORocCurvesSelector::passISO(UInt_t iMu, TString typeIso, float wp) {

  bool passIso = false;

  float PFRelIsoBeta = getISO(iMu, typeIso);
     
  if (PFRelIsoBeta <  wp)  passIso = true;	  

  return passIso;

}

//---------------------------------------------------------------------------------------------------------------------
// getISO: return, for the RECO muon of index 'iMu', the PF Relative Isolation indicated by 'typeIso' 
//---------------------------------------------------------------------------------------------------------------------
float MuonISORocCurvesSelector::getISO(UInt_t iMu, TString typeIso) {
  
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

  else if (Exists("T_Muon_neutralIsoPUPPIR03") && typeIso == "PUPPIR03") // PF Rel. ISO, dR=0.3 and PUPPI corrections
    PFRelIso = ( Get<float>("T_Muon_chargedHadronIsoR03",iMu) + Get<float>("T_Muon_neutralIsoPUPPIR03",iMu) )
      / pt;

  else if (Exists("T_Muon_neutralIsoPUPPIR04") && typeIso == "PUPPIR04") // PF Rel. ISO, dR=0.4 and PUPPI corrections
    PFRelIso = ( Get<float>("T_Muon_chargedHadronIsoR04",iMu) + Get<float>("T_Muon_neutralIsoPUPPIR04",iMu) )
      / pt;

  return PFRelIso;

}

//---------------------------------------------------------------------------------------------------------------------
// ISORocCurve: Fill histos to calculate Isolation Roc Curves for all muons
//---------------------------------------------------------------------------------------------------------------------
void MuonISORocCurvesSelector::ISORocCurve() {

  for (UInt_t i=0; i<G_RecoMuSize; i++) {

    float pt  = Get<float>("T_Muon_Pt", i); 
    float eta = Get<float>("T_Muon_Eta",i);

    if (pt < 20) continue;
    if (fabs(eta) > 2.4) continue;
    if (_Signal.Contains("DY") && !G_Muon_Matching[i]) continue;

    if (G_MuonID_Tight[i]) {

      h_RC_TightID_ISO03_AllMu          ->Fill(-1);
      h_RC_TightID_ISO04_AllMu          ->Fill(-1);
      h_RC_TightID_ISO03dBeta_AllMu     ->Fill(-1);
      h_RC_TightID_ISO04dBeta_AllMu     ->Fill(-1);
      h_RC_TightID_ISO03PFWeighted_AllMu->Fill(-1);
      h_RC_TightID_ISO04PFWeighted_AllMu->Fill(-1);
      h_RC_TightID_ISO03PUPPI_AllMu     ->Fill(-1);
      h_RC_TightID_ISO04PUPPI_AllMu     ->Fill(-1);

      for (int j = 0; j < 25; ++j) {

	if (passISO(i, "R03",           j*2.0/100.0)) h_RC_TightID_ISO03_AllMu          ->Fill(j);
	if (passISO(i, "dBetaR03",      j*2.0/100.0)) h_RC_TightID_ISO03dBeta_AllMu     ->Fill(j);
	if (passISO(i, "PFWeightedR03", j*2.0/100.0)) h_RC_TightID_ISO03PFWeighted_AllMu->Fill(j);
	if (passISO(i, "PUPPIR03",      j*2.0/100.0)) h_RC_TightID_ISO03PUPPI_AllMu     ->Fill(j);
	if (passISO(i, "R04",           j*2.0/100.0)) h_RC_TightID_ISO04_AllMu          ->Fill(j);
	if (passISO(i, "dBetaR04",      j*2.0/100.0)) h_RC_TightID_ISO04dBeta_AllMu     ->Fill(j);
	if (passISO(i, "PFWeightedR04", j*2.0/100.0)) h_RC_TightID_ISO04PFWeighted_AllMu->Fill(j);
	if (passISO(i, "PUPPIR04",      j*2.0/100.0)) h_RC_TightID_ISO04PUPPI_AllMu     ->Fill(j);

      }

    }

    if (G_MuonID_Medium[i]) {

      h_RC_MediumID_ISO03_AllMu          ->Fill(-1);
      h_RC_MediumID_ISO04_AllMu          ->Fill(-1);
      h_RC_MediumID_ISO03dBeta_AllMu     ->Fill(-1);
      h_RC_MediumID_ISO04dBeta_AllMu     ->Fill(-1);
      h_RC_MediumID_ISO03PFWeighted_AllMu->Fill(-1);
      h_RC_MediumID_ISO04PFWeighted_AllMu->Fill(-1);
      h_RC_MediumID_ISO03PUPPI_AllMu     ->Fill(-1);
      h_RC_MediumID_ISO04PUPPI_AllMu     ->Fill(-1);

      for (int j = 0; j < 25; ++j) {

	if (passISO(i, "R03",           j*2.0/100.0)) h_RC_MediumID_ISO03_AllMu          ->Fill(j);
	if (passISO(i, "dBetaR03",      j*2.0/100.0)) h_RC_MediumID_ISO03dBeta_AllMu     ->Fill(j);
	if (passISO(i, "PFWeightedR03", j*2.0/100.0)) h_RC_MediumID_ISO03PFWeighted_AllMu->Fill(j);
	if (passISO(i, "PUPPIR03",      j*2.0/100.0)) h_RC_MediumID_ISO03PUPPI_AllMu     ->Fill(j);
	if (passISO(i, "R04",           j*2.0/100.0)) h_RC_MediumID_ISO04_AllMu          ->Fill(j);
	if (passISO(i, "dBetaR04",      j*2.0/100.0)) h_RC_MediumID_ISO04dBeta_AllMu     ->Fill(j);
	if (passISO(i, "PFWeightedR04", j*2.0/100.0)) h_RC_MediumID_ISO04PFWeighted_AllMu->Fill(j);
	if (passISO(i, "PUPPIR04",      j*2.0/100.0)) h_RC_MediumID_ISO04PUPPI_AllMu     ->Fill(j);

      }

    }

  }
  
}

void MuonISORocCurvesSelector::Summary() {
  // Get Data Members at the client-master (after finishing the analysis at the workers nodes)
  // Only data members set here will be accesible at the client-master


  //ISO ROC Curves
  h_RC_TightID_ISO03_AllMu            = FindOutput<TH1F*>("h_RC_TightID_ISO03_AllMu");
  h_RC_TightID_ISO04_AllMu            = FindOutput<TH1F*>("h_RC_TightID_ISO04_AllMu");
  h_RC_TightID_ISO03dBeta_AllMu       = FindOutput<TH1F*>("h_RC_TightID_ISO03dBeta_AllMu");
  h_RC_TightID_ISO04dBeta_AllMu       = FindOutput<TH1F*>("h_RC_TightID_ISO04dBeta_AllMu");
  h_RC_TightID_ISO03PFWeighted_AllMu  = FindOutput<TH1F*>("h_RC_TightID_ISO03PFWeighted_AllMu");
  h_RC_TightID_ISO04PFWeighted_AllMu  = FindOutput<TH1F*>("h_RC_TightID_ISO04PFWeighted_AllMu");
  h_RC_TightID_ISO03PUPPI_AllMu       = FindOutput<TH1F*>("h_RC_TightID_ISO03PUPPI_AllMu");
  h_RC_TightID_ISO04PUPPI_AllMu       = FindOutput<TH1F*>("h_RC_TightID_ISO04PUPPI_AllMu");
  h_RC_MediumID_ISO03_AllMu           = FindOutput<TH1F*>("h_RC_MediumID_ISO03_AllMu");
  h_RC_MediumID_ISO04_AllMu           = FindOutput<TH1F*>("h_RC_MediumID_ISO04_AllMu");
  h_RC_MediumID_ISO03dBeta_AllMu      = FindOutput<TH1F*>("h_RC_MediumID_ISO03dBeta_AllMu");
  h_RC_MediumID_ISO04dBeta_AllMu      = FindOutput<TH1F*>("h_RC_MediumID_ISO04dBeta_AllMu");
  h_RC_MediumID_ISO03PFWeighted_AllMu = FindOutput<TH1F*>("h_RC_MediumID_ISO03PFWeighted_AllMu");
  h_RC_MediumID_ISO04PFWeighted_AllMu = FindOutput<TH1F*>("h_RC_MediumID_ISO04PFWeighted_AllMu");
  h_RC_MediumID_ISO03PUPPI_AllMu      = FindOutput<TH1F*>("h_RC_MediumID_ISO03PUPPI_AllMu");
  h_RC_MediumID_ISO04PUPPI_AllMu      = FindOutput<TH1F*>("h_RC_MediumID_ISO04PUPPI_AllMu");


}
