////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

/////            Muon Studies for CSA14                     ////

/////        A. Calderón (IFCA)   18 / 08 / 2014            ////

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////


#include "BasicmuonAnalyzer.h"

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

ClassImp(BasicmuonAnalyzer)

void BasicmuonAnalyzer::Initialise() {

  _Signal = GetParam<TString>("Signal");
  _IsDATA = GetParam<bool>("IsDATA");
  _NEvents = GetParam<int>("NEvents");
  _Luminosity = GetParam<float>("Luminosity");
  _XSection = GetParam<float>("XSection");
  _WhichRun = GetParam<int>("WhichRun");  

  //For debugging
  G_Debug_DefineAnalysisVariables = true;

 
//------------------------------------------------------------------------------
// Create histos
//------------------------------------------------------------------------------
 

  h_N_PV  = CreateH1F ("h_N_PV","h_N_PV",50,0,50); 
  h_N_PV->TH1::SetDefaultSumw2();
  h_N_PV2 = CreateH1F ("h_N_PV2","h_N_PV2",50,0,50); 
  h_N_PV3 = CreateH1F ("h_N_PV3","h_N_PV3",50,0,50); 


}


void BasicmuonAnalyzer::InsideLoop() {
 

 // The InsideLoop() function is called for each entry in the tree to be processed  


//------------------------------------------------------------------------------
// Define weights
//------------------------------------------------------------------------------
 
  //float pileupweight = 1.;

  //  if (!IsDATA)
    // pileupweight = fPUWeight->GetWeight(T_Event_nPU);

if (G_Debug_DefineAnalysisVariables) std::cout << "[DEBUG][Event "<< 
						 Get<int>("T_Event_EventNumber") <<"] Event start" << std::endl;

  float factN = 1;
  
  if (_XSection > 0) factN = _XSection * _Luminosity / _NEvents;

  //factN = factN*pileupweight;


  //------------------------------------------------------------------------------
  // Init variables
  //------------------------------------------------------------------------------

  ///** GEN INFORMATION
  G_GEN_PromptMuon_4vec.clear();
  G_GEN_Muon_4vec.clear();

  ///** MUONS
  G_Muon_4vec.clear();

  G_MuonID_Tight.clear();
  G_MuonID_Medium.clear();
  G_MuonID_HWW.clear();
  G_MuonID_MediumHWW.clear();
  G_MuonID_GLBorTRKArb.clear();
  G_MuonID_Fiducial.clear();

  G_MuonISO03.clear();
  G_MuonISO03_dBeta.clear();
  G_MuonISO03_PFWeights.clear();
  G_MuonISO04.clear();
  G_MuonISO04_dBeta.clear();
  G_MuonISO04_PFWeights.clear();

  G_isMuMu      = false;
  G_isMuTau     = false;
  G_isTauMu     = false;
  G_isTauTau    = false;
  G_isNonPrompt = false;


//------------------------------------------------------------------------------
// Get all RECO muons
//------------------------------------------------------------------------------

  UInt_t muonSize = 0;

  muonSize = Get<std::vector<float>*>("T_Muon_Px")->size();
   
  if ( muonSize > 0 ) {  // asking for at least one muon in the event 

    for (unsigned int i = 0; i < muonSize; i++) {
      
      //-->define the global 4D momentum for each muon. 
      G_Muon_4vec.push_back(TLorentzVector(Get<float>("T_Muon_Px",i), Get<float>("T_Muon_Py",i), 
					   Get<float>("T_Muon_Pz",i), Get<float>("T_Muon_Energy",i)));

    }
  }
 
  GetAllMuons(); 

//------------------------------------------------------------------------------
// Get all GEN prompt muon
//------------------------------------------------------------------------------
  
  SetGenInfo();

//------------------------------------------------------------------------------
// Get number of good vertex per event
//------------------------------------------------------------------------------

  unsigned int N_PV = 0;
  N_PV = Get<std::vector<float>*>("T_Vertex_z")->size();//GetSizeOf("T_Vertex_z");
  h_N_PV->Fill(N_PV, factN);
    
//------------------------------------------------------------------------------
// DILEPTONIC SELECTION
//------------------------------------------------------------------------------
  
  // Select the first two muons with higher pt (i==0 and i==1) 
  
  if ( Get<std::vector<float>*>("T_Muon_Px")->size() > 1 && (G_isMuMu || G_isMuTau || G_isTauMu || G_isTauTau)) { 
    
    if ( !_Signal.Contains("ISO")                                                            &&
	 //G_GEN_PromptMuon_4vec[0].Pt() > 10 &&  fabs(G_GEN_PromptMuon_4vec[0].Eta()) < 2.4   &&
	 //( ( !_Signal.Contains("Wjets") && G_GEN_PromptMuon_4vec[1].Pt() > 10                &&
	 //    fabs(G_GEN_PromptMuon_4vec[1].Eta()) < 2.4 ) || _Signal.Contains("Wjets") )     &&
	 G_MuonID_Fiducial[0]        && G_MuonID_Fiducial[1]                                 &&
	 G_MuonID_GLBorTRKArb[0]     && G_MuonID_GLBorTRKArb[1]                              ) {
 
      int ch1 = Get<int>("T_Muon_Charge",0);
      int ch2 = Get<int>("T_Muon_Charge",1);

      if ( ch1*ch2 < 0 ) {

	h_N_PV2->Fill(N_PV, factN);
	
	int Mu1 = 0; 
	int Mu2 = 1;

	bool isMatchGenReco = false;
	isMatchGenReco = MatchGenToReco(Mu1, Mu2);

	if (isMatchGenReco) { 

	  if (G_Debug_DefineAnalysisVariables) std::cout << "[DEBUG][Event "<< 
						 Get<int>("T_Event_EventNumber") <<"] Matched!" << std::endl;

	  h_N_PV3->Fill(N_PV, factN);

	}
     
      }
      
      
    }
  }
  
} // end inside Loop

////////////////////////////////////////////////////////////////////////////////
// MEMBER FUNCTIONS
//

//------------------------------------------------------------------------------
// Get All RECO muons (ID + ISO)
//------------------------------------------------------------------------------

void BasicmuonAnalyzer::GetAllMuons() {

 
  UInt_t _muonSize = 0;
  _muonSize = Get<std::vector<float>*>("T_Muon_Px")->size();

  int count = 0;
  
  if ( _muonSize > 0 ) {  // asking for at least one muon in the event 

    for (unsigned int i = 0; i < _muonSize; i++) {
      
      //  Define selection for numerator and denominator, and general Muon ID

      const int NFLAGS = 17;
      bool muon_sel[NFLAGS];   
      for (int j=0; j<NFLAGS; ++j) muon_sel[j] = false;

      muon_sel[0] = G_Muon_4vec[i].Pt() > 10. && fabs(G_Muon_4vec[i].Eta()) < 2.4;
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
      
      muon_sel[15] = (((muon_sel[2] && muon_sel[4] && muon_sel[5] && muon_sel[6]) || muon_sel[3]) &&
		     muon_sel[7] && muon_sel[8] && muon_sel[11] && muon_sel[12]);
      muon_sel[16] = (muon_sel[2] || muon_sel[3]);

      
      G_MuonID_Tight.push_back(muon_sel[13]);
      G_MuonID_Medium.push_back(muon_sel[14]);
      G_MuonID_HWW.push_back(muon_sel[15]);
      G_MuonID_MediumHWW.push_back(muon_sel[14] * muon_sel[11] * muon_sel[12]);
      G_MuonID_GLBorTRKArb.push_back(muon_sel[16]);
      G_MuonID_Fiducial.push_back(muon_sel[0]);


      G_MuonISO03.push_back(          passPFIso(i, "R03",          0.12));
      G_MuonISO03_dBeta.push_back(    passPFIso(i, "dBetaR03",     0.12));
      G_MuonISO03_PFWeights.push_back(passPFIso(i, "PFWeightsR03", 0.12));
      G_MuonISO04.push_back(          passPFIso(i, "R04",          0.12));
      G_MuonISO04_dBeta.push_back(    passPFIso(i, "dBetaR04",     0.12));
      G_MuonISO04_PFWeights.push_back(passPFIso(i, "PFWeightsR04", 0.12));

      if (muon_sel[0]) count++;

      
    } // end loop on muons

    if (G_Debug_DefineAnalysisVariables) std::cout << "[DEBUG][Event "<< Get<int>("T_Event_EventNumber") 
						   <<"] Got " << _muonSize <<" RECO muon(s)" << std::endl;

    if (G_Debug_DefineAnalysisVariables) std::cout << "[DEBUG][Event "<< Get<int>("T_Event_EventNumber") 
						   <<"] Got " << count <<" RECO fiducial muon(s)" << std::endl;

  } // end loop on at least one muon
  
}


bool BasicmuonAnalyzer::passMediumID(int iMu) {

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


//------------------------------------------------------------------------------
// PF isolation 
//-----------------------------------------------------------------------------

bool BasicmuonAnalyzer::passPFIso (int iMu, string typeIso, float wp) {

  bool passIso = false;

  float PFRelIsoBeta = getPFRelIso(iMu, typeIso);
     
  if ( PFRelIsoBeta <=  wp)   passIso = true;	  

  return passIso;

}

float BasicmuonAnalyzer::getPFRelIso (int iMu, string typeIso) {
  
  float PFRelIso = 999.9;
  float pt = Get<float>("T_Muon_Pt",iMu);

  if (typeIso == "R03")
    PFRelIso = ( Get<float>("T_Muon_chargedHadronIsoR03",iMu) + 
		 Get<float>("T_Muon_neutralHadronIsoR03",iMu) + 
		 Get<float>("T_Muon_photonIsoR03",iMu) )
      / pt;

  else if (typeIso == "R04")
    PFRelIso = ( Get<float>("T_Muon_chargedHadronIsoR04",iMu) + 
		 Get<float>("T_Muon_neutralHadronIsoR04",iMu) + 
		 Get<float>("T_Muon_photonIsoR04",iMu) )
      / pt;

  else if (typeIso == "dBetaR03") 
    PFRelIso = ( Get<float>("T_Muon_chargedHadronIsoR03",iMu) + max(0., Get<float>("T_Muon_neutralHadronIsoR03",iMu) + 
								    Get<float>("T_Muon_photonIsoR03",iMu) - 
								    0.5 * Get<float>("T_Muon_sumPUPtR03",iMu)) )
      / pt;

  else if (typeIso == "dBetaR04") 
    PFRelIso = ( Get<float>("T_Muon_chargedHadronIsoR04",iMu) + max(0., Get<float>("T_Muon_neutralHadronIsoR04",iMu) + 
								    Get<float>("T_Muon_photonIsoR04",iMu) - 
								    0.5 * Get<float>("T_Muon_sumPUPtR04",iMu)) )
      / pt;

  else if (typeIso == "PFWeightsR03")
    PFRelIso = ( Get<float>("T_Muon_chargedHadronIsoR03",iMu) + Get<float>("T_Muon_neutralIsoPFweightR03",iMu) )
      / pt;

  else if (typeIso == "PFWeightsR04")
    PFRelIso = ( Get<float>("T_Muon_chargedHadronIsoR04",iMu) + Get<float>("T_Muon_neutralIsoPFweightR04",iMu) )
      / pt;

  return PFRelIso;

}


// ---------------------------------------------------
// Set Generator level info
// ---------------------------------------------------

void BasicmuonAnalyzer::SetGenInfo() {
  
  UInt_t genPromptMuSize = 0;
  genPromptMuSize = Get<std::vector<float>*>("T_Gen_PromptMuon_Px")->size();
  
  UInt_t genPromptTauSize = 0;
  genPromptTauSize = Get<std::vector<float>*>("T_Gen_PromptTau_Px")->size();
  
  TLorentzVector p1 = TLorentzVector(0,0,0,0);
  TLorentzVector p2 = TLorentzVector(0,0,0,0);

  // if (G_Debug_DefineAnalysisVariables) std::cout << "[DEBUG][Event "<< Get<int>("T_Event_EventNumber") <<
  // 					 "] genPromptMuSize: " << genPromptMuSize << ", genPromptTauSize: " <<
  // 					 genPromptTauSize << std::endl;

  if ( _Signal.Contains("Wjets"))
  {
     
    UInt_t genNonPromptMuSize = 0;
    genNonPromptMuSize = Get<std::vector<float>*>("T_Gen_Muon_Px")->size();
    
    if ( genPromptMuSize == 1 && fabs(Get<int>("T_Gen_PromptMuon_MpdgId",0)) == 24) 
      
      G_isMuMu = true; 
    
    if ( genPromptMuSize < 1 && genPromptTauSize == 1 && 
  	 fabs(Get<int>("T_Gen_PromptTau_MpdgId",0))== 24 && 
  	 fabs(Get<int>("T_Gen_PromptTau_LepDec_pdgId",0)) == 13) 

      G_isTauMu = true; 
    
    
    if ( G_isMuMu ) {
      G_GEN_PromptMuon_4vec.push_back(TLorentzVector(Get<float>("T_Gen_PromptMuon_Px",0), 
  						     Get<float>("T_Gen_PromptMuon_Py",0),
  						     Get<float>("T_Gen_PromptMuon_Pz",0), 
  						     Get<float>("T_Gen_PromptMuon_Energy",0)));
      if ( genNonPromptMuSize > 0) {
  	G_isNonPrompt = true;
  	G_GEN_Muon_4vec.push_back(TLorentzVector(Get<float>("T_Gen_Muon_Px",0), 
  						 Get<float>("T_Gen_Muon_Py",0),
  						 Get<float>("T_Gen_Muon_Pz",0), 
  						 Get<float>("T_Gen_Muon_Energy",0)));
      }
    } 
    
    if ( G_isTauMu ) {
      G_GEN_PromptMuon_4vec.push_back(TLorentzVector(Get<float>("T_Gen_PromptTau_LepDec_Px",0),
  						     Get<float>("T_Gen_PromptTau_LepDec_Py",0),
  						     Get<float>("T_Gen_PromptTau_LepDec_Pz",0), 
  						     Get<float>("T_Gen_PromptTau_LepDec_Energy",0)));
      if ( genNonPromptMuSize > 1) {
  	G_isNonPrompt = true;
  	G_GEN_Muon_4vec.push_back(TLorentzVector(Get<float>("T_Gen_Muon_Px",1), 
  						 Get<float>("T_Gen_Muon_Py",1),
  						 Get<float>("T_Gen_Muon_Pz",1), 
  						 Get<float>("T_Gen_Muon_Energy",1)));
      }
    }
    
  }
    
  if (_Signal.Contains("GGHWW") || _Signal.Contains("TTbar"))
    {

      //if (Get<std::vector<float>*>("T_Gen_W_pt")->size() != 2) return;
      //if (!( _Signal.Contains("TTbar") 
      //	     && fabs(Get<int>("T_Gen_W_MpdgId",0))==6 && fabs(Get<int>("T_Gen_W_MpdgId",1))==6 )) return;
      
      if ( genPromptMuSize == 2 && fabs(Get<int>("T_Gen_PromptMuon_MpdgId",0)) == 24 && 
  	   fabs(Get<int>("T_Gen_PromptMuon_MpdgId",1)) == 24 &&
  	   (Get<int>("T_Gen_PromptMuon_pdgId",0)*Get<int>("T_Gen_PromptMuon_pdgId",1)) < 0) {

  	p1 = TLorentzVector(Get<float>("T_Gen_PromptMuon_Px",0), 
  			    Get<float>("T_Gen_PromptMuon_Py",0),
  			    Get<float>("T_Gen_PromptMuon_Pz",0), 
  			    Get<float>("T_Gen_PromptMuon_Energy",0));

  	p2 = TLorentzVector(Get<float>("T_Gen_PromptMuon_Px",1), 
  			    Get<float>("T_Gen_PromptMuon_Py",1),
  			    Get<float>("T_Gen_PromptMuon_Pz",1), 
  			    Get<float>("T_Gen_PromptMuon_Energy",1));

  	G_isMuMu = true; 

      }
      
      if ( genPromptMuSize == 1 && fabs(Get<int>("T_Gen_PromptMuon_MpdgId",0)) == 24 && 
  	   genPromptTauSize == 1 && fabs(Get<int>("T_Gen_PromptTau_MpdgId",0))== 24 && 
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

  	if (p1.Pt() >= p2.Pt()) G_isMuTau = true;
  	else                    G_isTauMu = true;

      }
      
      if ( genPromptMuSize < 1 && genPromptTauSize == 2 && 
  	   fabs(Get<int>("T_Gen_PromptTau_MpdgId",0))== 24 && 
  	   fabs(Get<int>("T_Gen_PromptTau_MpdgId",1))== 24 && 
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

  	G_isTauTau = true; 

      }
      
    }
  
  if (_Signal.Contains("DY"))
    {

      if ( genPromptMuSize == 2 && fabs(Get<int>("T_Gen_PromptMuon_MpdgId",0)) == 23 
      	   && fabs(Get<int>("T_Gen_PromptMuon_MpdgId",1)) == 23 &&
      	   (Get<int>("T_Gen_PromptMuon_pdgId",0)*Get<int>("T_Gen_PromptMuon_pdgId",1)) < 0) {

  	p1 = TLorentzVector(Get<float>("T_Gen_PromptMuon_Px",0), 
  			    Get<float>("T_Gen_PromptMuon_Py",0),
  			    Get<float>("T_Gen_PromptMuon_Pz",0), 
  			    Get<float>("T_Gen_PromptMuon_Energy",0));

  	p2 = TLorentzVector(Get<float>("T_Gen_PromptMuon_Px",1), 
  			    Get<float>("T_Gen_PromptMuon_Py",1),
  			    Get<float>("T_Gen_PromptMuon_Pz",1), 
  			    Get<float>("T_Gen_PromptMuon_Energy",1));

  	G_isMuMu = true; 

      }
      
      if ( genPromptMuSize == 1 && fabs(Get<int>("T_Gen_PromptMuon_MpdgId",0)) == 23 
      	   && genPromptTauSize == 1 && fabs(Get<int>("T_Gen_PromptTau_MpdgId",0))== 23 
      	   && fabs(Get<int>("T_Gen_PromptTau_LepDec_pdgId",0)) == 13 &&
      	   (Get<int>("T_Gen_PromptMuon_pdgId",0)*Get<int>("T_Gen_PromptTau_LepDec_pdgId",0)) < 0) {

  	p1 = TLorentzVector(Get<float>("T_Gen_PromptMuon_Px",0), 
  			    Get<float>("T_Gen_PromptMuon_Py",0),
  			    Get<float>("T_Gen_PromptMuon_Pz",0), 
  			    Get<float>("T_Gen_PromptMuon_Energy",0));

  	p2 = TLorentzVector(Get<float>("T_Gen_PromptTau_LepDec_Px",0), 
  			    Get<float>("T_Gen_PromptTau_LepDec_Py",0),
  			    Get<float>("T_Gen_PromptTau_LepDec_Pz",0), 
  			    Get<float>("T_Gen_PromptTau_LepDec_Energy",0));

  	if (p1.Pt() >= p2.Pt()) G_isMuTau = true;
  	else                    G_isTauMu = true;

      }
      
      if ( genPromptMuSize < 1 && genPromptTauSize == 2 && 
      	   fabs(Get<int>("T_Gen_PromptTau_MpdgId",0))== 23 && 
      	   fabs(Get<int>("T_Gen_PromptTau_MpdgId",1))== 23 && 
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

  	G_isTauTau = true; 

      }
      
    }

  if ((_Signal.Contains("DY") || _Signal.Contains("GGHWW") || _Signal.Contains("TTbar")) && 
      (G_isMuMu || G_isMuTau || G_isTauMu || G_isTauTau)) {

    if ( p1.Pt() >= p2.Pt() ) {

      G_GEN_PromptMuon_4vec.push_back(p1);
      G_GEN_PromptMuon_4vec.push_back(p2);

    }

    else {

      G_GEN_PromptMuon_4vec.push_back(p2);
      G_GEN_PromptMuon_4vec.push_back(p1);

    }

    if (G_Debug_DefineAnalysisVariables) std::cout << "[DEBUG][Event "<< 
						 Get<int>("T_Event_EventNumber") <<"] Got 2 GEN muons" << std::endl;

  }
        

}


//-----------------------------------------------------------
// Match Muons to Gen bosons
//-----------------------------------------------------------

bool BasicmuonAnalyzer::MatchGenToReco(int &mu1, int &mu2) {
  
  bool isMatched = false; 

  if (_Signal.Contains("GGHWW") || _Signal.Contains("DY") || _Signal.Contains("TTbar")) {

    mu1 = 0; mu2 = 1;

    float dR11 = 999.0;
    float dR12 = 999.0;
    float dR21 = 999.0;
    float dR22 = 999.0;
    
    if (G_isMuMu || G_isMuTau || G_isTauMu || G_isTauTau) { 

      //isMatched = true;

      dR11 =  G_Muon_4vec[0].DeltaR(G_GEN_PromptMuon_4vec[0]);
      dR22 =  G_Muon_4vec[1].DeltaR(G_GEN_PromptMuon_4vec[1]);

      dR12 =  G_Muon_4vec[0].DeltaR(G_GEN_PromptMuon_4vec[1]);
      dR21 =  G_Muon_4vec[1].DeltaR(G_GEN_PromptMuon_4vec[0]);

      //if (dR11 < 0.1) isMatched = true;

      if (dR11 < 0.1 && dR22 < 0.1 && dR11 < dR12 && dR22 < dR21) {
	isMatched = true;
	
	// for (UInt_t i = 0; i < GetSizeOf("T_Gen_Promptb_Px"); i++) {
	  
	//   TLorentzVector b(Get<float>("T_Gen_Promptb_Px",i), 
	// 		   Get<float>("T_Gen_Promptb_Py",i),
	// 		   Get<float>("T_Gen_Promptb_Pz",i), 
	// 		   Get<float>("T_Gen_Promptb_Energy",i));
	  
	//   for (UInt_t k=0; k<2; k++) {

	//     float dRbMu = b.DeltaR(G_Muon_4vec[k]);
	//     h_Dilep_dRbMu[k]->Fill(dRbMu);

	//     if ( dRbMu < 1.0) isMatched = false;

	//     if (G_Muon_4vec[k].Pt() < 20) continue;
	//     if (G_Muon_TightID[k])          h_Dilep_dRbMu_TightID[k]->Fill(dRbMu);
	//     if (G_Muon_TightISOR04[k])      h_Dilep_dRbMu_TightISOR04[k]->Fill(dRbMu);
	//     if (G_Muon_TightISOdBetaR04[k]) h_Dilep_dRbMu_TightISOdBetaR04[k]->Fill(dRbMu);

	//     // for (UInt_t j=0; j<T_JetAKCHS_Px->size(); j++) {
	//     //   if ( G_Jet_4vec[j].Pt()  < 30)  continue;
	//     //   if ( G_Jet_4vec[j].Eta() > 2.4) continue;
	//     //   if (!passJetID(j))              continue;
	//     //   if (G_Jet_4vec[j].DeltaR(G_Muon_4vec[k]) < 0.4) isMatched = false;
	//     // }

	//   }
	  
	// }
      }

      //else if (dR12 < 0.3 && dR21 < 0.3 && dR12 < dR11 && dR21 < dR22) { mu1 = 1; mu2 = 0; isMatched = true; }

    }

  }

  else if ( _Signal.Contains("Wjets")) {

    float dR1 = 999.0;
    float dR2 = 999.0;
	
    if ( (G_isMuMu || G_isTauMu)) {

      isMatched = false;
      
      if (G_isMuMu || G_isTauMu) {
	dR1 =  G_Muon_4vec[0].DeltaR(G_GEN_PromptMuon_4vec[0]);
	dR2 =  G_Muon_4vec[1].DeltaR(G_GEN_PromptMuon_4vec[0]);
      }
      
      if ( dR1 < 0.1 && dR1 < dR2 ) { 

	mu1 = 0; mu2 = 1; isMatched = true;

      }   
    
      // float dRfake = 999.0;
      // int isFakeMatched = 0;
      
      // if (G_isNonPrompt) {
      // 	dRfake = G_GEN_Muon_4vec[0].DeltaR(G_Muon_4vec[mu2]);
      // 	if (dRfake < 0.4) isFakeMatched = 1;	  
      // }
      
      // if ( isFakeMatched ) h_Dilep_Gen_Muon_MpdgId->Fill(Get<int>("T_Gen_Muon_MpdgId",0));
      
    }

  }
    
  return isMatched;

}


void BasicmuonAnalyzer::Summary() {
  // Get Data Members at the client-master (after finishing the analysis at the workers nodes)
  // Only data members set here will be accesible at the client-master


  ///*** 1D histos ***/// 

  h_N_PV  = FindOutput<TH1F*>("h_N_PV");
  h_N_PV2 = FindOutput<TH1F*>("h_N_PV2");
  h_N_PV3 = FindOutput<TH1F*>("h_N_PV3");
  

  cout << " ---------------------------------------------------" << endl;
  cout << " " << endl;  

  cout << " Number of Events::  " << _NEvents  << endl;

  float factN = 1.; 
  if (_XSection > 0) factN = _XSection * _Luminosity / _NEvents; //fractionoftotalevents;

  cout << " Normalization factor: " << factN << endl;
  cout << endl;
  cout << " ---------------------------------------------------" << endl;

}
