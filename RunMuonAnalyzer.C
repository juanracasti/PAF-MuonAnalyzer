///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////                                                                                             /////////////
/////////////                                    MUON ANALYSIS WITH PAF                                   /////////////
/////////////                                                                                             /////////////
/////////////                                  Juan R. Castiñeiras (IFCA)                                 /////////////
/////////////                                          Jun 2016                                           /////////////
/////////////                                                                                             /////////////
/////////////                              -> Adjust to a 120 width window <-                             /////////////
/////////////                                                                                             /////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "TSelector.h"

void RunMuonAnalyzer(const char* data, const TString environment) {
  // Function input parameters:
  // + data :        String to select the Trees to read
  //                 Possible values:
  //                     > "DR74X_50ns_MC_DY"
  //                     > "DR74X_50ns_MC_TTbar"
  //                     > "DR74X_50ns_MC_QCD"
  //                     > "PHYS14_PU20bx25_MC_GGHWW"
  //                     > "PHYS14_PU20bx25_MC_Wjets"
  //                     > "PHYS14_PU20bx25_MC_DY"
  //                     > "PHYS14_PU20bx25_MC_TTbar"
  //                     > "PHYS14_PU30bx50_MC_TTbar"
  //                     > "SingleMu_720"
  //                     > "SingleMu_740"
  // + environment : String to select the desired PAF Execution Environment
  //                 Possible values:
  //                     > "Sequential"
  //                     > "PROOFLite"
  //                     > "PROOFCluster"
 
  gSystem->Load("libPAF.so");   

  TString dataPath="/gpfs/csic_projects/tier3data/";
  TString Path_PHYS14 = "/gpfs/csic_projects/tier3data/TreesPHYS14/";
  TString Path_DR74X = "/gpfs/csic_projects/tier3data/TreesDR74X/";
  TString signal = data;

  // Manual Input Parameters
  float    luminosity       = 20000.0; //In pb-1
  bool     debug            = false;   //For verbose while debugging
  int      nEventsToProcess = -1;   //Number of events to be processed (-1 = all)
  bool     doReport         = false;   //Count events and print final report

  // Automatic Input Parameters (don't touch)
  bool     isdata             = false;
  int      whichRun           = 1;
  int      nEventsInTheSample = 1; //before skimming
  float    xSection           = 1.;
  
  //---------------------------------
  // INITIALISE PAF PROJECT
  //---------------------------------

  // Create Project
  PAFProject* myProject;

  // Initialize Project and set PAF Execution Environment
  if (environment == "Sequential")        myProject = new PAFProject( new PAFSequentialEnvironment() );
  else if (environment == "PROOFLite")    myProject = new PAFProject( new PAFPROOFLiteEnvironment(6) );
  else if (environment == "PROOFCluster") myProject = new PAFProject( new PAFPROOFClusterEnvironment(10,99) );
  else {
    cout << "[WARNING!] The selected PAF Execution Environment doesn't exist" << endl;
    cout << "[WARNING!] Switching to Sequential enviroment" << endl;
    myProject = new PAFProject( new PAFSequentialEnvironment() );
  }


  // Set default name and subdirectory of Trees
  myProject->SetDefaultTreeName("demo/Tree");

  ////////////////////////////////////////////////
  // ADD SAMPLES

  if (signal=="DR74X_50ns_MC_DY") {

    myProject->AddDataFile(Path_DR74X + "50ns/PuppiVar/Tree_ZJets_aMCatNLO_0.root");
    myProject->AddDataFile(Path_DR74X + "50ns/PuppiVar/Tree_ZJets_aMCatNLO_1.root");
    myProject->AddDataFile(Path_DR74X + "50ns/PuppiVar/Tree_ZJets_aMCatNLO_2.root");
    myProject->AddDataFile(Path_DR74X + "50ns/PuppiVar/Tree_ZJets_aMCatNLO_3.root");
    myProject->AddDataFile(Path_DR74X + "50ns/PuppiVar/Tree_ZJets_aMCatNLO_4.root");
   
    isdata             = false;
    nEventsInTheSample = 6810943; 
    xSection           = 6025.2;
    whichRun           = 2;

 }

  else if (signal=="DR74X_50ns_MC_TTbar") {

    myProject->AddDataFile(Path_DR74X + "50ns/PuppiVar/Tree_TTbar_Powheg_0.root");
   
    isdata             = false;
    nEventsInTheSample = 1079068; 
    xSection           = 6025.2;
    whichRun           = 2;

 }

  else if (signal=="DR74X_50ns_MC_QCD") {

    myProject->AddDataFile(Path_DR74X + "50ns/PuppiVar/Tree_QCDmuEnriched_0.root");
    myProject->AddDataFile(Path_DR74X + "50ns/PuppiVar/Tree_QCDmuEnriched_1.root");
   
    isdata             = false;
    nEventsInTheSample = 1848472; 
    xSection           = 1.;
    whichRun           = 2;

 }

  else if (signal=="PHYS14_PU20bx25_MC_GGHWW") {

    myProject->AddDataFile(Path_PHYS14 + "PU20bx25/Tree_HWW125.root");
    
    isdata             = false;
    nEventsInTheSample = 99555; 
    xSection           = 1.;
    whichRun           = 2; 
 
  }

  else if (signal=="PHYS14_PU20bx25_MC_Wjets") {

    myProject->AddDataFile(Path_PHYS14 + "PU20bx25/Tree_WJets_Madgraph.root");
      
    isdata             = false;
    nEventsInTheSample = 10017462; 
    xSection           = 20508.9;
    whichRun           = 2;

 }

  else if (signal=="PHYS14_PU20bx25_MC_DY") {

    myProject->AddDataFile(Path_PHYS14 + "PU20bx25/Tree_ZJets_Madgraph.root");
   
    isdata             = false;
    nEventsInTheSample = 2829164; 
    xSection           = 6025.2;
    whichRun           = 2;

 }

  else if (signal=="PHYS14_PU20bx25_MC_TTbar") {

    myProject->AddDataFile(Path_PHYS14 + "PU20bx25/Tree_TTJets_MadSpin_0.root");
    myProject->AddDataFile(Path_PHYS14 + "PU20bx25/Tree_TTJets_MadSpin_1.root");
    myProject->AddDataFile(Path_PHYS14 + "PU20bx25/Tree_TTJets_MadSpin_2.root");
    myProject->AddDataFile(Path_PHYS14 + "PU20bx25/Tree_TTJets_MadSpin_3.root");
    myProject->AddDataFile(Path_PHYS14 + "PU20bx25/Tree_TTJets_MadSpin_4.root");
    myProject->AddDataFile(Path_PHYS14 + "PU20bx25/Tree_TTJets_MadSpin_5.root");
    myProject->AddDataFile(Path_PHYS14 + "PU20bx25/Tree_TTJets_MadSpin_6.root");
    myProject->AddDataFile(Path_PHYS14 + "PU20bx25/Tree_TTJets_MadSpin_7.root");
   
    isdata             = false;
    nEventsInTheSample = 2829164; 
    xSection           = 6025.2;
    whichRun           = 2;

 }

  else if (signal=="PHYS14_PU30bx50_MC_TTbar") {

    myProject->AddDataFile(Path_PHYS14 + "PU30bx50/Tree_TTJets_MadSpin_0.root");
    myProject->AddDataFile(Path_PHYS14 + "PU30bx50/Tree_TTJets_MadSpin_1.root");
    myProject->AddDataFile(Path_PHYS14 + "PU30bx50/Tree_TTJets_MadSpin_2.root");
   
    isdata             = false;
    nEventsInTheSample = 2829164; 
    xSection           = 6025.2;
    whichRun           = 2;

 }

  else if (signal=="PHYS14_PU20bx25_MC_QCD") {

    myProject->AddDataFile(Path_PHYS14 + "PU20bx25/Tree_QCD_MuEnrichedPt15.root");
   
    isdata             = false;
    nEventsInTheSample = 787547; 
    xSection           = 1;
    whichRun           = 2;

 }

  else if (signal=="SingleMu_720") {

    myProject->AddDataFile(Path_PHYS14 + "PU20bx25/Tree_SingleMu_mu2012D_720.root");
   
    isdata             = true;
    nEventsInTheSample = 301812;  
    xSection           = 1;
    whichRun           = 1;

 }

  else if (signal=="SingleMu_740") {

    myProject->AddDataFile(Path_PHYS14 + "PU20bx25/Tree_SingleMu_mu2012D_740.root");
   
    isdata             = true;
    nEventsInTheSample = 299035; 
    xSection           = 1;
    whichRun           = 1;

 }

  //Number of events to process
  myProject->SetNEvents(nEventsToProcess);

  ///////////////////////////////
  // INPUT PARAMETERS
 
  myProject->SetInputParam("IsDATA",       isdata);
  myProject->SetInputParam("Signal",       data);
  myProject->SetInputParam("XSection",     xSection);
  myProject->SetInputParam("Luminosity",   luminosity);
  myProject->SetInputParam("NEvents",      nEventsInTheSample);
  myProject->SetInputParam("luminosityPU", 19468.3);  
  myProject->SetInputParam("WhichRun",     whichRun);
  myProject->SetInputParam("Debug",        debug);
  myProject->SetInputParam("Report",       doReport);

  ///////////////////////////////
  // OUTPUT FILE NAME
  // Specify the name of the file where you want your histograms to be saved

  myProject->SetOutputFile("files/"+signal+".root"); 

  ///////////////////////////////
  // SELECTOR AND PACKAGES

  //To Run in one step (Please, use this one for now!)
  myProject->AddSelectorPackage("MuonMasterSelector");

  //To Run by modules (CURRENTLY UNSUPPORTED, DO NOT USE!)
  //First add the core selector, this one is essential 
  // myProject->AddSelectorPackage("MuonCoreSelector");

  //Then, add the additional selectors if needed (and if available!)
  // myProject->AddSelectorPackage("MuonEfficienciesSelector");
  // myProject->AddSelectorPackage("MuonISORocCurvesSelector");

  // RUN THE ANALYSIS
  // ================

  myProject->Run();

  
  /////////////////////////////////////////////////////////////////////////


}
