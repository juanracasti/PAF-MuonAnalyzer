///////////////////////////////////////////////////////////////////////
//
//
/////////////////////////////////////////////////////////////
/////////////          MUON ANALYSIS            /////////////
/////////////////////////////////////////////////////////////

#include "TSelector.h"

void RunBasicmuonAnalyzer(const char* data) {
 
  gSystem->Load("libPAF.so");   

  float luminosity = 20000.0;
  TString dataPath="/gpfs/csic_projects/tier3data/TreesPHYS14/PU20bx25/";
  TString signal = data;
  bool debug =  false; //For verbose while debugging

  bool isdata = false;
  int whichRun = 1;
  int nEventsInTheSample = 1;
  float xSection = 1.;
  

  // **** INITIALISE PAF PROJECT

  PAFProject* myProject = new PAFProject( new PAFSequentialEnvironment() );

  myProject->SetDefaultTreeName("demo/Tree");



  //+++++++++++++++++++++++++++++++++
  // ADD SAMPLES

  if (signal=="MC_GGHWW_PU20bx25") {

    myProject->AddDataFile(dataPath + "Tree_HWW125.root");
    
    isdata = false;
    nEventsInTheSample = 99555; 
    xSection =  1.0 ;
    whichRun = 2; 
 
  }

  else if (signal=="MC_Wjets_PU20bx25") {

    myProject->AddDataFile(dataPath + "Tree_WJets_Madgraph.root");
      
    isdata = false;
    nEventsInTheSample = 10017462; 
    xSection =  20508.9;
    whichRun = 2;

 }

  else if (signal=="MC_DY_PU20bx25") {

    myProject->AddDataFile(dataPath + "Tree_ZJets_Madgraph.root");
   
    isdata = false;
    nEventsInTheSample = 2829164; 
    xSection =  6025.2;
    whichRun = 2;

 }

  else if (signal=="MC_TTbar_PU20bx25") {

    myProject->AddDataFile(dataPath + "Tree_TTJets_MadSpin_0.root");
    myProject->AddDataFile(dataPath + "Tree_TTJets_MadSpin_1.root");
    myProject->AddDataFile(dataPath + "Tree_TTJets_MadSpin_2.root");
    myProject->AddDataFile(dataPath + "Tree_TTJets_MadSpin_3.root");
    myProject->AddDataFile(dataPath + "Tree_TTJets_MadSpin_4.root");
    myProject->AddDataFile(dataPath + "Tree_TTJets_MadSpin_5.root");
    myProject->AddDataFile(dataPath + "Tree_TTJets_MadSpin_6.root");
    myProject->AddDataFile(dataPath + "Tree_TTJets_MadSpin_7.root");
   
    isdata = false;
    nEventsInTheSample = 2829164; 
    xSection =  6025.2;
    whichRun = 2;

 }


  else if (signal=="DY_ISO_PU20bx25") {

    myProject->AddDataFile(dataPath + "Tree_ZJets_Madgraph.root");
   
    isdata = false;
    nEventsInTheSample = 2829164; 
    xSection =  6025.2;
    whichRun = 2;

 }

  else if (signal=="TTbar_ISO_PU20bx25") {

    myProject->AddDataFile(dataPath + "Tree_TTJets_MadSpin_0.root");
    myProject->AddDataFile(dataPath + "Tree_TTJets_MadSpin_1.root");
    myProject->AddDataFile(dataPath + "Tree_TTJets_MadSpin_2.root");
    myProject->AddDataFile(dataPath + "Tree_TTJets_MadSpin_3.root");
    myProject->AddDataFile(dataPath + "Tree_TTJets_MadSpin_4.root");
    myProject->AddDataFile(dataPath + "Tree_TTJets_MadSpin_5.root");
    myProject->AddDataFile(dataPath + "Tree_TTJets_MadSpin_6.root");
    myProject->AddDataFile(dataPath + "Tree_TTJets_MadSpin_7.root");
   
    isdata = false;
    nEventsInTheSample = 2829164; 
    xSection =  6025.2;
    whichRun = 2;

 }

  else if (signal=="MC_TTbar_PU30bx50") {

    myProject->AddDataFile("/gpfs/csic_projects/tier3data/TreesPHYS14/PU30bx50/Tree_TTJets_MadSpin_0.root");
    myProject->AddDataFile("/gpfs/csic_projects/tier3data/TreesPHYS14/PU30bx50/Tree_TTJets_MadSpin_1.root");
    myProject->AddDataFile("/gpfs/csic_projects/tier3data/TreesPHYS14/PU30bx50/Tree_TTJets_MadSpin_2.root");
   
    isdata = false;
    nEventsInTheSample = 2829164; 
    xSection =  6025.2;
    whichRun = 2;

 }

  else if (signal=="QCD_ISO_PU20bx25") {

    myProject->AddDataFile(dataPath + "Tree_QCD_MuEnrichedPt15.root");
   
    isdata = false;
    nEventsInTheSample = 787547; 
    xSection =  1;
    whichRun = 2;

 }

  else if (signal=="SingleMu_720") {

    myProject->AddDataFile(dataPath + "Tree_SingleMu_mu2012D_720.root");
   
    isdata = true;
    nEventsInTheSample = 301812;  
    xSection =  1;
    whichRun = 1;

 }

  else if (signal=="SingleMu_740") {

    myProject->AddDataFile(dataPath + "Tree_SingleMu_mu2012D_740.root");
   
    isdata = true;
    nEventsInTheSample = 299035; 
    xSection =  1;
    whichRun = 1;

 }

  ///////////////////////////////
  // INPUT PARAMETERS
 
  myProject->SetInputParam("IsDATA", isdata);
  myProject->SetInputParam("Signal", data);
  myProject->SetInputParam("XSection", xSection);
  myProject->SetInputParam("Luminosity", luminosity);
  myProject->SetInputParam("NEvents", nEventsInTheSample); //before skimming
  myProject->SetInputParam("luminosityPU", 19468.3);  
  myProject->SetInputParam("WhichRun", whichRun);
  myProject->SetInputParam("Debug", debug);

  ///////////////////////////////
  // OUTPUT FILE NAME
  // Specify the name of the file where you want your histograms to be saved

  myProject->SetOutputFile("phys14_"+signal+".root"); 

  ///////////////////////////////
  // SELECTOR AND PACKAGES

  myProject->AddSelectorPackage("BasicmuonAnalyzer");

  // RUN THE ANALYSIS
  // ================

  myProject->Run();

  //
  /////////////////////////////////////////////////////////////////////////


}
